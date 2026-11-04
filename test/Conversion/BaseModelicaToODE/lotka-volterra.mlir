// Lotka-Volterra from Base Modelica through the whole flow: ode.system,
// composable rules, network, emulator JSON, simulation check.

// RUN: chemcomp-opt %s --convert-bmodelica-to-ode > %t.ode.mlir
// RUN: chemcomp-opt %t.ode.mlir --ode-to-crn | FileCheck %s
// RUN: chemcomp-opt %t.ode.mlir --ode-to-crn | chemcomp-translate -mlir-to-json \
// RUN:   > %t.crn.json
// RUN: FileCheck %s --check-prefix=JSON < %t.crn.json
// RUN: chemcomp-translate -ode-to-json %t.ode.mlir > %t.ode.json
// RUN: %if emulator %{ %chemcomp_python %check_realization --ode %t.ode.json \
// RUN:   --crn %t.crn.json --t-end 30 | FileCheck %s --check-prefix=SIM %}

!r = !bmodelica.real

bmodelica.model @LotkaVolterra {
    bmodelica.variable @rabbit : !bmodelica.variable<!r>
    bmodelica.variable @fox : !bmodelica.variable<!r>
    bmodelica.variable @alpha : !bmodelica.variable<!r, parameter>
    bmodelica.variable @beta : !bmodelica.variable<!r, parameter>
    bmodelica.variable @gamma : !bmodelica.variable<!r, parameter>

    bmodelica.start @rabbit {
        %0 = bmodelica.constant #bmodelica<real 2.0>
        bmodelica.yield %0 : !r
    } {each = false, fixed = true}
    bmodelica.start @fox {
        %0 = bmodelica.constant #bmodelica<real 1.0>
        bmodelica.yield %0 : !r
    } {each = false, fixed = true}
    bmodelica.binding_equation @alpha {
        %0 = bmodelica.constant #bmodelica<real 1.1>
        bmodelica.yield %0 : !r
    }
    bmodelica.binding_equation @beta {
        %0 = bmodelica.constant #bmodelica<real 0.4>
        bmodelica.yield %0 : !r
    }
    bmodelica.binding_equation @gamma {
        %0 = bmodelica.constant #bmodelica<real 0.4>
        bmodelica.yield %0 : !r
    }

    // der(rabbit) = alpha*rabbit - beta*rabbit*fox
    %t0 = bmodelica.equation_template inductions = [] {
        %x = bmodelica.variable.get @rabbit : !r
        %y = bmodelica.variable.get @fox : !r
        %dx = bmodelica.der %x : !r -> !r
        %alpha = bmodelica.variable.get @alpha : !r
        %beta = bmodelica.variable.get @beta : !r
        %0 = bmodelica.mul %alpha, %x : !r, !r -> !r
        %1 = bmodelica.mul %beta, %x : !r, !r -> !r
        %2 = bmodelica.mul %1, %y : !r, !r -> !r
        %3 = bmodelica.sub %0, %2 : !r, !r -> !r
        %l = bmodelica.equation_side %dx : tuple<!r>
        %r = bmodelica.equation_side %3 : tuple<!r>
        bmodelica.equation_sides %l, %r : tuple<!r>, tuple<!r>
    }

    // der(fox) = beta*rabbit*fox - gamma*fox
    %t1 = bmodelica.equation_template inductions = [] {
        %x = bmodelica.variable.get @rabbit : !r
        %y = bmodelica.variable.get @fox : !r
        %dy = bmodelica.der %y : !r -> !r
        %beta = bmodelica.variable.get @beta : !r
        %gamma = bmodelica.variable.get @gamma : !r
        %0 = bmodelica.mul %beta, %x : !r, !r -> !r
        %1 = bmodelica.mul %0, %y : !r, !r -> !r
        %2 = bmodelica.mul %gamma, %y : !r, !r -> !r
        %3 = bmodelica.sub %1, %2 : !r, !r -> !r
        %l = bmodelica.equation_side %dy : tuple<!r>
        %r = bmodelica.equation_side %3 : tuple<!r>
        bmodelica.equation_sides %l, %r : tuple<!r>, tuple<!r>
    }

    bmodelica.dynamic {
        bmodelica.equation_instance %t0
        bmodelica.equation_instance %t1
    }
}

// CHECK-LABEL: acr.crn @LotkaVolterra
// CHECK:         acr.crn_generic
// CHECK-SAME:      initial_values = array<f64: 2.000000e+00, 1.000000e+00>
// CHECK-SAME:      rules = ["predation", "autocatalysis", "degradation"]
// CHECK-SAME:      species_symbols = ["rabbit", "fox"]

// X -> 2X keeps its catalyst on the way to the emulator.
// JSON:      "species": [
// JSON-NEXT:   "rabbit",
// JSON-NEXT:   "fox"
// JSON:      "species": [
// JSON:      "reactants": [
// JSON-NEXT:   1,
// JSON-NEXT:   0
// JSON:      "products": [
// JSON-NEXT:   2,
// JSON-NEXT:   0
// JSON:      "rate": 1.1,

// SIM: PASS
