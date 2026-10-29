// The full ODE-to-emulator flow: a Base Modelica ODE model is mapped onto a
// known network, crn-db-match assigns concrete chemicals, and the JSON export
// carries names and rate constants to the reservoir emulator.

// RUN: chemcomp-opt %s --convert-bmodelica-to-acr \
// RUN:   | chemcomp-opt -crn-db-match='db=%S/Inputs/isomerase-db.json' \
// RUN:   | chemcomp-translate -mlir-to-json | FileCheck %s

!r = !bmodelica.real

// der(g) = kr*f - kf*g, der(f) = kf*g - kr*f
bmodelica.model @Isomerase {
    bmodelica.variable @g : !bmodelica.variable<!r>
    bmodelica.variable @f : !bmodelica.variable<!r>
    bmodelica.variable @kf : !bmodelica.variable<!r, parameter>
    bmodelica.variable @kr : !bmodelica.variable<!r, parameter>

    bmodelica.binding_equation @kf {
        %0 = bmodelica.constant #bmodelica<real 2.0>
        bmodelica.yield %0 : !r
    }
    bmodelica.binding_equation @kr {
        %0 = bmodelica.constant #bmodelica<real 0.5>
        bmodelica.yield %0 : !r
    }

    %t0 = bmodelica.equation_template inductions = [] {
        %g = bmodelica.variable.get @g : !r
        %f = bmodelica.variable.get @f : !r
        %dg = bmodelica.der %g : !r -> !r
        %kf = bmodelica.variable.get @kf : !r
        %kr = bmodelica.variable.get @kr : !r
        %0 = bmodelica.mul %kr, %f : !r, !r -> !r
        %1 = bmodelica.mul %kf, %g : !r, !r -> !r
        %2 = bmodelica.sub %0, %1 : !r, !r -> !r
        %l = bmodelica.equation_side %dg : tuple<!r>
        %r = bmodelica.equation_side %2 : tuple<!r>
        bmodelica.equation_sides %l, %r : tuple<!r>, tuple<!r>
    }

    %t1 = bmodelica.equation_template inductions = [] {
        %g = bmodelica.variable.get @g : !r
        %f = bmodelica.variable.get @f : !r
        %df = bmodelica.der %f : !r -> !r
        %kf = bmodelica.variable.get @kf : !r
        %kr = bmodelica.variable.get @kr : !r
        %0 = bmodelica.mul %kf, %g : !r, !r -> !r
        %1 = bmodelica.mul %kr, %f : !r, !r -> !r
        %2 = bmodelica.sub %0, %1 : !r, !r -> !r
        %l = bmodelica.equation_side %df : tuple<!r>
        %r = bmodelica.equation_side %2 : tuple<!r>
        bmodelica.equation_sides %l, %r : tuple<!r>, tuple<!r>
    }

    bmodelica.dynamic {
        bmodelica.equation_instance %t0
        bmodelica.equation_instance %t1
    }
}

// R1: G6P -> F6P (kf = 2), R2: F6P -> G6P (kr = 0.5).
// CHECK:      "reactions": [
// CHECK:        "species": [
// CHECK-NEXT:     "Glucose-6-phosphate",
// CHECK-NEXT:     "Fructose-6-phosphate"
// CHECK:        "reactants": [
// CHECK-NEXT:     1,
// CHECK-NEXT:     0
// CHECK:        "rate": 2,
// CHECK:        "reactants": [
// CHECK-NEXT:     0,
// CHECK-NEXT:     1
// CHECK:        "rate": 0.5,
