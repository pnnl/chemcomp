// The Lorenz system from Base Modelica. -x*z in der(y) is a negative
// cross-effect, so the pipeline dual-rail encodes it before composing.

// RUN: chemcomp-opt %s --convert-bmodelica-to-ode > %t.ode.mlir
// RUN: FileCheck %s --check-prefix=ODE < %t.ode.mlir
// RUN: chemcomp-opt %t.ode.mlir --ode-to-crn | FileCheck %s
// RUN: chemcomp-opt %t.ode.mlir --ode-to-crn | chemcomp-translate -mlir-to-json \
// RUN:   > %t.crn.json
// RUN: chemcomp-translate -ode-to-json %t.ode.mlir > %t.ode.json
// RUN: %if emulator %{ %chemcomp_python %check_realization --ode %t.ode.json \
// RUN:   --crn %t.crn.json --t-end 1 | FileCheck %s --check-prefix=SIM %}

!r = !bmodelica.real

bmodelica.model @Lorenz {
    bmodelica.variable @x : !bmodelica.variable<!r>
    bmodelica.variable @y : !bmodelica.variable<!r>
    bmodelica.variable @z : !bmodelica.variable<!r>
    bmodelica.variable @sigma : !bmodelica.variable<!r, parameter>
    bmodelica.variable @rho : !bmodelica.variable<!r, parameter>
    bmodelica.variable @beta : !bmodelica.variable<!r, parameter>

    bmodelica.start @x {
        %0 = bmodelica.constant #bmodelica<real 1.0>
        bmodelica.yield %0 : !r
    } {each = false, fixed = true}
    bmodelica.start @y {
        %0 = bmodelica.constant #bmodelica<real 1.0>
        bmodelica.yield %0 : !r
    } {each = false, fixed = true}
    bmodelica.start @z {
        %0 = bmodelica.constant #bmodelica<real 1.0>
        bmodelica.yield %0 : !r
    } {each = false, fixed = true}
    bmodelica.binding_equation @sigma {
        %0 = bmodelica.constant #bmodelica<real 10.0>
        bmodelica.yield %0 : !r
    }
    bmodelica.binding_equation @rho {
        %0 = bmodelica.constant #bmodelica<real 28.0>
        bmodelica.yield %0 : !r
    }
    bmodelica.binding_equation @beta {
        %0 = bmodelica.constant #bmodelica<real 2.6666666666666665>
        bmodelica.yield %0 : !r
    }

    // der(x) = sigma*(y - x)
    %t0 = bmodelica.equation_template inductions = [] {
        %x = bmodelica.variable.get @x : !r
        %y = bmodelica.variable.get @y : !r
        %dx = bmodelica.der %x : !r -> !r
        %sigma = bmodelica.variable.get @sigma : !r
        %0 = bmodelica.sub %y, %x : !r, !r -> !r
        %1 = bmodelica.mul %sigma, %0 : !r, !r -> !r
        %l = bmodelica.equation_side %dx : tuple<!r>
        %r = bmodelica.equation_side %1 : tuple<!r>
        bmodelica.equation_sides %l, %r : tuple<!r>, tuple<!r>
    }

    // der(y) = x*(rho - z) - y
    %t1 = bmodelica.equation_template inductions = [] {
        %x = bmodelica.variable.get @x : !r
        %y = bmodelica.variable.get @y : !r
        %z = bmodelica.variable.get @z : !r
        %dy = bmodelica.der %y : !r -> !r
        %rho = bmodelica.variable.get @rho : !r
        %0 = bmodelica.sub %rho, %z : !r, !r -> !r
        %1 = bmodelica.mul %x, %0 : !r, !r -> !r
        %2 = bmodelica.sub %1, %y : !r, !r -> !r
        %l = bmodelica.equation_side %dy : tuple<!r>
        %r = bmodelica.equation_side %2 : tuple<!r>
        bmodelica.equation_sides %l, %r : tuple<!r>, tuple<!r>
    }

    // der(z) = x*y - beta*z
    %t2 = bmodelica.equation_template inductions = [] {
        %x = bmodelica.variable.get @x : !r
        %y = bmodelica.variable.get @y : !r
        %z = bmodelica.variable.get @z : !r
        %dz = bmodelica.der %z : !r -> !r
        %beta = bmodelica.variable.get @beta : !r
        %0 = bmodelica.mul %x, %y : !r, !r -> !r
        %1 = bmodelica.mul %beta, %z : !r, !r -> !r
        %2 = bmodelica.sub %0, %1 : !r, !r -> !r
        %l = bmodelica.equation_side %dz : tuple<!r>
        %r = bmodelica.equation_side %2 : tuple<!r>
        bmodelica.equation_sides %l, %r : tuple<!r>, tuple<!r>
    }

    bmodelica.dynamic {
        bmodelica.equation_instance %t0
        bmodelica.equation_instance %t1
        bmodelica.equation_instance %t2
    }
}

// The products sigma*(y - x) and x*(rho - z) are expanded into monomials.
// ODE-LABEL: ode.system @Lorenz {
// ODE-NEXT:    ode.state @x initial 1.0
// ODE-NEXT:    ode.state @y initial 1.0
// ODE-NEXT:    ode.state @z initial 1.0
// ODE-NEXT:    ode.param @beta = 2.6666666666666665
// ODE-NEXT:    ode.param @rho = 28.0
// ODE-NEXT:    ode.param @sigma = 10.0
// ODE-NEXT:    ode.term @x += -1.0 params [@sigma] states [@x]
// ODE-NEXT:    ode.term @x += 1.0 params [@sigma] states [@y]
// ODE-NEXT:    ode.term @y += 1.0 params [@rho] states [@x]
// ODE-NEXT:    ode.term @y += -1.0 states [@x, @z]
// ODE-NEXT:    ode.term @y += -1.0 states [@y]
// ODE-NEXT:    ode.term @z += 1.0 states [@x, @y]
// ODE-NEXT:    ode.term @z += -1.0 params [@beta] states [@z]
// ODE-NEXT:  }

// CHECK-LABEL: acr.crn @Lorenz
// CHECK:         ode_state_names = ["x", "x", "y", "y", "z", "z"]
// CHECK-SAME:    species_symbols = ["x_p", "x_m", "y_p", "y_m", "z_p", "z_m"]

// SIM: PASS
