// Base Modelica to ode.system, before and after MARCO's own
// derivatives-materialization: both give the same system.

// RUN: chemcomp-opt %s --convert-bmodelica-to-ode | FileCheck %s
// RUN: chemcomp-opt %s --derivatives-materialization --convert-bmodelica-to-ode \
// RUN:   | FileCheck %s
// RUN: chemcomp-opt %s --convert-bmodelica-to-ode=keep-model \
// RUN:   | FileCheck %s --check-prefix=KEEP

!r = !bmodelica.real

// Modelica:
//   model Decay
//     Real x(start = 1.0, fixed = true);
//     parameter Real k = 0.5;
//   equation
//     der(x) = -k*x;
//   end Decay;
bmodelica.model @Decay {
    bmodelica.variable @x : !bmodelica.variable<!r>
    bmodelica.variable @k : !bmodelica.variable<!r, parameter>

    bmodelica.start @x {
        %0 = bmodelica.constant #bmodelica<real 1.0>
        bmodelica.yield %0 : !r
    } {each = false, fixed = true}

    bmodelica.binding_equation @k {
        %0 = bmodelica.constant #bmodelica<real 0.5>
        bmodelica.yield %0 : !r
    }

    %t0 = bmodelica.equation_template inductions = [] {
        %x = bmodelica.variable.get @x : !r
        %dx = bmodelica.der %x : !r -> !r
        %k = bmodelica.variable.get @k : !r
        %kx = bmodelica.mul %k, %x : !r, !r -> !r
        %rhs = bmodelica.neg %kx : !r -> !r
        %l = bmodelica.equation_side %dx : tuple<!r>
        %r = bmodelica.equation_side %rhs : tuple<!r>
        bmodelica.equation_sides %l, %r : tuple<!r>, tuple<!r>
    }

    bmodelica.dynamic {
        bmodelica.equation_instance %t0
    }
}

// CHECK-NOT:   bmodelica.model
// CHECK-LABEL: ode.system @Decay {
// CHECK-NEXT:    ode.state @x initial 1.0
// CHECK-NEXT:    ode.param @k = 0.5
// CHECK-NEXT:    ode.term @x += -1.0 params [@k] states [@x]
// CHECK-NEXT:  }

// KEEP: ode.system @Decay_ode
// KEEP: bmodelica.model @Decay
