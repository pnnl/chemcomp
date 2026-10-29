// With allow-unmatched, a model that matches no known network is reported as
// a warning and left in place, while the others are still converted.

// RUN: chemcomp-opt %s --convert-bmodelica-to-acr=allow-unmatched 2>&1 \
// RUN:   | FileCheck %s

!r = !bmodelica.real

// CHECK: warning: no known chemical reaction network matches the ODE system of model 'Quadratic'
// CHECK: note: ODE system: der(x) = -x^2
bmodelica.model @Quadratic {
    bmodelica.variable @x : !bmodelica.variable<!r>

    %t0 = bmodelica.equation_template inductions = [] {
        %x = bmodelica.variable.get @x : !r
        %dx = bmodelica.der %x : !r -> !r
        %0 = bmodelica.mul %x, %x : !r, !r -> !r
        %1 = bmodelica.neg %0 : !r -> !r
        %l = bmodelica.equation_side %dx : tuple<!r>
        %r = bmodelica.equation_side %1 : tuple<!r>
        bmodelica.equation_sides %l, %r : tuple<!r>, tuple<!r>
    }

    bmodelica.dynamic {
        bmodelica.equation_instance %t0
    }
}

bmodelica.model @Decay {
    bmodelica.variable @x : !bmodelica.variable<!r>

    %t0 = bmodelica.equation_template inductions = [] {
        %x = bmodelica.variable.get @x : !r
        %dx = bmodelica.der %x : !r -> !r
        %0 = bmodelica.neg %x : !r -> !r
        %l = bmodelica.equation_side %dx : tuple<!r>
        %r = bmodelica.equation_side %0 : tuple<!r>
        bmodelica.equation_sides %l, %r : tuple<!r>, tuple<!r>
    }

    bmodelica.dynamic {
        bmodelica.equation_instance %t0
    }
}

// CHECK: bmodelica.model @Quadratic
// CHECK: acr.crn @Decay
// CHECK: crn_rule = "decay"
// CHECK-SAME: rate_names = ["1"]
// CHECK-SAME: rate_values = array<f64: 1.000000e+00>
// CHECK-NOT: bmodelica.model @Decay
