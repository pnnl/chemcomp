// Reversible isomerization A <=> B, in the form produced by
// derivatives-materialization (der_a, der_b are variables). The equations
// appear in reverse order and with commuted products: matching works on the
// normalized polynomial system.

// RUN: chemcomp-opt %s --convert-bmodelica-to-acr | FileCheck %s

!r = !bmodelica.real

bmodelica.model @Isomerization der = [<@a, @der_a>, <@b, @der_b>] {
    bmodelica.variable @a : !bmodelica.variable<!r>
    bmodelica.variable @b : !bmodelica.variable<!r>
    bmodelica.variable @der_a : !bmodelica.variable<!r>
    bmodelica.variable @der_b : !bmodelica.variable<!r>
    bmodelica.variable @kf : !bmodelica.variable<!r, parameter>
    bmodelica.variable @kr : !bmodelica.variable<!r, parameter>

    bmodelica.binding_equation @kf {
        %0 = bmodelica.constant #bmodelica<real 3.0>
        bmodelica.yield %0 : !r
    }

    bmodelica.start @kr {
        %0 = bmodelica.constant #bmodelica<real 1.5>
        bmodelica.yield %0 : !r
    } {each = false, fixed = true}

    // der(a) = kr*b - kf*a
    %t0 = bmodelica.equation_template inductions = [] {
        %a = bmodelica.variable.get @a : !r
        %b = bmodelica.variable.get @b : !r
        %da = bmodelica.variable.get @der_a : !r
        %kf = bmodelica.variable.get @kf : !r
        %kr = bmodelica.variable.get @kr : !r
        %0 = bmodelica.mul %kr, %b : !r, !r -> !r
        %1 = bmodelica.mul %kf, %a : !r, !r -> !r
        %2 = bmodelica.sub %0, %1 : !r, !r -> !r
        %l = bmodelica.equation_side %da : tuple<!r>
        %r = bmodelica.equation_side %2 : tuple<!r>
        bmodelica.equation_sides %l, %r : tuple<!r>, tuple<!r>
    }

    // der(b) = a*kf - b*kr
    %t1 = bmodelica.equation_template inductions = [] {
        %a = bmodelica.variable.get @a : !r
        %b = bmodelica.variable.get @b : !r
        %db = bmodelica.variable.get @der_b : !r
        %kf = bmodelica.variable.get @kf : !r
        %kr = bmodelica.variable.get @kr : !r
        %0 = bmodelica.mul %a, %kf : !r, !r -> !r
        %1 = bmodelica.mul %b, %kr : !r, !r -> !r
        %2 = bmodelica.sub %0, %1 : !r, !r -> !r
        %l = bmodelica.equation_side %db : tuple<!r>
        %r = bmodelica.equation_side %2 : tuple<!r>
        bmodelica.equation_sides %l, %r : tuple<!r>, tuple<!r>
    }

    bmodelica.dynamic {
        bmodelica.equation_instance %t1
        bmodelica.equation_instance %t0
    }
}

// Reaction 0 is A -> B (rate kf), reaction 1 is B -> A (rate kr). The rate of
// kr comes from its start value, the one of kf from its binding equation.

// CHECK-LABEL: acr.crn @Isomerization(%{{.*}}: tensor<2x!acr.species>, %{{.*}}: tensor<2xf32>)
// CHECK:         acr.crn_generic
// CHECK-SAME:      crn_rule = "reversible-isomerization"
// CHECK-SAME:      ode_state_names = ["a", "b"]
// CHECK-SAME:      rate_names = ["kf", "kr"]
// CHECK-SAME:      rate_values = array<f64: 3.000000e+00, 1.500000e+00>
// CHECK-SAME:      stoich_matrix = dense<{{\[}}[-1, 1], [1, -1]]> : tensor<2x2xi32>
