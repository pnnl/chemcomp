// The same Lotka-Volterra system written differently: sides swapped, the
// derivative scaled, terms reordered and a numeric division. Matching works
// on the normalized polynomials, so the network is the one of
// lotka-volterra.mlir, with rates read off the coefficients.

// RUN: chemcomp-opt %s --convert-bmodelica-to-acr | FileCheck %s

!r = !bmodelica.real

bmodelica.model @Rewritten {
    bmodelica.variable @x : !bmodelica.variable<!r>
    bmodelica.variable @y : !bmodelica.variable<!r>
    bmodelica.variable @b : !bmodelica.variable<!r, parameter>

    // 2*x - b*(y*x) = der(x)
    %t0 = bmodelica.equation_template inductions = [] {
        %x = bmodelica.variable.get @x : !r
        %y = bmodelica.variable.get @y : !r
        %b = bmodelica.variable.get @b : !r
        %dx = bmodelica.der %x : !r -> !r
        %two = bmodelica.constant #bmodelica<real 2.0>
        %0 = bmodelica.mul %two, %x : !r, !r -> !r
        %1 = bmodelica.mul %y, %x : !r, !r -> !r
        %2 = bmodelica.mul %b, %1 : !r, !r -> !r
        %3 = bmodelica.sub %0, %2 : !r, !r -> !r
        %l = bmodelica.equation_side %3 : tuple<!r>
        %r = bmodelica.equation_side %dx : tuple<!r>
        bmodelica.equation_sides %l, %r : tuple<!r>, tuple<!r>
    }

    // 2*der(y) = -(y) + 2*(x*b)*y     i.e. der(y) = b*x*y - y/2
    %t1 = bmodelica.equation_template inductions = [] {
        %x = bmodelica.variable.get @x : !r
        %y = bmodelica.variable.get @y : !r
        %b = bmodelica.variable.get @b : !r
        %dy = bmodelica.der %y : !r -> !r
        %two = bmodelica.constant #bmodelica<int 2>
        %0 = bmodelica.mul %two, %dy : !bmodelica.int, !r -> !r
        %1 = bmodelica.neg %y : !r -> !r
        %2 = bmodelica.mul %x, %b : !r, !r -> !r
        %3 = bmodelica.mul %two, %2 : !bmodelica.int, !r -> !r
        %4 = bmodelica.mul %3, %y : !r, !r -> !r
        %5 = bmodelica.add %1, %4 : !r, !r -> !r
        %l = bmodelica.equation_side %0 : tuple<!r>
        %r = bmodelica.equation_side %5 : tuple<!r>
        bmodelica.equation_sides %l, %r : tuple<!r>, tuple<!r>
    }

    bmodelica.dynamic {
        bmodelica.equation_instance %t0
        bmodelica.equation_instance %t1
    }
}

// CHECK-LABEL: acr.crn @Rewritten
// CHECK:         acr.crn_generic
// CHECK-SAME:      crn_rule = "lotka-volterra"
// CHECK-SAME:      ode_state_names = ["x", "y"]
// CHECK-SAME:      rate_names = ["2", "b", "0.5"]
// CHECK-SAME:      stoich_matrix = dense<{{\[}}[1, 0], [-1, 1], [0, -1]]> : tensor<3x2xi32>
