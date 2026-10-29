// Birth-death, der(x) = 2 - kd*x: constant production 0 -> X and decay X -> 0.
// A numeric coefficient becomes the rate constant itself. kd has no binding,
// so the rates have names but no values.

// RUN: chemcomp-opt %s --convert-bmodelica-to-acr | FileCheck %s

!r = !bmodelica.real

bmodelica.model @BirthDeath {
    bmodelica.variable @x : !bmodelica.variable<!r>
    bmodelica.variable @kd : !bmodelica.variable<!r, parameter>

    %t0 = bmodelica.equation_template inductions = [] {
        %x = bmodelica.variable.get @x : !r
        %dx = bmodelica.der %x : !r -> !r
        %kd = bmodelica.variable.get @kd : !r
        %c = bmodelica.constant #bmodelica<real 2.0>
        %kdx = bmodelica.mul %kd, %x : !r, !r -> !r
        %rhs = bmodelica.sub %c, %kdx : !r, !r -> !r
        %l = bmodelica.equation_side %dx : tuple<!r>
        %r = bmodelica.equation_side %rhs : tuple<!r>
        bmodelica.equation_sides %l, %r : tuple<!r>, tuple<!r>
    }

    bmodelica.dynamic {
        bmodelica.equation_instance %t0
    }
}

// CHECK-LABEL: acr.crn @BirthDeath(%{{.*}}: tensor<1x!acr.species>, %{{.*}}: tensor<2xf32>)
// CHECK:         acr.crn_generic
// CHECK-SAME:      crn_rule = "birth-death"
// CHECK-SAME:      product_matrix = dense<{{\[}}[1], [0]]> : tensor<2x1xi32>
// CHECK-SAME:      rate_names = ["2", "kd"]
// CHECK-NOT:       rate_values
// CHECK-SAME:      reactant_matrix = dense<{{\[}}[0], [1]]> : tensor<2x1xi32>
// CHECK-SAME:      stoich_matrix = dense<{{\[}}[1], [-1]]> : tensor<2x1xi32>
