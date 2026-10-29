// Lotka-Volterra predator-prey, in mass-action form:
//   der(rabbit) = alpha*rabbit - beta*rabbit*fox
//   der(fox)    = beta*rabbit*fox - gamma*fox
// maps onto
//   R0: X -> 2X       (alpha)   prey reproduction, catalytic in X
//   R1: X + Y -> 2Y   (beta)    predation
//   R2: Y -> 0        (gamma)   predator death
// with X = rabbit and Y = fox. The predator is declared first, so species
// columns (model order: fox, rabbit) differ from the rule order (X, Y).

// RUN: chemcomp-opt %s --convert-bmodelica-to-acr | FileCheck %s

!r = !bmodelica.real

bmodelica.model @LotkaVolterra {
    bmodelica.variable @fox : !bmodelica.variable<!r>
    bmodelica.variable @rabbit : !bmodelica.variable<!r>
    bmodelica.variable @alpha : !bmodelica.variable<!r, parameter>
    bmodelica.variable @beta : !bmodelica.variable<!r, parameter>
    bmodelica.variable @gamma : !bmodelica.variable<!r, parameter>

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

// Columns are [fox, rabbit].
// CHECK-LABEL: acr.crn @LotkaVolterra(%{{.*}}: tensor<2x!acr.species>, %{{.*}}: tensor<3xf32>)
// CHECK:         acr.crn_generic
// CHECK-SAME:      crn_rule = "lotka-volterra"
// CHECK-SAME:      ode_state_names = ["fox", "rabbit"]
// CHECK-SAME:      product_matrix = dense<{{\[}}[0, 2], [2, 0], [0, 0]]> : tensor<3x2xi32>
// CHECK-SAME:      rate_names = ["alpha", "beta", "gamma"]
// CHECK-SAME:      rate_values = array<f64: 1.100000e+00, 4.000000e-01, 4.000000e-01>
// CHECK-SAME:      reactant_matrix = dense<{{\[}}[0, 1], [1, 1], [1, 0]]> : tensor<3x2xi32>
// The net stoichiometry loses the catalytic X in R0 (X -> 2X is just +X).
// CHECK-SAME:      stoich_matrix = dense<{{\[}}[0, 1], [1, -1], [-1, 0]]> : tensor<3x2xi32>
