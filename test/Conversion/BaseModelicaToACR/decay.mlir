// Exponential decay, der(x) = -k*x, maps onto the single reaction X -> 0.

// The model as emitted by the MARCO frontend (bmodelica.der) ...
// RUN: chemcomp-opt %s --convert-bmodelica-to-acr | FileCheck %s

// ... and after MARCO's own derivatives-materialization pass both convert to
// the same network.
// RUN: chemcomp-opt %s --derivatives-materialization --convert-bmodelica-to-acr \
// RUN:   | FileCheck %s

// keep-model leaves the source model next to the network, which is then named
// after the model with a _crn suffix.
// RUN: chemcomp-opt %s --convert-bmodelica-to-acr=keep-model \
// RUN:   | FileCheck %s --check-prefix=KEEP

!r = !bmodelica.real

bmodelica.model @Decay {
    bmodelica.variable @x : !bmodelica.variable<!r>
    bmodelica.variable @k : !bmodelica.variable<!r, parameter>

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
// CHECK-LABEL: acr.crn @Decay(%{{.*}}: tensor<1x!acr.species>, %{{.*}}: tensor<1xf32>) -> tensor<1x!acr.species>
// CHECK:         acr.crn_generic
// CHECK-SAME:      crn_rule = "decay"
// CHECK-SAME:      ode_state_names = ["x"]
// CHECK-SAME:      product_matrix = dense<0> : tensor<1x1xi32>
// CHECK-SAME:      rate_names = ["k"]
// CHECK-SAME:      rate_values = array<f64: 5.000000e-01>
// CHECK-SAME:      reactant_matrix = dense<1> : tensor<1x1xi32>
// CHECK-SAME:      stoich_matrix = dense<-1> : tensor<1x1xi32>
// CHECK:         acr.return
// CHECK-NOT:   bmodelica.model

// KEEP: acr.crn @Decay_crn
// KEEP: bmodelica.model @Decay
