// RUN: chemcomp-opt %s --acr-network-to-crn | FileCheck %s
// RUN: chemcomp-opt %s --acr-network-to-crn | chemcomp-translate -mlir-to-json \
// RUN:   | FileCheck %s --check-prefix=JSON

// The graph form lowers to the acr.crn / acr.crn_generic form, keeping
// species and reaction order. The full matrices keep the catalyst of
// X -> 2X, which the net stoich_matrix loses.

// CHECK-NOT:   acr.network
// CHECK-LABEL: acr.crn @LV(%{{.*}}: tensor<3x!acr.species>, %{{.*}}: tensor<3xf32>) -> tensor<3x!acr.species>
// CHECK:         acr.crn_generic
// CHECK-SAME:      initial_values = array<f64: 1.000000e+00, 0.000000e+00, 2.000000e+00>
// CHECK-SAME:      ode_state_names = ["x", "x", "y"]
// CHECK-SAME:      ode_state_signs = array<i64: 1, -1, 1>
// CHECK-SAME:      product_matrix = dense<{{\[}}[2, 0, 0], [0, 0, 2], [0, 0, 0]]> : tensor<3x3xi32>
// CHECK-SAME:      rate_names = ["a", "b", "0.4"]
// CHECK-SAME:      rate_values = array<f64: 1.100000e+00, 4.000000e-01, 4.000000e-01>
// CHECK-SAME:      reactant_matrix = dense<{{\[}}[1, 0, 0], [1, 0, 1], [0, 0, 1]]> : tensor<3x3xi32>
// CHECK-SAME:      rules = ["autocatalysis", "predation", ""]
// CHECK-SAME:      species_symbols = ["x_p", "x_m", "y"]
// CHECK-SAME:      stoich_matrix = dense<{{\[}}[1, 0, 0], [-1, 0, 1], [0, 0, -1]]> : tensor<3x3xi32>
// CHECK:         acr.return
acr.network @LV {
  acr.species @x_p initial 1.0 {ode_state = "x", rail_sign = 1 : i64}
  acr.species @x_m initial 0.0 {ode_state = "x", rail_sign = -1 : i64}
  acr.species @y initial 2.0
  acr.rxn @x_p -> 2 @x_p rate 1.1 {rate_name = "a", rule = "autocatalysis"}
  acr.rxn @x_p + @y -> 2 @y rate 0.4 {rate_name = "b", rule = "predation"}
  acr.rxn @y -> none rate 0.4
}

// JSON:      "species": [
// JSON-NEXT:   "x_p",
// JSON-NEXT:   "x_m",
// JSON-NEXT:   "y"
// JSON:      "populations": {
// JSON-NEXT:   "x_p": 1,
// JSON-NEXT:   "x_m": 0,
// JSON-NEXT:   "y": 2
// JSON:      "species_map": {
// JSON-NEXT:   "x_p": {
// JSON-NEXT:     "state": "x",
// JSON-NEXT:     "sign": 1
