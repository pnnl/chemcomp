// RUN: standalone-translate -mlir-to-json %s | FileCheck %s

// Tensor form with constant rates: 3 species, 2 reactions

acr.crn @tensor_net(%sv: tensor<3x!acr.species>, %rv: tensor<2xf32>)
    -> tensor<3x!acr.species> {
  %rates = arith.constant dense<[0.1, 0.05]> : tensor<2xf32>
  %out = acr.crn_generic
      {stoich_matrix = dense<[[-2, -1, 2], [0, -1, 1]]> : tensor<2x3xi32>}
      (%sv, %rates)
      : (tensor<3x!acr.species>, tensor<2xf32>) -> tensor<3x!acr.species>
  acr.return %out : tensor<3x!acr.species>
}

// CHECK:      {
// CHECK-NEXT:     "reactions": [
// CHECK-NEXT:         {
// CHECK-NEXT:             "species": [
// CHECK-NEXT:                 "S0",
// CHECK-NEXT:                 "S1",
// CHECK-NEXT:                 "S2"
// CHECK-NEXT:             ],
// CHECK-NEXT:             "reactants": [
// CHECK-NEXT:                 2,
// CHECK-NEXT:                 1,
// CHECK-NEXT:                 0
// CHECK-NEXT:             ],
// CHECK-NEXT:             "products": [
// CHECK-NEXT:                 0,
// CHECK-NEXT:                 0,
// CHECK-NEXT:                 2
// CHECK-NEXT:             ],
// CHECK-NEXT:             "rate": 0.1,
// CHECK-NEXT:             "name": "R1"
// CHECK-NEXT:         },
// CHECK-NEXT:         {
// CHECK-NEXT:             "species": [
// CHECK-NEXT:                 "S0",
// CHECK-NEXT:                 "S1",
// CHECK-NEXT:                 "S2"
// CHECK-NEXT:             ],
// CHECK-NEXT:             "reactants": [
// CHECK-NEXT:                 0,
// CHECK-NEXT:                 1,
// CHECK-NEXT:                 0
// CHECK-NEXT:             ],
// CHECK-NEXT:             "products": [
// CHECK-NEXT:                 0,
// CHECK-NEXT:                 0,
// CHECK-NEXT:                 1
// CHECK-NEXT:             ],
// CHECK-NEXT:             "rate": 0.05,
// CHECK-NEXT:             "name": "R2"
// CHECK-NEXT:         }
// CHECK-NEXT:     ],
// CHECK-NEXT:     "populations": {
// CHECK-NEXT:         "S0": 1,
// CHECK-NEXT:         "S1": 1,
// CHECK-NEXT:         "S2": 1
// CHECK-NEXT:     },
// CHECK-NEXT:     "tracked_species": []
// CHECK-NEXT: }
