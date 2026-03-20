// RUN: standalone-translate --split-input-file -mlir-to-json %s | FileCheck %s

// Chemical Reactions (2 species, 2 reactions):
//   2 S0 + S1 -> 2 S2    (stoich row: [-2, -1, 2, 0, 0])
//   S1 + S3 -> S4         (stoich row: [0, -1, 0, -1, 1])

acr.crn @net(%sv: tensor<5x!acr.species>, %rv: tensor<2xf32>)
    -> tensor<5x!acr.species> {
  %out = acr.crn_generic
      {stoich_matrix = dense<[[-2, -1, 2, 0, 0],
                               [0, -1, 0, -1, 1]]> : tensor<2x5xi32>}
      (%sv, %rv)
      : (tensor<5x!acr.species>, tensor<2xf32>) -> tensor<5x!acr.species>
  acr.return %out : tensor<5x!acr.species>
}

// CHECK:      {
// CHECK-NEXT:     "reactions": [
// CHECK-NEXT:         {
// CHECK-NEXT:             "species": [
// CHECK-NEXT:                 "S0",
// CHECK-NEXT:                 "S1",
// CHECK-NEXT:                 "S2",
// CHECK-NEXT:                 "S3",
// CHECK-NEXT:                 "S4"
// CHECK-NEXT:             ],
// CHECK-NEXT:             "reactants": [
// CHECK-NEXT:                 2,
// CHECK-NEXT:                 1,
// CHECK-NEXT:                 0,
// CHECK-NEXT:                 0,
// CHECK-NEXT:                 0
// CHECK-NEXT:             ],
// CHECK-NEXT:             "products": [
// CHECK-NEXT:                 0,
// CHECK-NEXT:                 0,
// CHECK-NEXT:                 2,
// CHECK-NEXT:                 0,
// CHECK-NEXT:                 0
// CHECK-NEXT:             ],
// CHECK-NEXT:             "rate": 0,
// CHECK-NEXT:             "name": "R1"
// CHECK-NEXT:         },
// CHECK-NEXT:         {
// CHECK-NEXT:             "species": [
// CHECK-NEXT:                 "S0",
// CHECK-NEXT:                 "S1",
// CHECK-NEXT:                 "S2",
// CHECK-NEXT:                 "S3",
// CHECK-NEXT:                 "S4"
// CHECK-NEXT:             ],
// CHECK-NEXT:             "reactants": [
// CHECK-NEXT:                 0,
// CHECK-NEXT:                 1,
// CHECK-NEXT:                 0,
// CHECK-NEXT:                 1,
// CHECK-NEXT:                 0
// CHECK-NEXT:             ],
// CHECK-NEXT:             "products": [
// CHECK-NEXT:                 0,
// CHECK-NEXT:                 0,
// CHECK-NEXT:                 0,
// CHECK-NEXT:                 0,
// CHECK-NEXT:                 1
// CHECK-NEXT:             ],
// CHECK-NEXT:             "rate": 0,
// CHECK-NEXT:             "name": "R2"
// CHECK-NEXT:         }
// CHECK-NEXT:     ],
// CHECK-NEXT:     "populations": {
// CHECK-NEXT:         "S0": 1,
// CHECK-NEXT:         "S1": 1,
// CHECK-NEXT:         "S2": 1,
// CHECK-NEXT:         "S3": 1,
// CHECK-NEXT:         "S4": 1
// CHECK-NEXT:     },
// CHECK-NEXT:     "tracked_species": []
// CHECK-NEXT: }
