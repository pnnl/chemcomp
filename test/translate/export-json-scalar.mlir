// RUN: standalone-translate -mlir-to-json %s | FileCheck %s

// Scalar form: 3 species, 2 reactions with constant rates

acr.crn @scalar_net(%s0: !acr.species, %s1: !acr.species, %s2: !acr.species,
                     %r0: f32, %r1: f32)
    -> (!acr.species, !acr.species, !acr.species) {
  %c0 = arith.constant 0.1 : f32
  %c1 = arith.constant 0.05 : f32
  %out0, %out1, %out2 = acr.crn_generic
      {stoich_matrix = dense<[[-2, -1, 2], [0, -1, 1]]> : tensor<2x3xi32>}
      (%s0, %s1, %s2, %c0, %c1)
          : (!acr.species, !acr.species, !acr.species, f32, f32)
         -> (!acr.species, !acr.species, !acr.species)
  acr.return %out0, %out1, %out2 : !acr.species, !acr.species, !acr.species
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
