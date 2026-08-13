// Elimination ordering and frontier computation on a chain CRN, where the
// optimal order is obvious by inspection.
//
// Four reactions over five species, each reaction sharing exactly one species
// with the next:
//
//   r0: s0 -> s1     r1: s1 -> s2     r2: s2 -> s3     r3: s3 -> s4
//
// A chain has treewidth 1, and processing the reactions in chain order keeps
// exactly one species live at each step. The rows below are deliberately
// declared out of chain order (r2, r0, r3, r1) so a correct ordering has to be
// derived rather than inherited from the input.

// RUN: standalone-opt %s -crn-db-match='db-slice=SubsetSlice dump-preprocessing=1' \
// RUN:   | FileCheck %s

acr.crn @chain(%sv: tensor<5x!acr.species>, %rv: tensor<4xf32>)
    -> tensor<5x!acr.species> {
  %out = acr.crn_generic
      {stoich_matrix = dense<[[ 0,  0, -1,  1,  0],
                              [-1,  1,  0,  0,  0],
                              [ 0,  0,  0, -1,  1],
                              [ 0, -1,  1,  0,  0]]> : tensor<4x5xi32>}
      (%sv, %rv) : (tensor<5x!acr.species>, tensor<4xf32>)
                -> tensor<5x!acr.species>
  acr.return %out : tensor<5x!acr.species>
}

// Declared order is r2, r0, r3, r1; the derived order is the chain order.
// CHECK: order: 1 3 0 2

// One species live per step: the chain's true treewidth.
// CHECK: max-frontier-width: 1
// CHECK-NEXT: frontier[0]:
// CHECK-NEXT: frontier[1]: 1
// CHECK-NEXT: frontier[2]: 2
// CHECK-NEXT: frontier[3]: 3
// CHECK-NEXT: frontier[4]:
