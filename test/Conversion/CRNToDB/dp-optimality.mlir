// The dynamic program must beat a greedy per-reaction matcher.
//
//   r0: s0 -> s1        r1: s0 -> s2
//
// Both reactions consume the shared species s0, so whatever s0 binds to must
// serve both. The database offers:
//
//   m01: X -> P     m02: Y -> Q     m03: Y -> R
//
// A greedy matcher processes r0 first and takes m01, the lowest motif id,
// binding s0 = X. Nothing then serves r1, because only m01 consumes X and its
// product P is already spoken for -- coverage 1 of 2.
//
// The dynamic program instead considers the consequences for r1 before
// committing, gives up the locally preferred m01, and reaches full coverage
// with m02 and m03 sharing s0 = Y.

// RUN: chemcomp-opt %s \
// RUN:   -crn-db-match='db=%S/Inputs/lookahead-db.json global-injectivity=1' \
// RUN:   | FileCheck %s

// Full coverage, so no diagnostics.
// RUN: chemcomp-opt %s \
// RUN:   -crn-db-match='db=%S/Inputs/lookahead-db.json global-injectivity=1' \
// RUN:   -verify-diagnostics

acr.crn @lookahead(%sv: tensor<3x!acr.species>, %rv: tensor<2xf32>)
    -> tensor<3x!acr.species> {
  %out = acr.crn_generic
      {stoich_matrix = dense<[[-1, 1, 0],
                              [-1, 0, 1]]> : tensor<2x3xi32>}
      (%sv, %rv) : (tensor<3x!acr.species>, tensor<2xf32>)
                -> tensor<3x!acr.species>
  acr.return %out : tensor<3x!acr.species>
}

// m01 is never chosen despite having the lowest id: full coverage wins.
// CHECK: acr.crn_generic
// CHECK-SAME: motif_ids = ["m02", "m03"]
// CHECK-SAME: species_names = ["Y", "Q", "R"]
