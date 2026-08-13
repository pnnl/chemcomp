// Documents a gap between the CF26-CRN2DB problem statement and its algorithm.
//
// Section 3.3 requires phi to be globally injective: no two abstract species
// may map to one chemical. But the DP state of section 4.2 cannot express
// that. trim() drops a binding as soon as its species leaves the frontier, so
// a chemical retired at step t can be bound again at step t+1 with nothing to
// notice.
//
// This CRN triggers it. s1 leaves the frontier after r0, so when r1 binds s2
// the DP no longer remembers that P is taken, and both s1 and s2 end up as P.
//
// Default behaviour matches the original algorithm: frontier-local
// injectivity, with a warning so the violation is visible rather than
// silent. global-injectivity=1
// carries retired chemicals in the state key and produces a genuinely injective
// assignment, at the cost of a state space exponential in the species count
// rather than the treewidth -- which is why it is not the default.

// RUN: standalone-opt %s -crn-db-match='db=%S/Inputs/lookahead-db.json' \
// RUN:   | FileCheck %s --check-prefix=LOCAL
// RUN: standalone-opt %s -crn-db-match='db=%S/Inputs/lookahead-db.json' \
// RUN:   -verify-diagnostics
// RUN: standalone-opt %s \
// RUN:   -crn-db-match='db=%S/Inputs/lookahead-db.json global-injectivity=1' \
// RUN:   | FileCheck %s --check-prefix=GLOBAL

acr.crn @reuse(%sv: tensor<3x!acr.species>, %rv: tensor<2xf32>)
    -> tensor<3x!acr.species> {
  // expected-warning @+1 {{assignment is not globally injective; reused: P (pass global-injectivity=1 to enforce injectivity)}}
  %out = acr.crn_generic
      {stoich_matrix = dense<[[-1, 1, 0],
                              [-1, 0, 1]]> : tensor<2x3xi32>}
      (%sv, %rv) : (tensor<3x!acr.species>, tensor<2xf32>)
                -> tensor<3x!acr.species>
  acr.return %out : tensor<3x!acr.species>
}

// Frontier-local: P is bound twice, and m01 is reused for both reactions.
// LOCAL: acr.crn_generic
// LOCAL-SAME: motif_ids = ["m01", "m01"]
// LOCAL-SAME: species_names = ["X", "P", "P"]

// Enforced: a genuinely injective assignment, no chemical reused.
// GLOBAL: acr.crn_generic
// GLOBAL-SAME: motif_ids = ["m02", "m03"]
// GLOBAL-SAME: species_names = ["Y", "Q", "R"]
