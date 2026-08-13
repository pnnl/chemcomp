// Tie-breaking and determinism.
//
// The paper specifies no tie-breaking rule, but FileCheck compares exact
// output, so one is imposed: higher score, then more matched reactions, then
// the lexicographically smaller motif-id sequence, then the smaller binding.
// Candidate motifs are always visited in ascending id order and every domain
// is a sorted vector, so nothing observable depends on hash order.
//
// AmbiguousSlice is five motifs sharing signature (-1,-1,1,1) -- the shape
// holding 11,421 of 43,455 motifs in the real database.

// Running twice in one pipeline must produce identical output; the second
// invocation re-matches an already-annotated operation.
// RUN: chemcomp-opt %s -crn-db-match='db-slice=AmbiguousSlice' \
// RUN:   -crn-db-match='db-slice=AmbiguousSlice' | FileCheck %s

// Scoring weights are reported the same way regardless of their values: with
// every weight zeroed all matches tie, and the tie-break still picks rxnamb01.
// RUN: chemcomp-opt %s -crn-db-match='db-slice=AmbiguousSlice wc=0 ws=0 wr=0' \
// RUN:   | FileCheck %s

acr.crn @ambiguous(%sv: tensor<4x!acr.species>, %rv: tensor<1xf32>)
    -> tensor<4x!acr.species> {
  %out = acr.crn_generic
      {stoich_matrix = dense<[[-1, -1, 1, 1]]> : tensor<1x4xi32>}
      (%sv, %rv) : (tensor<4x!acr.species>, tensor<1xf32>)
                -> tensor<4x!acr.species>
  acr.return %out : tensor<4x!acr.species>
}

// The lowest motif id wins the tie.
// CHECK: acr.crn_generic
// CHECK-SAME: motif_ids = ["rxnamb01"]
// CHECK-SAME: species_names = ["ATP", "H2O", "ADP", "Phosphate"]
