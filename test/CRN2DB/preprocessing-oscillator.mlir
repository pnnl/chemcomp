// The CF26-CRN2DB section 5.4 walkthrough, as a conformance test.
//
//     A + B -> C      C -> D      D + E -> A + F      F -> B + E
//
// Species columns are A B C D E F; reaction rows are the four above.
//
// REVISED: section 5.4 assumes this CRN has treewidth W = 2, reasoning that
// "species A, B, D, E each appear in two reactions". The frontier width
// defined in section 4.1 counts species simultaneously *live* -- seen so far
// and still to come -- not occurrences per species. Measuring against that
// definition, exhaustive search over all 24 reaction orderings gives a
// minimum achievable maximum frontier width of 3. The min-fill ordering below
// achieves that measured optimum.

// RUN: standalone-opt %s -crn-db-match='db-slice=OscillatorSlice dump-preprocessing=1' \
// RUN:   | FileCheck %s

acr.crn @oscillator(%sv: tensor<6x!acr.species>, %rv: tensor<4xf32>)
    -> tensor<6x!acr.species> {
  %out = acr.crn_generic
      {stoich_matrix = dense<[[-1, -1,  1,  0,  0,  0],
                              [ 0,  0, -1,  1,  0,  0],
                              [ 1,  0,  0, -1, -1,  1],
                              [ 0,  1,  0,  0,  1, -1]]> : tensor<4x6xi32>}
      (%sv, %rv) : (tensor<6x!acr.species>, tensor<4xf32>)
                -> tensor<6x!acr.species>
  acr.return %out : tensor<6x!acr.species>
}

// CHECK: order: 0 1 2 3
// CHECK: max-frontier-width: 3

// F_0 and F_m are empty by construction: nothing is live before the first
// reaction or after the last.
// CHECK-NEXT: frontier[0]:
// CHECK-NEXT: frontier[1]: 0 1 2
// CHECK-NEXT: frontier[2]: 0 1 3
// CHECK-NEXT: frontier[3]: 1 4 5
// CHECK-NEXT: frontier[4]:

// AC-3 prunes every domain to a single candidate, recovering exactly the
// reported assignment: A->ATP, B->ADP, C->AMP, D->Pi, E->H2O, F->H+.
// CHECK-NEXT: domain[0]: 1 ATP
// CHECK-NEXT: domain[1]: 1 ADP
// CHECK-NEXT: domain[2]: 1 AMP
// CHECK-NEXT: domain[3]: 1 Phosphate
// CHECK-NEXT: domain[4]: 1 H2O
// CHECK-NEXT: domain[5]: 1 H+

// Superset lookup reaches motifs carrying species beyond a reaction's support,
// so the candidate sets are wider than the eventual matches.
// CHECK-NEXT: candidates[0]: 2 rxnosc01 rxnosc03
// CHECK-NEXT: candidates[1]: 4 rxnosc01 rxnosc02 rxnosc03 rxnosc04
// CHECK-NEXT: candidates[2]: 1 rxnosc03
// CHECK-NEXT: candidates[3]: 2 rxnosc03 rxnosc04
