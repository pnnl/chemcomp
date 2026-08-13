// The headline conformance test: the paper's section 5.4 walkthrough.
//
//     A + B -> C      C -> D      D + E -> A + F      F -> B + E
//
// The paper reports the assignment A->ATP, B->ADP, C->AMP, D->Pi, E->H2O,
// F->H+ with full coverage, and the dynamic program must reproduce it.
//
// The paper is candid that its result is "a database-consistent assignment
// achieving full coverage, though not enforcing stoichiometric or
// thermodynamic feasibility"; the OscillatorSlice motifs inherit that caveat.

// RUN: standalone-opt %s -crn-db-match='db-slice=OscillatorSlice' | FileCheck %s

// Full coverage, so no diagnostics at all.
// RUN: standalone-opt %s -crn-db-match='db-slice=OscillatorSlice' \
// RUN:   -verify-diagnostics

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

// CHECK: acr.crn_generic
// CHECK-SAME: motif_ids = ["rxnosc01", "rxnosc02", "rxnosc03", "rxnosc04"]
// CHECK-SAME: species_names = ["ATP", "ADP", "AMP", "Phosphate", "H2O", "H+"]
