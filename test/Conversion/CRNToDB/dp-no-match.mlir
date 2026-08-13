// Partial and total match failure.
//
// Section 5.1 specifies a warning listing the unmatched reactions, after which
// the pass proceeds with the best partial matching so downstream passes can
// still consume a partially instantiated CRN. Unmatched entries are the empty
// string rather than absent, keeping the arrays positionally aligned with the
// stoichiometric matrix.

// RUN: chemcomp-opt %s -crn-db-match='db-slice=OscillatorSlice' -verify-diagnostics \
// RUN:   | FileCheck %s

// Row 0 has signature (-1,-1,1) and matches rxnosc01. Row 1 needs a species
// with coefficient -3, which no fixture motif offers, so it cannot match at
// all. The rows share no species, so row 1's failure is its own rather than a
// consequence of row 0's binding.
acr.crn @partial(%sv: tensor<5x!acr.species>, %rv: tensor<2xf32>)
    -> tensor<5x!acr.species> {
  // expected-warning @+1 {{unmatched reactions: 1}}
  %out = acr.crn_generic
      {stoich_matrix = dense<[[-1, -1,  1,  0,  0],
                              [ 0,  0,  0, -3,  1]]> : tensor<2x5xi32>}
      (%sv, %rv) : (tensor<5x!acr.species>, tensor<2xf32>)
                -> tensor<5x!acr.species>
  acr.return %out : tensor<5x!acr.species>
}

// The matched row keeps its motif and bindings; the unmatched one leaves empty
// strings, so both arrays stay aligned with the stoichiometric matrix.
// CHECK: acr.crn_generic
// CHECK-SAME: motif_ids = ["rxnosc01", ""]
// CHECK-SAME: species_names = ["ATP", "ADP", "AMP", "", ""]
