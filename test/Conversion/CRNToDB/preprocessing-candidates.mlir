// Candidate selection: signature lookup mode, the per-reaction cap, and the
// remark that makes truncation visible.

// RUN: chemcomp-opt %s -crn-db-match='dump-preprocessing=1' \
// RUN:   | FileCheck %s --check-prefix=SUPERSET
// RUN: chemcomp-opt %s -crn-db-match='strict-signature=1 dump-preprocessing=1' \
// RUN:   | FileCheck %s --check-prefix=STRICT
// RUN: chemcomp-opt %s \
// RUN:   -crn-db-match='db-slice=AmbiguousSlice max-candidates=2 dump-preprocessing=1' \
// RUN:   | FileCheck %s --check-prefix=CAPPED
// RUN: chemcomp-opt %s -crn-db-match='db-slice=AmbiguousSlice max-candidates=2' \
// RUN:   -verify-diagnostics
// RUN: chemcomp-opt %s \
// RUN:   -crn-db-match='db-slice=AmbiguousSlice max-per-signature=2 dump-preprocessing=1' \
// RUN:   | FileCheck %s --check-prefix=SLICED

acr.crn @one(%sv: tensor<4x!acr.species>, %rv: tensor<1xf32>)
    -> tensor<4x!acr.species> {
  // expected-remark @+1 {{candidate set for reaction 0 truncated from 5 to 2; the result is heuristic rather than optimal}}
  %out = acr.crn_generic
      {stoich_matrix = dense<[[-1, -1, 1, 1]]> : tensor<1x4xi32>}
      (%sv, %rv) : (tensor<4x!acr.species>, tensor<1xf32>)
                -> tensor<4x!acr.species>
  acr.return %out : tensor<4x!acr.species>
}

// Superset lookup reaches rxn00011, whose signature (-1,-1,1,1,1) properly
// contains the query's (-1,-1,1,1) -- the paper's coverage rule allows a motif
// to carry species beyond the query's support.
// SUPERSET: candidates[0]: 7 rxn00011 rxnamb01 rxnamb02 rxnamb03 rxnamb04 rxnamb05 rxnosc03

// Exact-signature lookup is paper-literal and misses it.
// STRICT: candidates[0]: 6 rxnamb01 rxnamb02 rxnamb03 rxnamb04 rxnamb05 rxnosc03

// The cap keeps the lowest motif ids, so truncation is reproducible, and the
// discarded count is reported rather than lost.
// CAPPED: candidates[0]: 2 rxnamb01 rxnamb02
// CAPPED: truncated[0]: 5 -> 2

// max-per-signature is a different mechanism: it shrinks the database at load
// time, so the discarded motifs are never candidates and nothing is reported
// as truncated. Keeping the two separate is what stops a load-time cut from
// silently suppressing the remark above.
// SLICED: candidates[0]: 2 rxnamb01 rxnamb02
// SLICED-NOT: truncated
