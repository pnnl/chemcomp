// Slice selection and the SliceSpec filters.

// RUN: chemcomp-opt %s -crn-db-match='db-slice=OscillatorSlice dump-db=1' \
// RUN:   | FileCheck %s --check-prefix=OSC
// RUN: chemcomp-opt %s -crn-db-match='db-slice=AmbiguousSlice dump-db=1' \
// RUN:   | FileCheck %s --check-prefix=AMB
// RUN: chemcomp-opt %s -crn-db-match='dump-db=1' \
// RUN:   | FileCheck %s --check-prefix=ALL
// RUN: chemcomp-opt %s \
// RUN:   -crn-db-match='db-slice=AmbiguousSlice max-per-signature=2 dump-db=1' \
// RUN:   | FileCheck %s --check-prefix=CAPPED
// RUN: chemcomp-opt %s \
// RUN:   -crn-db-match='filter-species=ATP dump-db=1' \
// RUN:   | FileCheck %s --check-prefix=SPECIES
// RUN: chemcomp-opt %s \
// RUN:   -crn-db-match='filter-source=KEGG: dump-db=1' \
// RUN:   | FileCheck %s --check-prefix=SOURCE

module {}

// The four motifs backing the paper's section 5.4 walkthrough.
// OSC: motifs: 4
// OSC: motif rxnosc01
// OSC-NEXT: motif rxnosc02
// OSC-NEXT: motif rxnosc03
// OSC-NEXT: motif rxnosc04

// Every ambiguity-set motif shares one signature -- the shape that holds
// 11,421 of 43,455 motifs in the real database.
// AMB: motifs: 5
// AMB: signatures: 1
// AMB: sig=(-1, -1, 1, 1)

// An empty db-slice loads the whole fixture.
// ALL: motifs: 22

// Capping keeps the lowest motif ids, so truncation is reproducible.
// CAPPED: motifs: 2
// CAPPED: motif rxnamb01
// CAPPED-NEXT: motif rxnamb02
// CAPPED-NOT: motif rxnamb03

// Species filtering keeps only motifs mentioning ATP.
// SPECIES: motifs: 3
// SPECIES: motif rxnamb01
// SPECIES-NEXT: motif rxnosc01
// SPECIES-NEXT: motif rxnosc03

// Source filtering by prefix. Only the paper's Table 1 records carry a KEGG
// cross-reference; the CSVs have no such column, so this dimension is nearly
// single-valued in practice.
// SOURCE: motifs: 10
// SOURCE: source=KEGG:R00004|ModelSEED:rxn00001
