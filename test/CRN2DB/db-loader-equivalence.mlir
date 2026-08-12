// Both loaders must produce indistinguishable databases: the compiled-in
// TableGen fixture and the runtime JSON path feed the same MotifDatabase, and
// the matcher must not be able to tell which one it got.
//
// Inputs/tiny-db.json mirrors the fixture's TinySlice record for record, so
// the two dumps below are compared against one identical set of expectations.

// RUN: standalone-opt %s -crn-db-match='db-slice=TinySlice dump-db=1' | FileCheck %s
// RUN: standalone-opt %s -crn-db-match='db=%S/Inputs/tiny-db.json dump-db=1' | FileCheck %s

module {}

// CHECK: motifs: 3
// CHECK: species: 8
// CHECK: signatures: 3

// Motifs are emitted in id order, with species sorted by their interned id and
// the stoichiometry map printed alongside.
// CHECK-NEXT: motif rxn00001 sig=(-1, -1, 1, 2) {H2O:-1, PPi:-1, Phosphate:2, H+:1} deltaG=-3.46 rate=[0.01,1] source=KEGG:R00004|ModelSEED:rxn00001

// No rate_range: kappa_j stays unknown rather than defaulting to zero.
// CHECK-NEXT: motif rxn00004 sig=(-1, 2) {Parapyruvate:-1, Pyruvate:2} deltaG=4.49 rate=unknown source=KEGG:R00008|ModelSEED:rxn00004

// CHECK-NEXT: motif rxn00006 sig=(-2, 1, 2) {H2O:2, H2O2:-2, O2:1} deltaG=-46.06 rate=[0.1,10] source=KEGG:R00009|ModelSEED:rxn00006
