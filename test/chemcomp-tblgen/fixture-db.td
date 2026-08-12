// RUN: chemcomp-tblgen -gen-motif-db -I %S/../../include \
// RUN:   %S/../../include/Standalone/CRN2DB/MotifsFixture.td | FileCheck %s

// Guards the committed fixture itself. Deliberately checks structure and a few
// anchor records rather than every row, so adding motifs does not churn the
// test, but a schema change or a corrupted record does fail it.

// The paper's Table 1 records, transcribed from the source CSV. Note rxn00003
// is stored reactant-first: its equation in reactions_reduced.csv uses a `<=`
// arrow, which the paper transcribes without flipping.
// CHECK-DAG: {"rxn00001", {"H2O", "PPi", "Phosphate", "H+"}, {-1, -1, 2, 1}, -3.46,
// CHECK-DAG: {"rxn00003", {"Pyruvate", "H+", "CO2", "ALCTT"}, {-2, -1, 1, 1}, 8.27,
// CHECK-DAG: {"rxn00006", {"H2O2", "H2O", "O2"}, {-2, 2, 1}, -46.06,

// rxn00005 is genuinely absent from ModelSEED; the numbering gap is real.
// CHECK-NOT: "rxn00005"

// The oscillator motifs backing the paper's section 5.4 walkthrough.
// CHECK-DAG: {"rxnosc01", {"ATP", "ADP", "AMP"}, {-1, -1, 1},
// CHECK-DAG: {"rxnosc04", {"H+", "ADP", "H2O"}, {-1, 1, 1},

// rxnamb05 carries no thermodynamic data at all.
// CHECK-DAG: {"rxnamb05", {"NAD", "NADH", "FAD", "FADH2"}, {-1, -1, 1, 1}, kUnknown, kUnknown, 0.0, 0.0, false,

// Coverage probes: rxnsub02 is rxnsub01 plus one species, so it is reachable
// only through superset signature lookup.
// CHECK-DAG: {"rxnsub01", {"Glucose", "Fructose"}, {-1, 1},
// CHECK-DAG: {"rxnsub02", {"Glucose", "Fructose", "H2O"}, {-1, 1, 1},

// Slices are emitted in name order.
// CHECK-LABEL: const FixtureSlice kFixtureSlices[] = {
// CHECK: {"AllFixture",
// CHECK: {"AmbiguousSlice", {"rxnamb01", "rxnamb02", "rxnamb03", "rxnamb04", "rxnamb05"}},
// CHECK: {"OscillatorSlice", {"rxnosc01", "rxnosc02", "rxnosc03", "rxnosc04"}},
// CHECK: {"PaperTable1",
// CHECK: {"SolubleSlice", {"rxn00006", "rxnosc01", "rxnosc02", "rxnosc03", "rxnosc04"}},
// CHECK: {"SubsetSlice", {"rxnsub01", "rxnsub02", "rxnsub03"}},
// CHECK: {"TinySlice", {"rxn00001", "rxn00004", "rxn00006"}},
