// RUN: chemcomp-tblgen -gen-motif-db -I %S/../../include %s | FileCheck %s

// Exercises the -gen-motif-db backend on a self-contained record set, rather
// than on the real fixture, so the expectations stay readable and stable when
// motifs are added to MotifsFixture.td.

include "Standalone/CRN2DB/MotifBase.td"

//===----------------------------------------------------------------------===//
// Counts, emitted ahead of the tables.
//===----------------------------------------------------------------------===//

// CHECK: kNumFixtureMotifs = 3;
// CHECK: kNumFixtureSlices = 2;

//===----------------------------------------------------------------------===//
// A fully populated record.
//===----------------------------------------------------------------------===//

def RxnFull : Motif<"rxn00002", ["H2O", "PPi", "Phosphate"], [-1, -1, 2]> {
  let DeltaG = "-3.46";
  let DeltaGErr = "0.05";
  let RateRange = ["0.01", "1.0"];
  let Source = "KEGG:R00004|ModelSEED:rxn00002";
}

//===----------------------------------------------------------------------===//
// Absent thermodynamic data and absent rate range must stay distinguishable
// from zero: 43.1% of ModelSEED rows carry no deltaG, and scoring a missing
// measurement as 0 would fabricate one.
//===----------------------------------------------------------------------===//

def RxnUnknown : Motif<"rxn00003", ["A", "B"], [-1, 1]> {
  let Source = "ModelSEED:rxn00003";
}

//===----------------------------------------------------------------------===//
// The Soluble property bit, and the sort-by-id normalization: this record is
// declared last but carries the lowest id, so it must be emitted first.
//===----------------------------------------------------------------------===//

def RxnSoluble : SolubleMotif<"rxn00001", ["H2O2", "H2O", "O2"], [-2, 2, 1]> {
  let DeltaG = "-46.06";
  let Source = "ModelSEED:rxn00001";
}

// One ordered block covers every row: the table is emitted in motif-id order
// regardless of declaration order, so pinning the sequence also proves the
// sort. The trailing 1 on rxn00001 is kMotifFlagSoluble (bit 0), set on that
// record alone.
// CHECK-LABEL: const FixtureMotif kFixtureMotifs[] = {
// CHECK-NEXT:  {"rxn00001", {"H2O2", "H2O", "O2"}, {-2, 2, 1}, -46.06, kUnknown, 0.0, 0.0, false, "ModelSEED:rxn00001", 1},
// CHECK-NEXT:  {"rxn00002", {"H2O", "PPi", "Phosphate"}, {-1, -1, 2}, -3.46, 0.05, 0.01, 1.0, true, "KEGG:R00004|ModelSEED:rxn00002", 0},
// CHECK-NEXT:  {"rxn00003", {"A", "B"}, {-1, 1}, kUnknown, kUnknown, 0.0, 0.0, false, "ModelSEED:rxn00003", 0},

//===----------------------------------------------------------------------===//
// Slices: emitted in name order, with membership sorted by motif id and
// de-duplicated.
//===----------------------------------------------------------------------===//

def ZSlice : Slice<"ZSlice", [RxnFull, RxnSoluble, RxnFull]> {
  let Summary = "Duplicate membership collapses.";
}

def ASlice : Slice<"ASlice", [RxnUnknown]>;

// CHECK-LABEL: const FixtureSlice kFixtureSlices[] = {
// CHECK-NEXT:  {"ASlice", {"rxn00003"}},
// CHECK-NEXT:  // Duplicate membership collapses.
// CHECK-NEXT:  {"ZSlice", {"rxn00001", "rxn00002"}},
