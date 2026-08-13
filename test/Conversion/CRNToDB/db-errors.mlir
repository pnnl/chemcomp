// Loader failures must be diagnosed, not crash and not silently yield an
// empty database -- an empty database looks like "nothing matched", which
// sends you debugging the wrong thing.

// RUN: not chemcomp-opt %s -crn-db-match='db=%S/Inputs/does-not-exist.json' 2>&1 \
// RUN:   | FileCheck %s --check-prefix=MISSING
// RUN: not chemcomp-opt %s -crn-db-match='db=%S/Inputs/malformed-db.json' 2>&1 \
// RUN:   | FileCheck %s --check-prefix=MALFORMED
// RUN: not chemcomp-opt %s -crn-db-match='db=%S/Inputs/mismatched-db.json' 2>&1 \
// RUN:   | FileCheck %s --check-prefix=MISMATCH
// RUN: not chemcomp-opt %s -crn-db-match='db-slice=NoSuchSlice' 2>&1 \
// RUN:   | FileCheck %s --check-prefix=BADSLICE
// RUN: not chemcomp-opt %s \
// RUN:   -crn-db-match='db=%S/Inputs/tiny-db.json db-slice=TinySlice' 2>&1 \
// RUN:   | FileCheck %s --check-prefix=BOTH

module {}

// MISSING: could not open motif database
// MISSING-SAME: does-not-exist.json

// MALFORMED: expected a top-level array of motif objects

// A species/stoich length mismatch names the offending motif, which matters in
// a 43k-record file.
// MISMATCH: motif 'rxnbad01' has 3 species but 2 stoichiometric coefficients

// An unknown slice lists what is actually available rather than failing bare.
// BADSLICE: unknown fixture slice 'NoSuchSlice'
// BADSLICE-SAME: available slices: AllFixture, AmbiguousSlice, OscillatorSlice, PaperTable1, SolubleSlice, SubsetSlice, TinySlice

// db-slice names a fixture slice, so pairing it with a JSON database is a
// mistaken invocation rather than something to silently ignore.
// BOTH: 'db-slice' selects a slice of the compiled-in fixture and cannot be combined with 'db'
