// RUN: not chemcomp-tblgen -gen-motif-db -I %S/../../include %s 2>&1 | FileCheck %s

// Species and Stoich are positionally parallel, so a length mismatch makes the
// stoichiometry map ambiguous. TableGen cannot express the constraint, so the
// backend enforces it -- with a source location, which is why this is caught
// here rather than by MotifDatabase::addMotif returning false at load time.

include "Standalone/CRN2DB/MotifBase.td"

def RxnMismatched : Motif<"rxn00001", ["H2O", "PPi", "Phosphate"], [-1, 1]>;

// CHECK: error:
// CHECK-SAME: motif 'rxn00001' has 3 species but 2 stoichiometric coefficients
