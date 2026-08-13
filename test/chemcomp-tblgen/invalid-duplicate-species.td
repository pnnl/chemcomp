// RUN: not chemcomp-tblgen -gen-motif-db -I %S/../../include %s 2>&1 | FileCheck %s

// A species repeated inside one motif leaves nu_j(s) with two values.

include "Standalone/Conversion/CRNToDB/MotifBase.td"

def RxnRepeated : Motif<"rxn00001", ["H2O", "PPi", "H2O"], [-1, -1, 2]>;

// CHECK: error:
// CHECK-SAME: motif 'rxn00001' repeats species 'H2O'
