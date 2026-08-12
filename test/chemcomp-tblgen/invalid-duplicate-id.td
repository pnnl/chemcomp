// RUN: not chemcomp-tblgen -gen-motif-db -I %S/../../include %s 2>&1 | FileCheck %s

// Two records claiming one motif id would make slice membership ambiguous and
// silently shadow one another in the generated table.

include "Standalone/CRN2DB/MotifBase.td"

def RxnFirst : Motif<"rxn00001", ["A", "B"], [-1, 1]>;
def RxnSecond : Motif<"rxn00001", ["C", "D"], [-1, 1]>;

// CHECK: error:
// CHECK-SAME: duplicate motif id 'rxn00001'
