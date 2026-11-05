// RUN: chemcomp-opt %s --acr-network-merge -split-input-file | FileCheck %s

// Same reactants, same rate: X + Y -> Y and X + Y -> X + 2Y act together
// like X + Y -> 2Y (net change -X + Y at rate b*x*y either way).
// CHECK-LABEL: acr.network @Predation {
// CHECK-NEXT:    acr.species @x
// CHECK-NEXT:    acr.species @y
// CHECK-NEXT:    acr.rxn @x + @y -> 2 @y rate 0.4 {from = "der(x): -b*x*y; der(y): b*x*y", rate_name = "b", rule = "canonical"}
// CHECK-NEXT:  }
acr.network @Predation {
  acr.species @x
  acr.species @y
  acr.rxn @x + @y -> @y rate 0.4 {rate_name = "b", rule = "canonical", from = "der(x): -b*x*y"}
  acr.rxn @x + @y -> @x + 2 @y rate 0.4 {rate_name = "b", rule = "canonical", from = "der(y): b*x*y"}
}

// -----

// Identical reactions add their rates; a reaction that changes nothing is
// dropped.
// CHECK-LABEL: acr.network @Sum {
// CHECK-NEXT:    acr.species @x
// CHECK-NEXT:    acr.rxn @x -> none rate 0.75 {rate_name = "k1 + k2"}
// CHECK-NEXT:  }
acr.network @Sum {
  acr.species @x
  acr.rxn @x -> none rate 0.5 {rate_name = "k1"}
  acr.rxn @x -> none rate 0.25 {rate_name = "k2"}
}

// -----

// X -> 0 and X -> Y at the same rate would need X -> Y - X: not a reaction.
// Different rates are never merged.
// CHECK-LABEL: acr.network @Refused {
// CHECK-NEXT:    acr.species @x
// CHECK-NEXT:    acr.species @y
// CHECK-NEXT:    acr.rxn @x -> none rate 1.0 {rate_name = "k"}
// CHECK-NEXT:    acr.rxn @x -> @y rate 1.0 {rate_name = "k"}
// CHECK-NEXT:    acr.rxn @y -> none rate 1.0 {rate_name = "k"}
// CHECK-NEXT:    acr.rxn @y -> 2 @y rate 2.0 {rate_name = "k2"}
// CHECK-NEXT:  }
acr.network @Refused {
  acr.species @x
  acr.species @y
  acr.rxn @x -> none rate 1.0 {rate_name = "k"}
  acr.rxn @x -> @y rate 1.0 {rate_name = "k"}
  acr.rxn @y -> none rate 1.0 {rate_name = "k"}
  acr.rxn @y -> 2 @y rate 2.0 {rate_name = "k2"}
}
