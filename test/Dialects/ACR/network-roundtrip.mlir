// RUN: chemcomp-opt %s | chemcomp-opt | FileCheck %s

// CHECK-LABEL: acr.network @LotkaVolterra {
// CHECK-NEXT:    acr.species @x initial 2.0
// CHECK-NEXT:    acr.species @y initial 0.5 {ode_state = "y"}
// CHECK-NEXT:    acr.rxn @x -> 2 @x rate 1.1 {rate_name = "a", rule = "autocatalysis"}
// CHECK-NEXT:    acr.rxn @x + @y -> 2 @y rate 0.4 {rate_name = "b", rule = "predation"}
// CHECK-NEXT:    acr.rxn @y -> none rate 0.4
// CHECK-NEXT:    acr.rxn none -> @x
// CHECK-NEXT:    acr.rxn 2 @x + @y -> 3 @x rate 1.0e-06
// CHECK-NEXT:  }
acr.network @LotkaVolterra {
  acr.species @x initial 2.0
  acr.species @y initial 0.5 {ode_state = "y"}
  acr.rxn @x -> 2 @x rate 1.1 {rate_name = "a", rule = "autocatalysis"}
  acr.rxn @y + @x -> 2 @y rate 0.4 {rate_name = "b", rule = "predation"}
  acr.rxn @y -> none rate 0.4
  acr.rxn none -> @x
  acr.rxn 2 @x + @y -> 3 @x rate 0.000001
}
