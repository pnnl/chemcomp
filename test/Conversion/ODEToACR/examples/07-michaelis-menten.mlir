// Michaelis-Menten enzyme kinetics with full mass action, E + S <=> C -> E + P:
//   der(e) = -k1*e*s + km*c + kc*c
//   der(s) = -k1*e*s + km*c
//   der(c) =  k1*e*s - km*c - kc*c
//   der(p) =  kc*c
// 4 states, 9 terms, degree 2, kinetic.

// RUN: chemcomp-opt %s --ode-normalize --ode-dual-rail --convert-ode-to-acr \
// RUN:   --acr-network-merge | FileCheck %s
// RUN: chemcomp-opt %s --ode-to-crn=canonical-fallback=false \
// RUN:   | FileCheck %s --check-prefix=CRN
// RUN: chemcomp-opt %s --ode-to-crn | chemcomp-translate -mlir-to-json > %t.crn.json
// RUN: chemcomp-translate -ode-to-json %s > %t.ode.json
// RUN: %if emulator %{ %chemcomp_python %check_realization --ode %t.ode.json \
// RUN:   --crn %t.crn.json --t-end 20 | FileCheck %s --check-prefix=SIM %}

ode.system @MichaelisMenten {
  ode.state @e initial 1.0
  ode.state @s initial 2.0
  ode.state @c initial 0.0
  ode.state @p initial 0.0
  ode.param @k1 = 1.0
  ode.param @km = 0.5
  ode.param @kc = 0.3
  ode.term @e += -1.0 params [@k1] states [@e, @s]
  ode.term @e += 1.0 params [@km] states [@c]
  ode.term @e += 1.0 params [@kc] states [@c]
  ode.term @s += -1.0 params [@k1] states [@e, @s]
  ode.term @s += 1.0 params [@km] states [@c]
  ode.term @c += 1.0 params [@k1] states [@e, @s]
  ode.term @c += -1.0 params [@km] states [@c]
  ode.term @c += -1.0 params [@kc] states [@c]
  ode.term @p += 1.0 params [@kc] states [@c]
}

// CHECK-LABEL: acr.network @MichaelisMenten {
// CHECK-NEXT:    acr.species @e initial 1.0
// CHECK-NEXT:    acr.species @s initial 2.0
// CHECK-NEXT:    acr.species @c initial 0.0
// CHECK-NEXT:    acr.species @p initial 0.0
// CHECK-NEXT:    acr.rxn @e + @s -> @c rate 1.0 {from = "der(e): -k1*e*s; der(s): -k1*e*s; der(c): k1*e*s", rate_name = "k1", rule = "binding"}
// CHECK-NEXT:    acr.rxn @c -> @e + @s rate 0.5 {from = "der(c): -km*c; der(e): km*c; der(s): km*c", rate_name = "km", rule = "dissociation"}
// CHECK-NEXT:    acr.rxn @c -> @e + @p rate 0.3 {from = "der(c): -kc*c; der(e): kc*c; der(p): kc*c", rate_name = "kc", rule = "dissociation"}
// CHECK-NEXT:  }

// Composition only (no canonical fallback) gives the same reactions.
// CRN: acr.crn_generic
// CRN-SAME: rules = ["binding", "dissociation", "dissociation"]

// The emulator simulates the network; its states match the ODE solution.
// SIM: PASS
