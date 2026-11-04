// Reversible isomerization A <=> B:
//   der(a) = -kf*a + kr*b
//   der(b) =  kf*a - kr*b
// 2 states, 4 terms, degree 1, kinetic.

// RUN: chemcomp-opt %s --ode-normalize --ode-dual-rail --convert-ode-to-acr \
// RUN:   --acr-network-merge | FileCheck %s
// RUN: chemcomp-opt %s --ode-to-crn=canonical-fallback=false \
// RUN:   | FileCheck %s --check-prefix=CRN
// RUN: chemcomp-opt %s --ode-to-crn | chemcomp-translate -mlir-to-json > %t.crn.json
// RUN: chemcomp-translate -ode-to-json %s > %t.ode.json
// RUN: %if emulator %{ %chemcomp_python %check_realization --ode %t.ode.json \
// RUN:   --crn %t.crn.json --t-end 4 | FileCheck %s --check-prefix=SIM %}

ode.system @Isomerization {
  ode.state @a initial 1.0
  ode.state @b initial 1.0
  ode.param @kf = 2.0
  ode.param @kr = 0.5
  ode.term @a += -1.0 params [@kf] states [@a]
  ode.term @a += 1.0 params [@kr] states [@b]
  ode.term @b += 1.0 params [@kf] states [@a]
  ode.term @b += -1.0 params [@kr] states [@b]
}

// CHECK-LABEL: acr.network @Isomerization {
// CHECK-NEXT:    acr.species @a initial 1.0
// CHECK-NEXT:    acr.species @b initial 1.0
// CHECK-NEXT:    acr.rxn @a -> @b rate 2.0 {from = "der(a): -kf*a; der(b): kf*a", rate_name = "kf", rule = "conversion"}
// CHECK-NEXT:    acr.rxn @b -> @a rate 0.5 {from = "der(b): -kr*b; der(a): kr*b", rate_name = "kr", rule = "conversion"}
// CHECK-NEXT:  }

// Composition only (no canonical fallback) gives the same reactions.
// CRN: acr.crn_generic
// CRN-SAME: rules = ["conversion", "conversion"]

// The emulator simulates the network; its states match the ODE solution.
// SIM: PASS
