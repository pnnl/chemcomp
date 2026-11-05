// Second-order decay: der(x) = -k*x^2. The full-match table of
// convert-bmodelica-to-acr has no network for it; a single rule does.
// 1 state, 1 term, degree 2, kinetic.

// RUN: chemcomp-opt %s --ode-normalize --ode-dual-rail --convert-ode-to-acr \
// RUN:   --acr-network-merge | FileCheck %s
// RUN: chemcomp-opt %s --ode-to-crn=canonical-fallback=false \
// RUN:   | FileCheck %s --check-prefix=CRN
// RUN: chemcomp-opt %s --ode-to-crn | chemcomp-translate -mlir-to-json > %t.crn.json
// RUN: chemcomp-translate -ode-to-json %s > %t.ode.json
// RUN: %if emulator %{ %chemcomp_python %check_realization --ode %t.ode.json \
// RUN:   --crn %t.crn.json --t-end 10 | FileCheck %s --check-prefix=SIM %}

ode.system @Dimerization {
  ode.state @x initial 2.0
  ode.param @k = 0.5
  ode.term @x += -1.0 params [@k] states [@x, @x]
}

// CHECK-LABEL: acr.network @Dimerization {
// CHECK-NEXT:    acr.species @x initial 2.0
// CHECK-NEXT:    acr.rxn 2 @x -> @x rate 0.5 {from = "der(x): -k*x^2", rate_name = "k", rule = "second-order-decay"}
// CHECK-NEXT:  }

// Composition only (no canonical fallback) gives the same reactions.
// CRN: acr.crn_generic
// CRN-SAME: rules = ["second-order-decay"]

// The emulator simulates the network; its states match the ODE solution.
// SIM: PASS
