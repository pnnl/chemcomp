// Exponential decay: der(x) = -k*x.
// 1 state, 1 term, degree 1, kinetic.

// RUN: chemcomp-opt %s --ode-normalize --ode-dual-rail --convert-ode-to-acr \
// RUN:   --acr-network-merge | FileCheck %s
// RUN: chemcomp-opt %s --ode-to-crn=canonical-fallback=false \
// RUN:   | FileCheck %s --check-prefix=CRN
// RUN: chemcomp-opt %s --ode-to-crn | chemcomp-translate -mlir-to-json > %t.crn.json
// RUN: chemcomp-translate -ode-to-json %s > %t.ode.json
// RUN: %if emulator %{ %chemcomp_python %check_realization --ode %t.ode.json \
// RUN:   --crn %t.crn.json --t-end 10 | FileCheck %s --check-prefix=SIM %}

ode.system @Decay {
  ode.state @x initial 1.0
  ode.param @k = 0.5
  ode.term @x += -1.0 params [@k] states [@x]
}

// CHECK-LABEL: acr.network @Decay {
// CHECK-NEXT:    acr.species @x initial 1.0
// CHECK-NEXT:    acr.rxn @x -> none rate 0.5 {from = "der(x): -k*x", rate_name = "k", rule = "degradation"}
// CHECK-NEXT:  }

// Composition only (no canonical fallback) gives the same reactions.
// CRN: acr.crn_generic
// CRN-SAME: rules = ["degradation"]

// The emulator simulates the network; its states match the ODE solution.
// SIM: PASS
