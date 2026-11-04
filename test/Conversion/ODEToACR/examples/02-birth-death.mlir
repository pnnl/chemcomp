// Birth-death: der(x) = kb - kd*x.
// 1 state, 2 terms, degrees 0 and 1, kinetic.

// RUN: chemcomp-opt %s --ode-normalize --ode-dual-rail --convert-ode-to-acr \
// RUN:   --acr-network-merge | FileCheck %s
// RUN: chemcomp-opt %s --ode-to-crn=canonical-fallback=false \
// RUN:   | FileCheck %s --check-prefix=CRN
// RUN: chemcomp-opt %s --ode-to-crn | chemcomp-translate -mlir-to-json > %t.crn.json
// RUN: chemcomp-translate -ode-to-json %s > %t.ode.json
// RUN: %if emulator %{ %chemcomp_python %check_realization --ode %t.ode.json \
// RUN:   --crn %t.crn.json --t-end 10 | FileCheck %s --check-prefix=SIM %}

ode.system @BirthDeath {
  ode.state @x initial 0.0
  ode.param @kb = 2.0
  ode.param @kd = 0.5
  ode.term @x += 1.0 params [@kb]
  ode.term @x += -1.0 params [@kd] states [@x]
}

// CHECK-LABEL: acr.network @BirthDeath {
// CHECK-NEXT:    acr.species @x initial 0.0
// CHECK-NEXT:    acr.rxn @x -> none rate 0.5 {from = "der(x): -kd*x", rate_name = "kd", rule = "degradation"}
// CHECK-NEXT:    acr.rxn none -> @x rate 2.0 {from = "der(x): kb", rate_name = "kb", rule = "production"}
// CHECK-NEXT:  }

// Composition only (no canonical fallback) gives the same reactions.
// CRN: acr.crn_generic
// CRN-SAME: rules = ["degradation", "production"]

// The emulator simulates the network; its states match the ODE solution.
// SIM: PASS
