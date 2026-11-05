// SIR epidemic:
//   der(s) = -beta*s*i
//   der(i) =  beta*s*i - gamma*i
//   der(r) =  gamma*i
// 3 states, 4 terms, degree 2, kinetic.

// RUN: chemcomp-opt %s --ode-normalize --ode-dual-rail --convert-ode-to-acr \
// RUN:   --acr-network-merge | FileCheck %s
// RUN: chemcomp-opt %s --ode-to-crn=canonical-fallback=false \
// RUN:   | FileCheck %s --check-prefix=CRN
// RUN: chemcomp-opt %s --ode-to-crn | chemcomp-translate -mlir-to-json > %t.crn.json
// RUN: chemcomp-translate -ode-to-json %s > %t.ode.json
// RUN: %if emulator %{ %chemcomp_python %check_realization --ode %t.ode.json \
// RUN:   --crn %t.crn.json --t-end 60 | FileCheck %s --check-prefix=SIM %}

ode.system @SIR {
  ode.state @s initial 0.99
  ode.state @i initial 0.01
  ode.state @r initial 0.0
  ode.param @beta = 0.5
  ode.param @gamma = 0.1
  ode.term @s += -1.0 params [@beta] states [@i, @s]
  ode.term @i += 1.0 params [@beta] states [@i, @s]
  ode.term @i += -1.0 params [@gamma] states [@i]
  ode.term @r += 1.0 params [@gamma] states [@i]
}

// CHECK-LABEL: acr.network @SIR {
// CHECK-NEXT:    acr.species @s initial 0.99
// CHECK-NEXT:    acr.species @i initial 0.01
// CHECK-NEXT:    acr.species @r initial 0.0
// CHECK-NEXT:    acr.rxn @i + @s -> 2 @i rate 0.5 {from = "der(s): -beta*i*s; der(i): beta*i*s", rate_name = "beta", rule = "predation"}
// CHECK-NEXT:    acr.rxn @i -> @r rate 0.1 {from = "der(i): -gamma*i; der(r): gamma*i", rate_name = "gamma", rule = "conversion"}
// CHECK-NEXT:  }

// Composition only (no canonical fallback) gives the same reactions.
// CRN: acr.crn_generic
// CRN-SAME: rules = ["predation", "conversion"]

// The emulator simulates the network; its states match the ODE solution.
// SIM: PASS
