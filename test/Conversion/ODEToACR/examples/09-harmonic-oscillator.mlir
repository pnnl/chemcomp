// Harmonic oscillator:
//   der(x) = -w*y
//   der(y) =  w*x
// 2 states, 2 terms, degree 1, NOT kinetic: -w*y in der(x) does not contain x
// (and x, y go negative). Dual-rail encoding makes it compile.

// RUN: chemcomp-opt %s --ode-normalize --ode-dual-rail --convert-ode-to-acr \
// RUN:   --acr-network-merge | FileCheck %s
// RUN: chemcomp-opt %s --ode-to-crn=canonical-fallback=false \
// RUN:   | FileCheck %s --check-prefix=CRN
// RUN: chemcomp-opt %s --ode-to-crn | chemcomp-translate -mlir-to-json > %t.crn.json
// RUN: chemcomp-translate -ode-to-json %s > %t.ode.json
// RUN: %if emulator %{ %chemcomp_python %check_realization --ode %t.ode.json \
// RUN:   --crn %t.crn.json --t-end 10 | FileCheck %s --check-prefix=SIM %}

ode.system @Oscillator {
  ode.state @x initial 1.0
  ode.state @y initial 0.0
  ode.param @w = 1.0
  ode.term @x += -1.0 params [@w] states [@y]
  ode.term @y += 1.0 params [@w] states [@x]
}

// CHECK-LABEL: acr.network @Oscillator {
// CHECK-NEXT:    acr.species @x_p initial 1.0 {ode_state = "x", rail_sign = 1 : i64}
// CHECK-NEXT:    acr.species @x_m initial 0.0 {ode_state = "x", rail_sign = -1 : i64}
// CHECK-NEXT:    acr.species @y_p initial 0.0 {ode_state = "y", rail_sign = 1 : i64}
// CHECK-NEXT:    acr.species @y_m initial 0.0 {ode_state = "y", rail_sign = -1 : i64}
// CHECK-NEXT:    acr.rxn @x_m + @x_p -> none rate 1.0 {from = "der(x_p): -x_m*x_p; der(x_m): -x_m*x_p", rate_name = "1", rule = "annihilation"}
// CHECK-NEXT:    acr.rxn @y_m + @y_p -> none rate 1.0 {from = "der(y_p): -y_m*y_p; der(y_m): -y_m*y_p", rate_name = "1", rule = "annihilation"}
// CHECK-NEXT:    acr.rxn @x_p -> @x_p + @y_p rate 1.0 {from = "der(y_p): w*x_p", rate_name = "w", rule = "catalytic-production"}
// CHECK-NEXT:    acr.rxn @x_m -> @x_m + @y_m rate 1.0 {from = "der(y_m): w*x_m", rate_name = "w", rule = "catalytic-production"}
// CHECK-NEXT:    acr.rxn @y_p -> @x_m + @y_p rate 1.0 {from = "der(x_m): w*y_p", rate_name = "w", rule = "catalytic-production"}
// CHECK-NEXT:    acr.rxn @y_m -> @x_p + @y_m rate 1.0 {from = "der(x_p): w*y_m", rate_name = "w", rule = "catalytic-production"}
// CHECK-NEXT:  }

// Composition only (no canonical fallback) gives the same reactions.
// CRN: acr.crn_generic
// CRN-SAME: rules = ["annihilation", "annihilation", "catalytic-production", "catalytic-production", "catalytic-production", "catalytic-production"]

// The emulator simulates the network; its states match the ODE solution.
// SIM: PASS
