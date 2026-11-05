// Brusselator (limit cycle for b > 1 + a^2):
//   der(x) = a - (b + 1)*x + x^2*y
//   der(y) = b*x - x^2*y
// 2 states, 6 terms, degrees 0 to 3, kinetic.

// RUN: chemcomp-opt %s --ode-normalize --ode-dual-rail --convert-ode-to-acr \
// RUN:   --acr-network-merge | FileCheck %s
// RUN: chemcomp-opt %s --ode-to-crn=canonical-fallback=false \
// RUN:   | FileCheck %s --check-prefix=CRN
// RUN: chemcomp-opt %s --ode-to-crn | chemcomp-translate -mlir-to-json > %t.crn.json
// RUN: chemcomp-translate -ode-to-json %s > %t.ode.json
// RUN: %if emulator %{ %chemcomp_python %check_realization --ode %t.ode.json \
// RUN:   --crn %t.crn.json --t-end 20 | FileCheck %s --check-prefix=SIM %}

ode.system @Brusselator {
  ode.state @x initial 1.0
  ode.state @y initial 1.0
  ode.param @a = 1.0
  ode.param @b = 3.0
  ode.term @x += 1.0 params [@a]
  ode.term @x += -1.0 params [@b] states [@x]
  ode.term @x += -1.0 states [@x]
  ode.term @x += 1.0 states [@x, @x, @y]
  ode.term @y += 1.0 params [@b] states [@x]
  ode.term @y += -1.0 states [@x, @x, @y]
}

// CHECK-LABEL: acr.network @Brusselator {
// CHECK-NEXT:    acr.species @x initial 1.0
// CHECK-NEXT:    acr.species @y initial 1.0
// CHECK-NEXT:    acr.rxn 2 @x + @y -> 3 @x rate 1.0 {from = "der(x): x^2*y; der(y): -x^2*y", rate_name = "1", rule = "cubic-autocatalysis"}
// CHECK-NEXT:    acr.rxn @x -> @y rate 3.0 {from = "der(x): -b*x; der(y): b*x", rate_name = "b", rule = "conversion"}
// CHECK-NEXT:    acr.rxn @x -> none rate 1.0 {from = "der(x): -x", rate_name = "1", rule = "degradation"}
// CHECK-NEXT:    acr.rxn none -> @x rate 1.0 {from = "der(x): a", rate_name = "a", rule = "production"}
// CHECK-NEXT:  }

// Composition only (no canonical fallback) gives the same reactions.
// CRN: acr.crn_generic
// CRN-SAME: rules = ["cubic-autocatalysis", "conversion", "degradation", "production"]

// The emulator simulates the network; its states match the ODE solution.
// SIM: PASS
