// Lotka-Volterra predator-prey:
//   der(x) = a*x - b*x*y      (x: prey)
//   der(y) = b*x*y - c*y      (y: predator)
// 2 states, 4 terms, degree 2, kinetic.

// RUN: chemcomp-opt %s --ode-normalize --ode-dual-rail --convert-ode-to-acr \
// RUN:   --acr-network-merge | FileCheck %s
// RUN: chemcomp-opt %s --ode-to-crn=canonical-fallback=false \
// RUN:   | FileCheck %s --check-prefix=CRN
// RUN: chemcomp-opt %s --ode-to-crn | chemcomp-translate -mlir-to-json > %t.crn.json
// RUN: chemcomp-translate -ode-to-json %s > %t.ode.json
// RUN: %if emulator %{ %chemcomp_python %check_realization --ode %t.ode.json \
// RUN:   --crn %t.crn.json --t-end 30 | FileCheck %s --check-prefix=SIM %}

ode.system @LotkaVolterra {
  ode.state @x initial 2.0
  ode.state @y initial 2.0
  ode.param @a = 1.1
  ode.param @b = 0.4
  ode.param @c = 0.4
  ode.term @x += 1.0 params [@a] states [@x]
  ode.term @x += -1.0 params [@b] states [@x, @y]
  ode.term @y += 1.0 params [@b] states [@x, @y]
  ode.term @y += -1.0 params [@c] states [@y]
}

// CHECK-LABEL: acr.network @LotkaVolterra {
// CHECK-NEXT:    acr.species @x initial 2.0
// CHECK-NEXT:    acr.species @y initial 2.0
// CHECK-NEXT:    acr.rxn @x + @y -> 2 @y rate 0.4 {from = "der(x): -b*x*y; der(y): b*x*y", rate_name = "b", rule = "predation"}
// CHECK-NEXT:    acr.rxn @x -> 2 @x rate 1.1 {from = "der(x): a*x", rate_name = "a", rule = "autocatalysis"}
// CHECK-NEXT:    acr.rxn @y -> none rate 0.4 {from = "der(y): -c*y", rate_name = "c", rule = "degradation"}
// CHECK-NEXT:  }

// Composition only (no canonical fallback) gives the same reactions.
// CRN: acr.crn_generic
// CRN-SAME: rules = ["predation", "autocatalysis", "degradation"]

// The emulator simulates the network; its states match the ODE solution.
// SIM: PASS
