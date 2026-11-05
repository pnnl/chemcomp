// RUN: chemcomp-translate -ode-to-json %s | FileCheck %s

// The reference ODE for checking compiled networks: parameter values are
// folded into the coefficients.

// CHECK:      "name": "Decay",
// CHECK:      "states": [
// CHECK-NEXT:   {
// CHECK-NEXT:     "name": "x",
// CHECK-NEXT:     "initial": 1
// CHECK:      "terms": [
// CHECK-NEXT:   {
// CHECK-NEXT:     "state": "x",
// CHECK-NEXT:     "coeff": -0.25,
// CHECK-NEXT:     "monomial": {
// CHECK-NEXT:       "x": 1
ode.system @Decay {
  ode.state @x initial 1.0
  ode.param @k = 0.5
  ode.term @x += -0.5 params [@k] states [@x]
}
