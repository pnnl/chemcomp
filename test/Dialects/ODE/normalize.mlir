// RUN: chemcomp-opt %s --ode-normalize | FileCheck %s

// Like terms are summed, cancelling terms are dropped, and terms are ordered
// by state (declaration order), then by monomial.

// CHECK-LABEL: ode.system @S {
// CHECK:         ode.term @x += -1.5 params [@k] states [@x]
// CHECK-NEXT:    ode.term @x += 1.0 states [@y]
// CHECK-NEXT:    ode.term @y += 3.0
// CHECK-NOT:     ode.term
// CHECK:       }
ode.system @S {
  ode.state @x
  ode.state @y
  ode.param @k = 1.0
  ode.term @y += 1.0
  ode.term @x += 1.0 states [@y]
  ode.term @x += -1.0 params [@k] states [@x]
  ode.term @y += 2.0
  ode.term @x += -0.5 params [@k] states [@x]
  ode.term @y += 0.25 states [@x, @y]
  ode.term @y += -0.25 states [@y, @x]
}
