// RUN: chemcomp-opt %s | chemcomp-opt | FileCheck %s

// Lotka-Volterra: der(x) = a*x - b*x*y, der(y) = b*x*y - c*y.

// CHECK-LABEL: ode.system @LotkaVolterra {
// CHECK-NEXT:    ode.state @x initial 2.0
// CHECK-NEXT:    ode.state @y initial 0.5
// CHECK-NEXT:    ode.param @a = 1.1
// CHECK-NEXT:    ode.param @b = 0.4
// CHECK-NEXT:    ode.param @c
// CHECK-NEXT:    ode.term @x += 1.0 params [@a] states [@x]
// CHECK-NEXT:    ode.term @x += -1.0 params [@b] states [@x, @y]
// CHECK-NEXT:    ode.term @y += 1.0 params [@b] states [@x, @y]
// CHECK-NEXT:    ode.term @y += -1.0 params [@c] states [@y]
// CHECK-NEXT:  }
ode.system @LotkaVolterra {
  ode.state @x initial 2.0
  ode.state @y initial 0.5
  ode.param @a = 1.1
  ode.param @b = 0.4
  ode.param @c
  ode.term @x += 1.0 params [@a] states [@x]
  ode.term @x += -1.0 params [@b] states [@x, @y]
  ode.term @y += 1.0 params [@b] states [@x, @y]
  ode.term @y += -1.0 params [@c] states [@y]
}

// Constant terms, powers (repeated symbols), small numbers and rail
// attributes.
// CHECK-LABEL: ode.system @Misc attributes {dual_rail = ["x"]} {
// CHECK-NEXT:    ode.state @x_p {rail_of = "x", rail_sign = 1 : i64}
// CHECK-NEXT:    ode.state @x_m {rail_of = "x", rail_sign = -1 : i64}
// CHECK-NEXT:    ode.term @x_p += 2.0
// CHECK-NEXT:    ode.term @x_p += -1.0e-06 states [@x_p, @x_p, @x_m]
ode.system @Misc attributes {dual_rail = ["x"]} {
  ode.state @x_p {rail_of = "x", rail_sign = 1 : i64}
  ode.state @x_m {rail_of = "x", rail_sign = -1 : i64}
  ode.term @x_p += 2.0
  ode.term @x_p += -0.000001 states [@x_p, @x_p, @x_m]
}
