// RUN: chemcomp-opt %s --ode-dual-rail -split-input-file | FileCheck %s
// RUN: chemcomp-opt %s --ode-dual-rail=mode=all -split-input-file \
// RUN:   | FileCheck %s --check-prefix=ALL
// RUN: chemcomp-opt %s --ode-dual-rail=annihilation-rate=5 -split-input-file \
// RUN:   | FileCheck %s --check-prefix=ETA

// Kinetic systems are left alone in auto mode.
// CHECK-LABEL: ode.system @Decay {
// CHECK-NEXT:    ode.state @x initial 1.0
// CHECK-NEXT:    ode.term @x += -1.0 states [@x]
// CHECK-NEXT:  }

// With mode=all, -x stays a degradation of the plus rail and the minus rail
// only degrades: x_m starts at 0 and stays there.
// ALL-LABEL: ode.system @Decay attributes {dual_rail = ["x"]} {
// ALL-NEXT:    ode.state @x_p initial 1.0 {rail_of = "x", rail_sign = 1 : i64}
// ALL-NEXT:    ode.state @x_m initial 0.0 {rail_of = "x", rail_sign = -1 : i64}
// ALL-NEXT:    ode.term @x_p += -1.0 states [@x_m, @x_p]
// ALL-NEXT:    ode.term @x_p += -1.0 states [@x_p]
// ALL-NEXT:    ode.term @x_m += -1.0 states [@x_m]
// ALL-NEXT:    ode.term @x_m += -1.0 states [@x_m, @x_p]
// ALL-NEXT:  }
ode.system @Decay {
  ode.state @x initial 1.0
  ode.term @x += -1.0 states [@x]
}

// -----

// Harmonic oscillator der(x) = -w*y, der(y) = w*x. Only y has a cross-effect
// at first; splitting it turns -w*y into -w*y_p + w*y_m in der(x), so x gets
// split too. Every rail is produced by another rail and the pairs annihilate.
// x = -1 starts on the minus rail.
// CHECK-LABEL: ode.system @Oscillator attributes {dual_rail = ["x", "y"]} {
// CHECK-NEXT:    ode.state @x_p initial 0.0 {rail_of = "x", rail_sign = 1 : i64}
// CHECK-NEXT:    ode.state @x_m initial 1.0 {rail_of = "x", rail_sign = -1 : i64}
// CHECK-NEXT:    ode.state @y_p initial 0.5 {rail_of = "y", rail_sign = 1 : i64}
// CHECK-NEXT:    ode.state @y_m initial 0.0 {rail_of = "y", rail_sign = -1 : i64}
// CHECK-NEXT:    ode.param @w = 2.0
// CHECK-NEXT:    ode.term @x_p += -1.0 states [@x_m, @x_p]
// CHECK-NEXT:    ode.term @x_p += 1.0 params [@w] states [@y_m]
// CHECK-NEXT:    ode.term @x_m += -1.0 states [@x_m, @x_p]
// CHECK-NEXT:    ode.term @x_m += 1.0 params [@w] states [@y_p]
// CHECK-NEXT:    ode.term @y_p += 1.0 params [@w] states [@x_p]
// CHECK-NEXT:    ode.term @y_p += -1.0 states [@y_m, @y_p]
// CHECK-NEXT:    ode.term @y_m += 1.0 params [@w] states [@x_m]
// CHECK-NEXT:    ode.term @y_m += -1.0 states [@y_m, @y_p]
// CHECK-NEXT:  }

// ETA-LABEL: ode.system @Oscillator
// ETA:         ode.term @x_p += -5.0 states [@x_m, @x_p]
ode.system @Oscillator {
  ode.state @x initial -1.0
  ode.state @y initial 0.5
  ode.param @w = 2.0
  ode.term @x += -1.0 params [@w] states [@y]
  ode.term @y += 1.0 params [@w] states [@x]
}

// -----

// Squares expand binomially: der(z) = -x^2 with x split gives
// -(x_p^2 - 2*x_p*x_m + x_m^2), routed term by term.
// CHECK-LABEL: ode.system @Square attributes {dual_rail = ["x", "z"]} {
// CHECK:         ode.term @z_p += 2.0 states [@x_m, @x_p]
// CHECK:         ode.term @z_m += 1.0 states [@x_m, @x_m]
// CHECK:         ode.term @z_m += 1.0 states [@x_p, @x_p]
ode.system @Square {
  ode.state @x
  ode.state @z
  ode.term @x += -1.0
  ode.term @z += -1.0 states [@x, @x]
}
