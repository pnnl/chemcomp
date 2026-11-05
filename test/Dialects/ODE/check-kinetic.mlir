// RUN: chemcomp-opt %s --ode-check-kinetic -split-input-file -verify-diagnostics
// RUN: not chemcomp-opt %s --ode-check-kinetic=require -split-input-file 2>&1 \
// RUN:   | FileCheck %s --check-prefix=REQUIRE

// Hungarian Lemma: every negative term of der(x) must contain x.

// expected-remark @below {{system 'LotkaVolterra' is kinetic: every negative term of der(x) contains x}}
ode.system @LotkaVolterra {
  ode.state @x
  ode.state @y
  ode.term @x += 1.0 states [@x]
  ode.term @x += -1.0 states [@x, @y]
  ode.term @y += 1.0 states [@x, @y]
  ode.term @y += -1.0 states [@y]
}

// -----

// The harmonic oscillator der(x) = -y, der(y) = x: -y in der(x) would have to
// consume x at a rate that does not depend on x.
ode.system @Oscillator {
  ode.state @x
  ode.state @y
  // expected-remark @below {{negative cross-effect: der(x) has -y, which does not contain x}}
  ode.term @x += -1.0 states [@y]
  ode.term @y += 1.0 states [@x]
}

// REQUIRE: error: negative cross-effect: der(x) has -y, which does not contain x
