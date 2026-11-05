// RUN: chemcomp-opt %s --convert-ode-to-acr -verify-diagnostics

// Without dual-rail encoding a system with negative cross-effects has no
// mass-action network (Hungarian Lemma), and the conversion says why.

// expected-error @below {{system 'Oscillator' has no mass-action network: der(x) has -y, which does not contain x}}
// expected-note @below {{run ode-dual-rail first (the ode-to-crn pipeline does), or rewrite the model}}
ode.system @Oscillator {
  ode.state @x
  ode.state @y
  ode.term @x += -1.0 states [@y]
  ode.term @y += 1.0 states [@x]
}
