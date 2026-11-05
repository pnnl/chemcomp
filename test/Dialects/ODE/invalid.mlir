// RUN: chemcomp-opt %s -split-input-file -verify-diagnostics

ode.system @S {
  ode.state @x
  // expected-error @below {{derivative of 'k', which is not an ode.state of this system}}
  ode.term @k += 1.0
  ode.param @k = 1.0
}

// -----

ode.system @S {
  ode.state @x
  // expected-error @below {{'@y' is not an ode.state of this system}}
  ode.term @x += 1.0 states [@y]
}

// -----

ode.system @S {
  ode.state @x
  ode.state @k
  // expected-error @below {{'@k' is not an ode.param of this system}}
  ode.term @x += 1.0 params [@k]
}

// -----

ode.system @S {
  // expected-error @below {{parameter values are rate constants and must not be negative}}
  ode.param @k = -1.0
}

// -----

ode.system @S {
  // expected-error @below {{'arith.constant' op is not allowed in an ode.system}}
  %0 = arith.constant 1 : i32
}
