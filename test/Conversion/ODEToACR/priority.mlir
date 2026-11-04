// RUN: chemcomp-opt %s --convert-ode-to-acr | FileCheck %s

// der(x) = -k*x, der(y) = k*x can be covered by conversion X -> Y (one
// reaction) or by degradation X -> 0 plus catalytic production X -> X + Y.
// Rules are tried by priority, so the larger motif wins.

// CHECK-LABEL: acr.network @Conversion {
// CHECK:         acr.rxn @x -> @y rate 0.5 {from = "der(x): -k*x; der(y): k*x", rate_name = "k", rule = "conversion"}
// CHECK-NOT:     acr.rxn
// CHECK:       }
ode.system @Conversion {
  ode.state @x initial 1.0
  ode.state @y initial 0.0
  ode.param @k = 0.5
  ode.term @x += -1.0 params [@k] states [@x]
  ode.term @y += 1.0 params [@k] states [@x]
}

// With different rate constants the two terms are not one conversion: each
// is claimed by its own one-term rule.
// CHECK-LABEL: acr.network @NotAConversion {
// CHECK:         acr.rxn @x -> @x + @y rate 0.25 {from = "der(y): k2*x", rate_name = "k2", rule = "catalytic-production"}
// CHECK:         acr.rxn @x -> none rate 0.5 {from = "der(x): -k1*x", rate_name = "k1", rule = "degradation"}
ode.system @NotAConversion {
  ode.state @x initial 1.0
  ode.state @y initial 0.0
  ode.param @k1 = 0.5
  ode.param @k2 = 0.25
  ode.term @x += -1.0 params [@k1] states [@x]
  ode.term @y += 1.0 params [@k2] states [@x]
}

// Numeric coefficients become the rate; parameters without a value leave the
// rate out but keep its expression.
// CHECK-LABEL: acr.network @Symbolic {
// CHECK:         acr.rxn @x -> none {from = "der(x): -0.5*kd*x", rate_name = "0.5*kd", rule = "degradation"}
// CHECK:         acr.rxn none -> @x rate 2.0 {from = "der(x): 2", rate_name = "2", rule = "production"}
ode.system @Symbolic {
  ode.state @x
  ode.param @kd
  ode.term @x += 2.0
  ode.term @x += -0.5 params [@kd] states [@x]
}
