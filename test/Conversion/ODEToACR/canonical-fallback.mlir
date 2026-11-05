// RUN: chemcomp-opt %s --convert-ode-to-acr | FileCheck %s
// RUN: not chemcomp-opt %s --convert-ode-to-acr=canonical-fallback=false 2>&1 \
// RUN:   | FileCheck %s --check-prefix=ONLY

// Lotka-Volterra with different predation rates b1 != b2 is still a
// mass-action system, but no rule claims -b1*x*y and b2*x*y together.
// Rules take the other terms; the canonical rule covers the two left over,
// one reaction each: -b1*x*y in der(x) is X + Y -> Y, b2*x*y in der(y) is
// X + Y -> X + 2Y.

// CHECK-LABEL: acr.network @Unbalanced {
// CHECK:         acr.rxn @x -> 2 @x rate 1.0 {{.*}}rule = "autocatalysis"}
// CHECK:         acr.rxn @y -> none rate 1.0 {{.*}}rule = "degradation"}
// CHECK:         acr.rxn @x + @y -> @y rate 0.3 {from = "der(x): -b1*x*y", rate_name = "b1", rule = "canonical"}
// CHECK:         acr.rxn @x + @y -> @x + 2 @y rate 0.5 {from = "der(y): b2*x*y", rate_name = "b2", rule = "canonical"}

// ONLY: error: no composable rule covers 2 term(s) of system 'Unbalanced': der(x): -b1*x*y; der(y): b2*x*y
// ONLY: note: enable canonical-fallback to cover them with one reaction per term
ode.system @Unbalanced {
  ode.state @x initial 1.0
  ode.state @y initial 1.0
  ode.param @b1 = 0.3
  ode.param @b2 = 0.5
  ode.term @x += 1.0 states [@x]
  ode.term @x += -1.0 params [@b1] states [@x, @y]
  ode.term @y += 1.0 params [@b2] states [@x, @y]
  ode.term @y += -1.0 states [@y]
}
