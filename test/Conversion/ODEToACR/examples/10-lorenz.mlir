// Lorenz system (chaotic):
//   der(x) = sigma*(y - x)
//   der(y) = x*(rho - z) - y
//   der(z) = x*y - beta*z
// 3 states, 7 terms, degree 2, NOT kinetic: -x*z in der(y) does not contain y.

// RUN: chemcomp-opt %s --ode-normalize --ode-dual-rail --convert-ode-to-acr \
// RUN:   --acr-network-merge | FileCheck %s
// RUN: not chemcomp-opt %s --ode-to-crn=canonical-fallback=false 2>&1 \
// RUN:   | FileCheck %s --check-prefix=ONLY
// RUN: chemcomp-opt %s --ode-to-crn | FileCheck %s --check-prefix=CRN
// RUN: chemcomp-opt %s --ode-to-crn | chemcomp-translate -mlir-to-json > %t.crn.json
// RUN: chemcomp-translate -ode-to-json %s > %t.ode.json
// RUN: %if emulator %{ %chemcomp_python %check_realization --ode %t.ode.json \
// RUN:   --crn %t.crn.json --t-end 1 | FileCheck %s --check-prefix=SIM %}

ode.system @Lorenz {
  ode.state @x initial 1.0
  ode.state @y initial 1.0
  ode.state @z initial 1.0
  ode.param @sigma = 10.0
  ode.param @rho = 28.0
  ode.param @beta = 2.6666666666666665
  ode.term @x += 1.0 params [@sigma] states [@y]
  ode.term @x += -1.0 params [@sigma] states [@x]
  ode.term @y += 1.0 params [@rho] states [@x]
  ode.term @y += -1.0 states [@x, @z]
  ode.term @y += -1.0 states [@y]
  ode.term @z += 1.0 states [@x, @y]
  ode.term @z += -1.0 params [@beta] states [@z]
}

// 21 reactions x 6 species is past MLIR's 100-element hex limit; chemcomp-opt
// still prints the matrices as numbers.
// CRN:     acr.crn_generic
// CRN-NOT: 0x
// CRN-SAME: product_matrix = dense<{{\[\[}}
// CRN-SAME: reactant_matrix = dense<{{\[\[}}
// CRN-SAME: stoich_matrix = dense<{{\[\[}}

// CHECK-LABEL: acr.network @Lorenz {
// CHECK-NEXT:    acr.species @x_p initial 1.0 {ode_state = "x", rail_sign = 1 : i64}
// CHECK-NEXT:    acr.species @x_m initial 0.0 {ode_state = "x", rail_sign = -1 : i64}
// CHECK-NEXT:    acr.species @y_p initial 1.0 {ode_state = "y", rail_sign = 1 : i64}
// CHECK-NEXT:    acr.species @y_m initial 0.0 {ode_state = "y", rail_sign = -1 : i64}
// CHECK-NEXT:    acr.species @z_p initial 1.0 {ode_state = "z", rail_sign = 1 : i64}
// CHECK-NEXT:    acr.species @z_m initial 0.0 {ode_state = "z", rail_sign = -1 : i64}
// CHECK-NEXT:    acr.rxn @x_m + @x_p -> none rate 1.0 {from = "der(x_p): -x_m*x_p; der(x_m): -x_m*x_p", rate_name = "1", rule = "annihilation"}
// CHECK-NEXT:    acr.rxn @y_m + @y_p -> none rate 1.0 {from = "der(y_p): -y_m*y_p; der(y_m): -y_m*y_p", rate_name = "1", rule = "annihilation"}
// CHECK-NEXT:    acr.rxn @z_m + @z_p -> none rate 1.0 {from = "der(z_p): -z_m*z_p; der(z_m): -z_m*z_p", rate_name = "1", rule = "annihilation"}
// CHECK-NEXT:    acr.rxn @x_p -> @x_p + @y_p rate 28.0 {from = "der(y_p): rho*x_p", rate_name = "rho", rule = "catalytic-production"}
// CHECK-NEXT:    acr.rxn @x_m -> @x_m + @y_m rate 28.0 {from = "der(y_m): rho*x_m", rate_name = "rho", rule = "catalytic-production"}
// CHECK-NEXT:    acr.rxn @y_p -> @x_p + @y_p rate 10.0 {from = "der(x_p): sigma*y_p", rate_name = "sigma", rule = "catalytic-production"}
// CHECK-NEXT:    acr.rxn @y_m -> @x_m + @y_m rate 10.0 {from = "der(x_m): sigma*y_m", rate_name = "sigma", rule = "catalytic-production"}
// CHECK-NEXT:    acr.rxn @x_p -> none rate 10.0 {from = "der(x_p): -sigma*x_p", rate_name = "sigma", rule = "degradation"}
// CHECK-NEXT:    acr.rxn @x_m -> none rate 10.0 {from = "der(x_m): -sigma*x_m", rate_name = "sigma", rule = "degradation"}
// CHECK-NEXT:    acr.rxn @y_p -> none rate 1.0 {from = "der(y_p): -y_p", rate_name = "1", rule = "degradation"}
// CHECK-NEXT:    acr.rxn @y_m -> none rate 1.0 {from = "der(y_m): -y_m", rate_name = "1", rule = "degradation"}
// CHECK-NEXT:    acr.rxn @z_p -> none rate 2.6666666666666665 {from = "der(z_p): -beta*z_p", rate_name = "beta", rule = "degradation"}
// CHECK-NEXT:    acr.rxn @z_m -> none rate 2.6666666666666665 {from = "der(z_m): -beta*z_m", rate_name = "beta", rule = "degradation"}
// CHECK-NEXT:    acr.rxn @x_m + @z_p -> @x_m + @y_p + @z_p rate 1.0 {from = "der(y_p): x_m*z_p", rate_name = "1", rule = "canonical"}
// CHECK-NEXT:    acr.rxn @x_p + @z_m -> @x_p + @y_p + @z_m rate 1.0 {from = "der(y_p): x_p*z_m", rate_name = "1", rule = "canonical"}
// CHECK-NEXT:    acr.rxn @x_m + @z_m -> @x_m + @y_m + @z_m rate 1.0 {from = "der(y_m): x_m*z_m", rate_name = "1", rule = "canonical"}
// CHECK-NEXT:    acr.rxn @x_p + @z_p -> @x_p + @y_m + @z_p rate 1.0 {from = "der(y_m): x_p*z_p", rate_name = "1", rule = "canonical"}
// CHECK-NEXT:    acr.rxn @x_m + @y_m -> @x_m + @y_m + @z_p rate 1.0 {from = "der(z_p): x_m*y_m", rate_name = "1", rule = "canonical"}
// CHECK-NEXT:    acr.rxn @x_p + @y_p -> @x_p + @y_p + @z_p rate 1.0 {from = "der(z_p): x_p*y_p", rate_name = "1", rule = "canonical"}
// CHECK-NEXT:    acr.rxn @x_m + @y_p -> @x_m + @y_p + @z_m rate 1.0 {from = "der(z_m): x_m*y_p", rate_name = "1", rule = "canonical"}
// CHECK-NEXT:    acr.rxn @x_p + @y_m -> @x_p + @y_m + @z_m rate 1.0 {from = "der(z_m): x_p*y_m", rate_name = "1", rule = "canonical"}
// CHECK-NEXT:  }

// Composition only: no rule covers the bimolecular catalytic terms that
// x*z and x*y become on the rails.
// ONLY: error: no composable rule covers 8 term(s) of system 'Lorenz': der(y_p): x_m*z_p; der(y_p): x_p*z_m; der(y_m): x_m*z_m; der(y_m): x_p*z_p; der(z_p): x_m*y_m; der(z_p): x_p*y_p; der(z_m): x_m*y_p; der(z_m): x_p*y_m

// The emulator simulates the network; its states match the ODE solution.
// SIM: PASS
