// RUN: chemcomp-translate -mlir-to-json %s | FileCheck %s

// Catalytic reactions must survive the export. X -> 2X has net change +1 for
// X, which split by sign would read back as 0 -> X (a constant source instead
// of autocatalysis). With reactant_matrix / product_matrix recorded, the export
// uses them. initial_values become the populations, rate_values keep their
// double precision and ode_state_names / ode_state_signs give species_map.

acr.crn @lv(%sv: tensor<2x!acr.species>, %rv: tensor<3xf32>)
    -> tensor<2x!acr.species> {
  %out = acr.crn_generic
      {stoich_matrix = dense<[[1, 0], [-1, 1], [0, -1]]> : tensor<3x2xi32>,
       reactant_matrix = dense<[[1, 0], [1, 1], [0, 1]]> : tensor<3x2xi32>,
       product_matrix = dense<[[2, 0], [0, 2], [0, 0]]> : tensor<3x2xi32>,
       rate_values = array<f64: 1.1, 0.4, 0.123456789>,
       initial_values = array<f64: 2.0, 0.5>,
       ode_state_names = ["x", "x"],
       ode_state_signs = array<i64: 1, -1>}
      (%sv, %rv) : (tensor<2x!acr.species>, tensor<3xf32>)
      -> tensor<2x!acr.species>
  acr.return %out : tensor<2x!acr.species>
}

// R1: X -> 2X
// CHECK:      "reactants": [
// CHECK-NEXT:   1,
// CHECK-NEXT:   0
// CHECK:      "products": [
// CHECK-NEXT:   2,
// CHECK-NEXT:   0
// CHECK:      "rate": 1.1,
// R2: X + Y -> 2Y
// CHECK:      "reactants": [
// CHECK-NEXT:   1,
// CHECK-NEXT:   1
// CHECK:      "products": [
// CHECK-NEXT:   0,
// CHECK-NEXT:   2
// CHECK:      "rate": 0.4,
// CHECK:      "rate": 0.123456789,
// CHECK:      "populations": {
// CHECK-NEXT:   "S0": 2,
// CHECK-NEXT:   "S1": 0.5
// CHECK:      "species_map": {
// CHECK-NEXT:   "S0": {
// CHECK-NEXT:     "state": "x",
// CHECK-NEXT:     "sign": 1
// CHECK:        "S1": {
// CHECK-NEXT:     "state": "x",
// CHECK-NEXT:     "sign": -1
