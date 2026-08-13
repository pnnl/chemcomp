// RUN: chemcomp-opt -split-input-file -verify-diagnostics %s

//===----------------------------------------------------------------------===//
// acr.react — tensor form: mismatched stoich dimension
//===----------------------------------------------------------------------===//

func.func @react_tensor_stoich_mismatch(%sv: tensor<3x!acr.species>,
                                         %stoich: tensor<4xf32>, %k: f32) {
  // expected-error @+1 {{stoichiometry tensor dimension (4) must match species tensor dimension (3)}}
  %out = acr.react(%sv, %stoich, %k)
      : (tensor<3x!acr.species>, tensor<4xf32>, f32) -> tensor<3x!acr.species>
  return
}

// -----

//===----------------------------------------------------------------------===//
// acr.react — tensor form: wrong number of inputs
//===----------------------------------------------------------------------===//

func.func @react_tensor_extra_input(%sv: tensor<3x!acr.species>,
                                     %stoich: tensor<3xf32>, %k: f32,
                                     %extra: f32) {
  // expected-error @+1 {{tensor form expects exactly 3 inputs}}
  %out = acr.react(%sv, %stoich, %k, %extra)
      : (tensor<3x!acr.species>, tensor<3xf32>, f32, f32)
     -> tensor<3x!acr.species>
  return
}

// -----

//===----------------------------------------------------------------------===//
// acr.react — scalar form: wrong number of results
//===----------------------------------------------------------------------===//

func.func @react_scalar_result_mismatch(%h2: !acr.species, %o2: !acr.species,
                                         %hs: f32, %os: f32, %k: f32) {
  // expected-error @+1 {{scalar form expects 2 results, got 1}}
  %out = acr.react(%h2, %o2, %hs, %os, %k)
      : (!acr.species, !acr.species, f32, f32, f32) -> !acr.species
  return
}

// -----

//===----------------------------------------------------------------------===//
// acr.crn_generic — missing stoich_matrix attribute
//===----------------------------------------------------------------------===//

func.func @crn_generic_missing_stoich(%sv: tensor<3x!acr.species>,
                                       %rv: tensor<2xf32>) {
  // expected-error @+1 {{requires a 'stoich_matrix' dense integer attribute}}
  %out = acr.crn_generic
      (%sv, %rv) : (tensor<3x!acr.species>, tensor<2xf32>) -> tensor<3x!acr.species>
  return
}

// -----

//===----------------------------------------------------------------------===//
// acr.crn_generic — stoich_matrix must be 2-D
//===----------------------------------------------------------------------===//

func.func @crn_generic_stoich_rank1(%sv: tensor<3x!acr.species>,
                                     %rv: tensor<2xf32>) {
  // expected-error @+1 {{'stoich_matrix' must be a 2-D tensor}}
  %out = acr.crn_generic
      {stoich_matrix = dense<[-2, -1, 2]> : tensor<3xi32>}
      (%sv, %rv) : (tensor<3x!acr.species>, tensor<2xf32>) -> tensor<3x!acr.species>
  return
}

// -----

//===----------------------------------------------------------------------===//
// acr.crn_generic — tensor form: species dimension mismatch with stoich_matrix
//===----------------------------------------------------------------------===//

func.func @crn_generic_species_dim_mismatch(%sv: tensor<4x!acr.species>,
                                             %rv: tensor<2xf32>) {
  // expected-error @+1 {{species tensor dimension (4) must match stoich_matrix columns (3)}}
  %out = acr.crn_generic
      {stoich_matrix = dense<[[-2, -1, 2], [0, -1, 1]]> : tensor<2x3xi32>}
      (%sv, %rv) : (tensor<4x!acr.species>, tensor<2xf32>) -> tensor<4x!acr.species>
  return
}

// -----

//===----------------------------------------------------------------------===//
// acr.crn_generic — tensor form: rate dimension mismatch with stoich_matrix
//===----------------------------------------------------------------------===//

func.func @crn_generic_rate_dim_mismatch(%sv: tensor<3x!acr.species>,
                                          %rv: tensor<5xf32>) {
  // expected-error @+1 {{rate tensor dimension (5) must match stoich_matrix rows (2)}}
  %out = acr.crn_generic
      {stoich_matrix = dense<[[-2, -1, 2], [0, -1, 1]]> : tensor<2x3xi32>}
      (%sv, %rv) : (tensor<3x!acr.species>, tensor<5xf32>) -> tensor<3x!acr.species>
  return
}

// -----

//===----------------------------------------------------------------------===//
// acr.crn_generic — scalar form: wrong number of inputs
//===----------------------------------------------------------------------===//

func.func @crn_generic_scalar_input_mismatch(%s0: !acr.species, %s1: !acr.species,
                                              %r0: f32) {
  // expected-error @+1 {{scalar form expects 5 inputs (3 species + 2 rates), got 3}}
  %out0, %out1, %out2 = acr.crn_generic
      {stoich_matrix = dense<[[-2, -1, 2], [0, -1, 1]]> : tensor<2x3xi32>}
      (%s0, %s1, %r0)
          : (!acr.species, !acr.species, f32)
         -> (!acr.species, !acr.species, !acr.species)
  return
}
