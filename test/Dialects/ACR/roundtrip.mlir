// RUN: standalone-opt %s -split-input-file | FileCheck %s

//===----------------------------------------------------------------------===//
// acr.stoich_coeff
//===----------------------------------------------------------------------===//

// CHECK-LABEL: func.func @test_stoich_coeff
func.func @test_stoich_coeff() {
  // CHECK: acr.stoich_coeff -2 : f32
  %h2_stoich  = acr.stoich_coeff -2 : f32
  // CHECK: acr.stoich_coeff -1 : f32
  %o2_stoich  = acr.stoich_coeff -1 : f32
  // CHECK: acr.stoich_coeff 2 : f32
  %h2o_stoich = acr.stoich_coeff  2 : f32
  return
}

// -----

//===----------------------------------------------------------------------===//
// acr.react — scalar form (individual species + stoich f32 values)
//===----------------------------------------------------------------------===//

// CHECK-LABEL: func.func @test_react_scalar
func.func @test_react_scalar(%h2: !acr.species, %o2: !acr.species,
                              %h2o: !acr.species, %k: f32) {
  %hs = acr.stoich_coeff -2 : f32
  %os = acr.stoich_coeff -1 : f32
  %ws = acr.stoich_coeff  2 : f32
  // CHECK: acr.react
  // CHECK-SAME: !acr.species, !acr.species, !acr.species, f32, f32, f32, f32
  // CHECK-SAME: -> (!acr.species, !acr.species, !acr.species)
  %h2_o, %o2_o, %h2o_o = acr.react(%h2, %o2, %h2o, %hs, %os, %ws, %k)
      : (!acr.species, !acr.species, !acr.species, f32, f32, f32, f32)
     -> (!acr.species, !acr.species, !acr.species)
  return
}

// -----

//===----------------------------------------------------------------------===//
// acr.reaction — definition with body
//===----------------------------------------------------------------------===//

// CHECK-LABEL: acr.reaction @reaction_1
// CHECK-SAME: (%{{.*}}: !acr.species, %{{.*}}: !acr.species, %{{.*}}: !acr.species, %{{.*}}: f32)
acr.reaction @reaction_1(%h2: !acr.species, %o2: !acr.species,
                          %h2o: !acr.species, %k: f32)
    -> (!acr.species, !acr.species, !acr.species) {
  %hs = acr.stoich_coeff -2 : f32
  %os = acr.stoich_coeff -1 : f32
  %ws = acr.stoich_coeff  2 : f32
  // CHECK: acr.react
  %h2_o, %o2_o, %h2o_o = acr.react(%h2, %o2, %h2o, %hs, %os, %ws, %k)
      : (!acr.species, !acr.species, !acr.species, f32, f32, f32, f32)
     -> (!acr.species, !acr.species, !acr.species)
  // CHECK: acr.return
  acr.return %h2_o, %o2_o, %h2o_o : !acr.species, !acr.species, !acr.species
}

// -----

//===----------------------------------------------------------------------===//
// acr.reaction — private (external) declaration
//===----------------------------------------------------------------------===//

// CHECK-LABEL: acr.reaction private @reaction_1
acr.reaction private @reaction_1(!acr.species, !acr.species, !acr.species, f32)
    -> (!acr.species, !acr.species, !acr.species)

// -----

//===----------------------------------------------------------------------===//
// acr.crn — with acr.reaction_instance calls
//===----------------------------------------------------------------------===//

acr.reaction private @r1(!acr.species, !acr.species, !acr.species, f32)
    -> (!acr.species, !acr.species, !acr.species)

acr.reaction private @r2(!acr.species, !acr.species, !acr.species, f32)
    -> (!acr.species, !acr.species, !acr.species)

// CHECK-LABEL: acr.crn @reaction_network
// CHECK-SAME: (%{{.*}}: !acr.species, %{{.*}}: !acr.species
acr.crn @reaction_network(%h2: !acr.species, %o2: !acr.species,
                           %co: !acr.species, %h2o: !acr.species,
                           %co2: !acr.species, %k1: f32, %k2: f32)
    -> (!acr.species, !acr.species, !acr.species, !acr.species, !acr.species) {
  // CHECK: acr.reaction_instance @r1
  %h2_o, %o2_o, %h2o_o =
      acr.reaction_instance @r1(%h2, %o2, %h2o, %k1)
          : (!acr.species, !acr.species, !acr.species, f32)
         -> (!acr.species, !acr.species, !acr.species)

  // CHECK: acr.reaction_instance @r2
  %co_o, %o2_f, %co2_o =
      acr.reaction_instance @r2(%co, %o2_o, %co2, %k2)
          : (!acr.species, !acr.species, !acr.species, f32)
         -> (!acr.species, !acr.species, !acr.species)

  // CHECK: acr.return
  acr.return %h2_o, %o2_f, %h2o_o, %co_o, %co2_o
      : !acr.species, !acr.species, !acr.species, !acr.species, !acr.species
}

// -----

//===----------------------------------------------------------------------===//
// acr.crn_generic — compact stoichiometric matrix form
//===----------------------------------------------------------------------===//

// CHECK-LABEL: acr.crn @generic_network
acr.crn @generic_network(%s0: !acr.species, %s1: !acr.species, %s2: !acr.species,
                          %r0: f32, %r1: f32)
    -> (!acr.species, !acr.species, !acr.species) {
  // CHECK: acr.crn_generic
  // CHECK-SAME: stoich_matrix
  // CHECK-SAME: tensor<2x3xi32>
  %out0, %out1, %out2 = acr.crn_generic
      {stoich_matrix = dense<[[-2, -1, 2], [0, -1, 1]]> : tensor<2x3xi32>}
      (%s0, %s1, %s2, %r0, %r1)
          : (!acr.species, !acr.species, !acr.species, f32, f32)
         -> (!acr.species, !acr.species, !acr.species)
  acr.return %out0, %out1, %out2 : !acr.species, !acr.species, !acr.species
}

// -----

//===----------------------------------------------------------------------===//
// acr.react — tensor form
//===----------------------------------------------------------------------===//

// CHECK-LABEL: func.func @test_react_tensor
func.func @test_react_tensor(%sv: tensor<3x!acr.species>, %stoich: tensor<3xf32>,
                              %k: f32) {
  // CHECK: acr.react
  // CHECK-SAME: tensor<3x!acr.species>, tensor<3xf32>, f32
  // CHECK-SAME: -> tensor<3x!acr.species>
  %sv_out = acr.react(%sv, %stoich, %k)
      : (tensor<3x!acr.species>, tensor<3xf32>, f32) -> tensor<3x!acr.species>
  return
}

// -----

//===----------------------------------------------------------------------===//
// acr.crn_generic — tensor form
//===----------------------------------------------------------------------===//

// CHECK-LABEL: acr.crn @generic_network_tensor
acr.crn @generic_network_tensor(%sv: tensor<3x!acr.species>, %rv: tensor<2xf32>)
    -> tensor<3x!acr.species> {
  // CHECK: acr.crn_generic
  // CHECK-SAME: stoich_matrix
  // CHECK-SAME: tensor<2x3xi32>
  // CHECK-SAME: tensor<3x!acr.species>, tensor<2xf32>
  // CHECK-SAME: -> tensor<3x!acr.species>
  %out = acr.crn_generic
      {stoich_matrix = dense<[[-2, -1, 2], [0, -1, 1]]> : tensor<2x3xi32>}
      (%sv, %rv) : (tensor<3x!acr.species>, tensor<2xf32>) -> tensor<3x!acr.species>
  acr.return %out : tensor<3x!acr.species>
}

// -----

//===----------------------------------------------------------------------===//
// acr.reaction — tensor form definition
//===----------------------------------------------------------------------===//

// CHECK-LABEL: acr.reaction @reaction_tensor
// CHECK-SAME: (%{{.*}}: tensor<3x!acr.species>, %{{.*}}: f32)
acr.reaction @reaction_tensor(%sv: tensor<3x!acr.species>, %k: f32)
    -> tensor<3x!acr.species> {
  %stoich = arith.constant dense<[-2.0, -1.0, 2.0]> : tensor<3xf32>
  // CHECK: acr.react
  %sv_out = acr.react(%sv, %stoich, %k)
      : (tensor<3x!acr.species>, tensor<3xf32>, f32) -> tensor<3x!acr.species>
  // CHECK: acr.return
  acr.return %sv_out : tensor<3x!acr.species>
}
