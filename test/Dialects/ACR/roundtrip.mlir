// RUN: standalone-opt %s -split-input-file | FileCheck %s

// CHECK-LABEL: func @test_acr_operations
func.func @test_acr_operations() {

  %rate = arith.constant 1.0 : f64
  // CHECK: acr.chemreaction
  %r1 = acr.chemreaction [1,0,0], [0,1,0], %rate : 
    (vector<3xi32>, vector<3xi32>, f64) -> cr<3>
  %r2 = acr.chemreaction [0,1,0], [0,0,1], %rate : 
    (vector<3xi32>, vector<3xi32>, f64) -> cr<3>
  %r3 = acr.chemreaction [0,0,1], [1,0,0], %rate : 
    (vector<3xi32>, vector<3xi32>, f64) -> cr<3>

  return
}