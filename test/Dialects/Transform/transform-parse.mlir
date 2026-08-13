// Checks that transform ops from LoopExtensionOps and SCFTransformOps can be
// used together.

// CHECK-LABEL: func @test_mixed_loop_extension_scf_transform(
func.func @test_mixed_loop_extension_scf_transform(%arg: tensor<?xf32>) -> tensor<?xf32> {
  %lb = "test.foo"() : () -> (index)
  %ub = "test.foo"() : () -> (index)
  %step = "test.foo"() : () -> (index)
  // CHECK: scf.for
  // CHECK: scf.for
  %0 = scf.for %iv = %lb to %ub step %step iter_args(%t = %arg) -> (tensor<?xf32>) {
    %u = arith.constant 1.0 : f32
    %1 = "test.foo"(%t) : (tensor<?xf32>) -> (tensor<?xf32>)
    scf.yield %1 : tensor<?xf32>
  }
  return %0 : tensor<?xf32>
}

module attributes {transform.with_named_sequence} {
  transform.named_sequence @__transform_main(%arg0: !transform.any_op {transform.readonly}) {
    %0 = transform.structured.match ops{["scf.for"]} in %arg0 : (!transform.any_op) -> !transform.any_op
    // transform.debug.emit_remark_at %0, "Found a for loop" : !transform.any_op
    transform.loop.hoist_loop_invariant_subsets %0 : !transform.any_op
    transform.loop.unroll %0 { factor = 4 } : !transform.any_op
    transform.yield
  }
}