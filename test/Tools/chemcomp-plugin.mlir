// RUN: mlir-opt %s --load-dialect-plugin=%standalone_libs/ChemCompPlugin%shlibext --pass-pipeline="builtin.module(acr-switch-bar-foo)" | FileCheck %s

//===----------------------------------------------------------------------===//
// Dialect plugin load path.
//
// Upstream mlir-opt knows nothing about the acr dialect, so parsing the
// acr.stoich_coeff op below only succeeds if the plugin registered the dialect.
// The pipeline additionally proves the plugin's pass registration took effect.
//===----------------------------------------------------------------------===//

module {
  // CHECK-LABEL: func @foo()
  func.func @bar() {
    return
  }

  // The result is returned rather than dropped: acr.stoich_coeff is Pure, so
  // an unused result would be folded away before FileCheck sees it.
  // CHECK-LABEL: func @acr_types
  // CHECK: acr.stoich_coeff -2 : f32
  func.func @acr_types() -> f32 {
    %0 = acr.stoich_coeff -2 : f32
    return %0 : f32
  }
}
