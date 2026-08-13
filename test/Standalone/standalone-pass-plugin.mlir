// RUN: mlir-opt %s --load-pass-plugin=%standalone_libs/StandalonePlugin%shlibext --pass-pipeline="builtin.module(acr-switch-bar-foo)" | FileCheck %s

//===----------------------------------------------------------------------===//
// Pass plugin load path.
//
// acr-switch-bar-foo is not linked into upstream mlir-opt, so the pipeline
// resolves only if the plugin registered it. Unlike the dialect-plugin test
// this input parses without any custom dialect, isolating pass registration.
//===----------------------------------------------------------------------===//

module {
  // CHECK-LABEL: func @foo()
  func.func @bar() {
    return
  }

  // CHECK-LABEL: func @abar()
  func.func @abar() {
    return
  }
}
