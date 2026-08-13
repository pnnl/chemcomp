// RUN: chemcomp-opt %s -acr-switch-bar-foo | FileCheck %s

//===----------------------------------------------------------------------===//
// acr-switch-bar-foo
//
// Guards ACR pass registration: chemcomp-opt must call
// mlir::acr::registerPasses(), otherwise this pass is unreachable from the
// command line and the RUN line fails with "no such pass". Later ACR-adjacent
// passes (crn-db-match) depend on that registration call existing.
//===----------------------------------------------------------------------===//

// CHECK-LABEL: func.func @foo
// CHECK-NOT: func.func @bar
func.func @bar() {
  return
}

// An unrelated symbol must be left alone.
// CHECK-LABEL: func.func @untouched
func.func @untouched() {
  return
}
