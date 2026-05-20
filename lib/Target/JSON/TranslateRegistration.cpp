//===- TranslateRegistration.cpp - Register JSON translations --*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#include "Standalone/Dialect/ACR/ACRDialect.h"
#include "Standalone/Target/JSON/JSONEmitter.h"
#include "mlir/Dialect/Arith/IR/Arith.h"
#include "mlir/IR/BuiltinOps.h"
#include "mlir/IR/Dialect.h"
#include "mlir/Tools/mlir-translate/Translation.h"

using namespace mlir;

void mlir::acr::registerToJSONTranslation() {
  TranslateFromMLIRRegistration reg(
      "mlir-to-json", "export ACR CRN to JSON format",
      [](Operation *op, raw_ostream &output) {
        return translateToJSON(op, output);
      },
      [](DialectRegistry &registry) {
        registry.insert<acr::ACRDialect, arith::ArithDialect>();
      });
}

void mlir::acr::registerFromJSONTranslation() {
  TranslateToMLIRRegistration reg(
      "json-to-mlir", "import JSON to ACR CRN MLIR",
      [](llvm::SourceMgr &sourceMgr,
         MLIRContext *context) -> OwningOpRef<Operation *> {
        return translateFromJSON(sourceMgr, context);
      },
      [](DialectRegistry &registry) {
        registry.insert<acr::ACRDialect, arith::ArithDialect>();
      });
}
