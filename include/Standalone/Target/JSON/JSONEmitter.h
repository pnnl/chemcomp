//===- JSONEmitter.h - ACR CRN JSON translation --------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#ifndef STANDALONE_TARGET_JSON_JSONEMITTER_H
#define STANDALONE_TARGET_JSON_JSONEMITTER_H

#include "mlir/IR/BuiltinOps.h"
#include "mlir/IR/OwningOpRef.h"
#include "mlir/Support/LogicalResult.h"
#include "llvm/Support/SourceMgr.h"

namespace mlir {
class MLIRContext;
class Operation;

namespace acr {

/// Translate ACR CRN MLIR operations to CRN JSON format.
LogicalResult translateToJSON(Operation *op, llvm::raw_ostream &os);

/// Translate CRN JSON format to ACR CRN MLIR operations.
OwningOpRef<Operation *> translateFromJSON(llvm::SourceMgr &sourceMgr,
                                           MLIRContext *context);

/// Registration functions for mlir-translate.
void registerToJSONTranslation();
void registerFromJSONTranslation();

} // namespace acr
} // namespace mlir

#endif // STANDALONE_TARGET_JSON_JSONEMITTER_H
