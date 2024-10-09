//===- StandaloneTypes.cpp - Standalone dialect types -----------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#include "Standalone/Dialect/Standalone/StandaloneTypes.h"

#include "Standalone/Dialect/Standalone/StandaloneDialect.h"
#include "mlir/IR/Builders.h"
#include "mlir/IR/DialectImplementation.h"
#include "llvm/ADT/TypeSwitch.h"

using namespace mlir::standalone;

#define GET_TYPEDEF_CLASSES
#include "Standalone/Dialect/Standalone/StandaloneOpsTypes.cpp.inc"

void StandaloneDialect::registerTypes() {
  addTypes<
#define GET_TYPEDEF_LIST
#include "Standalone/Dialect/Standalone/StandaloneOpsTypes.cpp.inc"
      >();
}
