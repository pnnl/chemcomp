//===- ACRTypes.cpp - ACR dialect types -----------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#include "Standalone/Dialect/ACR/ACRTypes.h"

#include "Standalone/Dialect/ACR/ACRDialect.h"
#include "mlir/IR/Builders.h"
#include "mlir/IR/DialectImplementation.h"
#include "llvm/ADT/TypeSwitch.h"

using namespace mlir::acr;

#define GET_TYPEDEF_CLASSES
#include "Standalone/Dialect/ACR/ACROpsTypes.cpp.inc"

void ACRDialect::registerTypes()
{
  addTypes<
#define GET_TYPEDEF_LIST
#include "Standalone/Dialect/ACR/ACROpsTypes.cpp.inc"
      >();
}
