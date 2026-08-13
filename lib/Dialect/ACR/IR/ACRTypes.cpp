//===- ACRTypes.cpp - ACR dialect types -----------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#include "ChemComp/Dialect/ACR/ACRTypes.h"

#include "ChemComp/Dialect/ACR/ACRDialect.h"
#include "mlir/IR/Builders.h"
#include "mlir/IR/DialectImplementation.h"
#include "llvm/ADT/TypeSwitch.h"

using namespace mlir::cc::acr;

#define GET_TYPEDEF_CLASSES
#include "ChemComp/Dialect/ACR/ACROpsTypes.cpp.inc"

void ACRDialect::registerTypes() {
  addTypes<
#define GET_TYPEDEF_LIST
#include "ChemComp/Dialect/ACR/ACROpsTypes.cpp.inc"
      >();
}
