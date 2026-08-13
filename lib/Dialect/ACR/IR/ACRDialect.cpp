//===- ACRDialect.cpp - ACR dialect ---------------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#include "ChemComp/Dialect/ACR/ACRDialect.h"
#include "ChemComp/Dialect/ACR/ACROps.h"
#include "ChemComp/Dialect/ACR/ACRTypes.h"

using namespace mlir;
using namespace mlir::cc::acr;

#include "ChemComp/Dialect/ACR/ACROpsDialect.cpp.inc"

//===----------------------------------------------------------------------===//
// ACR dialect.
//===----------------------------------------------------------------------===//

void ACRDialect::initialize() {
  addOperations<
#define GET_OP_LIST
#include "ChemComp/Dialect/ACR/ACROps.cpp.inc"
      >();
  registerTypes();
}
