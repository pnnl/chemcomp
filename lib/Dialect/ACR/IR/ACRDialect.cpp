//===- ACRDialect.cpp - ACR dialect ---------------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#include "Standalone/Dialect/ACR/ACRDialect.h"
#include "Standalone/Dialect/ACR/ACROps.h"
#include "Standalone/Dialect/ACR/ACRTypes.h"

using namespace mlir;
using namespace mlir::acr;

#include "Standalone/Dialect/ACR/ACROpsDialect.cpp.inc"

//===----------------------------------------------------------------------===//
// ACR dialect.
//===----------------------------------------------------------------------===//

void ACRDialect::initialize()
{
  addOperations<
#define GET_OP_LIST
#include "Standalone/Dialect/ACR/ACROps.cpp.inc"
      >();
  registerTypes();
}
