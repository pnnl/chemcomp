//===- ODEDialect.cpp - ODE dialect ---------------------------------------===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#include "ChemComp/Dialect/ODE/ODEDialect.h"
#include "ChemComp/Dialect/ODE/ODEOps.h"

using namespace mlir;
using namespace mlir::cc::ode;

#include "ChemComp/Dialect/ODE/ODEOpsDialect.cpp.inc"

void ODEDialect::initialize() {
  addOperations<
#define GET_OP_LIST
#include "ChemComp/Dialect/ODE/ODEOps.cpp.inc"
      >();
}
