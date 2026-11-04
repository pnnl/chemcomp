//===- ODEOps.h - ODE dialect operations ------------------------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#ifndef CHEMCOMP_DIALECT_ODE_ODEOPS_H
#define CHEMCOMP_DIALECT_ODE_ODEOPS_H

#include "mlir/Bytecode/BytecodeOpInterface.h"
#include "mlir/IR/BuiltinAttributes.h"
#include "mlir/IR/Dialect.h"
#include "mlir/IR/OpDefinition.h"
#include "mlir/IR/RegionKindInterface.h"
#include "mlir/IR/SymbolTable.h"

#include "ChemComp/Dialect/ODE/ODEDialect.h"

#define GET_OP_CLASSES
#include "ChemComp/Dialect/ODE/ODEOps.h.inc"

#endif // CHEMCOMP_DIALECT_ODE_ODEOPS_H
