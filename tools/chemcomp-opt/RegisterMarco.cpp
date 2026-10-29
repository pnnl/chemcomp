//===- RegisterMarco.cpp - MARCO registration for chemcomp-opt ------------===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#include "RegisterMarco.h"

#include "ChemComp/Conversion/BaseModelicaToACR/Passes.h"

#include "marco/Dialect/BaseModelica/IR/Dialect.h"
#include "marco/Dialect/BaseModelica/Transforms/AllInterfaces.h"
#include "marco/Dialect/BaseModelica/Transforms/Passes.h"
#include "marco/Dialect/Modeling/IR/Dialect.h"

void mlir::cc::registerMarco(DialectRegistry &registry) {
  registry.insert<mlir::bmodelica::BaseModelicaDialect,
                  mlir::modeling::ModelingDialect>();
  mlir::bmodelica::registerAllDialectInterfaceImplementations(registry);
  mlir::bmodelica::registerBaseModelicaPasses();
  mlir::cc::registerBaseModelicaToACRPasses();
}
