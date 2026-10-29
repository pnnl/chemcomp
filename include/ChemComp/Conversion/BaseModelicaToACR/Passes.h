//===- Passes.h - Base Modelica to ACR conversion passes --------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#ifndef CHEMCOMP_CONVERSION_BASEMODELICATOACR_PASSES_H
#define CHEMCOMP_CONVERSION_BASEMODELICATOACR_PASSES_H

#include "ChemComp/Dialect/ACR/ACRDialect.h"
#include "mlir/IR/BuiltinOps.h"
#include "mlir/Pass/Pass.h"

namespace mlir {
namespace cc {

#define GEN_PASS_DECL_CONVERTBASEMODELICATOACR
#define GEN_PASS_REGISTRATION
#include "ChemComp/Conversion/BaseModelicaToACR/Passes.h.inc"

} // namespace cc
} // namespace mlir

#endif // CHEMCOMP_CONVERSION_BASEMODELICATOACR_PASSES_H
