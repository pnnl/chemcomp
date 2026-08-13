//===- CRNToDB.h - Abstract CRN to database matching ------------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#ifndef CHEMCOMP_CONVERSION_CRNTODB_CRNTODB_H
#define CHEMCOMP_CONVERSION_CRNTODB_CRNTODB_H

#include "ChemComp/Dialect/ACR/ACRDialect.h"
#include "ChemComp/Dialect/ACR/ACROps.h"
#include "mlir/IR/BuiltinOps.h"
#include "mlir/Pass/Pass.h"
#include <memory>

namespace mlir {
namespace cc {

#define GEN_PASS_DECL_CRNDBMATCH
#include "ChemComp/Conversion/Passes.h.inc"

} // namespace cc
} // namespace mlir

#endif // CHEMCOMP_CONVERSION_CRNTODB_CRNTODB_H
