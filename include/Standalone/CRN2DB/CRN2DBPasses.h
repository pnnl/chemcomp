//===- CRN2DBPasses.h - CRN2DB passes ---------------------------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#ifndef STANDALONE_CRN2DB_CRN2DBPASSES_H
#define STANDALONE_CRN2DB_CRN2DBPASSES_H

#include "Standalone/Dialect/ACR/ACRDialect.h"
#include "Standalone/Dialect/ACR/ACROps.h"
#include "mlir/IR/BuiltinOps.h"
#include "mlir/Pass/Pass.h"
#include <memory>

namespace mlir {
namespace crn2db {
#define GEN_PASS_DECL
#include "Standalone/CRN2DB/CRN2DBPasses.h.inc"

#define GEN_PASS_REGISTRATION
#include "Standalone/CRN2DB/CRN2DBPasses.h.inc"
} // namespace crn2db
} // namespace mlir

#endif // STANDALONE_CRN2DB_CRN2DBPASSES_H
