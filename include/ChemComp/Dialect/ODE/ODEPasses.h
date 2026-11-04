//===- ODEPasses.h - ODE dialect passes -------------------------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#ifndef CHEMCOMP_DIALECT_ODE_ODEPASSES_H
#define CHEMCOMP_DIALECT_ODE_ODEPASSES_H

#include "ChemComp/Dialect/ODE/ODEDialect.h"
#include "ChemComp/Dialect/ODE/ODEOps.h"
#include "mlir/IR/BuiltinOps.h"
#include "mlir/Pass/Pass.h"

#include <memory>

namespace mlir {
namespace cc {
namespace ode {

#define GEN_PASS_DECL
#include "ChemComp/Dialect/ODE/ODEPasses.h.inc"

#define GEN_PASS_REGISTRATION
#include "ChemComp/Dialect/ODE/ODEPasses.h.inc"

} // namespace ode
} // namespace cc
} // namespace mlir

#endif // CHEMCOMP_DIALECT_ODE_ODEPASSES_H
