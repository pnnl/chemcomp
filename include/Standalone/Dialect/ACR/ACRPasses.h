//===- ACRPasses.h - ACR passes  ------------------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
#ifndef ACR_ACRPASSES_H
#define ACR_ACRPASSES_H

#include "Standalone/Dialect/ACR/ACRDialect.h"
#include "Standalone/Dialect/ACR/ACROps.h"
#include "mlir/Pass/Pass.h"
#include <memory>

namespace mlir {
namespace acr {
#define GEN_PASS_DECL
#include "Standalone/Dialect/ACR/ACRPasses.h.inc"

#define GEN_PASS_REGISTRATION
#include "Standalone/Dialect/ACR/ACRPasses.h.inc"
} // namespace acr
} // namespace mlir

#endif
