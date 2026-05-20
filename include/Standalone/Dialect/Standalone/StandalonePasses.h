//===- StandalonePasses.h - Standalone passes  ------------------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
#ifndef STANDALONE_STANDALONEPASSES_H
#define STANDALONE_STANDALONEPASSES_H

#include "Standalone/Dialect/Standalone/StandaloneDialect.h"
#include "Standalone/Dialect/Standalone/StandaloneOps.h"
#include "mlir/Pass/Pass.h"
#include <memory>

namespace mlir {
namespace standalone {
#define GEN_PASS_DECL
#include "Standalone/Dialect/Standalone/StandalonePasses.h.inc"

#define GEN_PASS_REGISTRATION
#include "Standalone/Dialect/Standalone/StandalonePasses.h.inc"
} // namespace standalone
} // namespace mlir

#endif
