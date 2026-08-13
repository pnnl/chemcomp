//===- Passes.h - Conversion Pass Construction and Registration -----------===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#ifndef CHEMCOMP_CONVERSION_PASSES_H
#define CHEMCOMP_CONVERSION_PASSES_H

#include "Standalone/Conversion/CRNToDB/CRNToDB.h"

namespace mlir {

/// Generate the code for registering conversion passes.
#define GEN_PASS_REGISTRATION
#include "Standalone/Conversion/Passes.h.inc"

} // namespace mlir

#endif // CHEMCOMP_CONVERSION_PASSES_H
