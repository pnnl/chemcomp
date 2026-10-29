//===- RegisterMarco.h - MARCO registration for chemcomp-opt ----*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#ifndef CHEMCOMP_TOOLS_CHEMCOMP_OPT_REGISTERMARCO_H
#define CHEMCOMP_TOOLS_CHEMCOMP_OPT_REGISTERMARCO_H

namespace mlir {
class DialectRegistry;

namespace cc {

/// Register the MARCO dialects, external models and passes, plus the
/// conversions that depend on them. Does nothing when ChemComp is configured
/// with CHEMCOMP_ENABLE_MARCO=OFF.
void registerMarco(DialectRegistry &registry);

} // namespace cc
} // namespace mlir

#endif // CHEMCOMP_TOOLS_CHEMCOMP_OPT_REGISTERMARCO_H
