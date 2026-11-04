//===- ODERulesGen.h - ODE rule table generator -----------------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#ifndef CHEMCOMP_TBLGEN_ODERULESGEN_H
#define CHEMCOMP_TBLGEN_ODERULESGEN_H

namespace llvm {
class raw_ostream;
class RecordKeeper;
} // namespace llvm

/// Emit the composable ODE rule table. Returns true on error.
bool emitODERules(llvm::raw_ostream &os, const llvm::RecordKeeper &records);

#endif // CHEMCOMP_TBLGEN_ODERULESGEN_H
