//===- MotifDBGen.h - Motif database generator ------------------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#ifndef CHEMCOMP_TBLGEN_MOTIFDBGEN_H
#define CHEMCOMP_TBLGEN_MOTIFDBGEN_H

namespace llvm {
class raw_ostream;
class RecordKeeper;
} // namespace llvm

/// Emit the static motif fixture table. Returns true on error.
bool emitMotifDB(llvm::raw_ostream &os, llvm::RecordKeeper &records);

#endif // CHEMCOMP_TBLGEN_MOTIFDBGEN_H
