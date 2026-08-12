//===- chemcomp-tblgen.cpp - ChemComp TableGen driver -----------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// TableGen driver for ChemComp-specific backends.
//
// Separate from mlir-tblgen because these backends emit plain data tables
// rather than MLIR IR constructs, and so need nothing from MLIR's TableGen
// library. Dialect and pass definitions continue to use mlir-tblgen.
//
//===----------------------------------------------------------------------===//

#include "llvm/Support/CommandLine.h"
#include "llvm/Support/PrettyStackTrace.h"
#include "llvm/Support/Signals.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/TableGen/Main.h"
#include "llvm/TableGen/Record.h"

#include "MotifDBGen.h"

using namespace llvm;

namespace {

enum ActionType { GenMotifDB };

cl::opt<ActionType> action(
    cl::desc("Action to perform:"),
    cl::values(clEnumValN(GenMotifDB, "gen-motif-db",
                          "Generate the reaction motif fixture database")));

bool chemcompTableGenMain(raw_ostream &os, RecordKeeper &records) {
  switch (action) {
  case GenMotifDB:
    return emitMotifDB(os, records);
  }
  llvm_unreachable("unhandled action");
}

} // namespace

int main(int argc, char **argv) {
  sys::PrintStackTraceOnErrorSignal(argv[0]);
  PrettyStackTraceProgram x(argc, argv);
  cl::ParseCommandLineOptions(argc, argv);

  return TableGenMain(argv[0], &chemcompTableGenMain);
}
