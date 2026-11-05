//===- chemcomp-opt.cpp -----------------------------------------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#include "mlir/Dialect/Arith/IR/Arith.h"
#include "mlir/Dialect/Func/IR/FuncOps.h"
#include "mlir/IR/MLIRContext.h"
#include "mlir/InitAllDialects.h"
#include "mlir/InitAllPasses.h"
#include "mlir/Support/FileUtilities.h"
#include "mlir/Tools/mlir-opt/MlirOptMain.h"

#include "ChemComp/Dialect/ACR/ACRDialect.h"
#include "ChemComp/Dialect/ACR/ACRPasses.h"
#include "ChemComp/Dialect/ODE/ODEDialect.h"
#include "ChemComp/Dialect/ODE/ODEPasses.h"

#include "ChemComp/Conversion/ODEToACR/Passes.h"
#include "ChemComp/Conversion/Passes.h"

#include "RegisterMarco.h"

#include <vector>

int main(int argc, char **argv) {
  mlir::registerAllPasses();
  mlir::cc::acr::registerPasses();
  mlir::cc::ode::registerPasses();
  mlir::cc::registerConversionPasses();
  mlir::cc::registerODEToACRPasses();
  mlir::cc::registerODEToCRNPipeline();

  mlir::DialectRegistry registry;
  registry.insert<mlir::arith::ArithDialect, mlir::func::FuncDialect,
                  mlir::cc::acr::ACRDialect, mlir::cc::ode::ODEDialect>();
  mlir::cc::registerMarco(registry);
  // Add the following to include *all* MLIR Core dialects, or selectively
  // include what you need like above. You only need to register dialects that
  // will be *parsed* by the tool, not the one generated
  // registerAllDialects(registry);

  // Stoichiometry matrices are meant to be read: print them as numbers even
  // past MLIR's default hex limit of 100 elements (e.g. 21 reactions x 6
  // species). An explicit --mlir-print-elementsattrs-with-hex-if-larger on the
  // command line comes later and wins.
  std::vector<char *> args(argv, argv + argc);
  char noHex[] = "--mlir-print-elementsattrs-with-hex-if-larger=-1";
  args.insert(args.begin() + 1, noHex);

  return mlir::asMainReturnCode(
      mlir::MlirOptMain(static_cast<int>(args.size()), args.data(),
                        "ChemComp optimizer driver\n", registry));
}
