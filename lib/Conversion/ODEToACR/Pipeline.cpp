//===- Pipeline.cpp - The ode-to-crn pipeline -----------------------------===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#include "ChemComp/Conversion/ODEToACR/Passes.h"

#include "ChemComp/Dialect/ACR/ACRPasses.h"
#include "ChemComp/Dialect/ODE/ODEPasses.h"

#include "mlir/Pass/PassManager.h"
#include "mlir/Pass/PassOptions.h"
#include "mlir/Pass/PassRegistry.h"

using namespace mlir;

namespace {

struct ODEToCRNOptions : public PassPipelineOptions<ODEToCRNOptions> {
  Option<std::string> dualRail{
      *this, "dual-rail",
      llvm::cl::desc("Dual-rail encoding: 'auto' (only states that need it), "
                     "'all', or 'none'"),
      llvm::cl::init("auto")};
  Option<double> annihilationRate{
      *this, "annihilation-rate",
      llvm::cl::desc("Rate of the annihilation x_p + x_m -> 0"),
      llvm::cl::init(1.0)};
  Option<bool> canonicalFallback{
      *this, "canonical-fallback",
      llvm::cl::desc("Cover terms no rule claims with the canonical "
                     "Hungarian-Lemma reaction"),
      llvm::cl::init(true)};
  Option<bool> merge{*this, "merge", llvm::cl::desc("Run acr-network-merge"),
                     llvm::cl::init(true)};
};

void buildODEToCRNPipeline(OpPassManager &pm, const ODEToCRNOptions &options) {
  pm.addPass(cc::ode::createODENormalize());
  if (options.dualRail != "none") {
    cc::ode::ODEDualRailOptions dualRail;
    dualRail.mode = options.dualRail;
    dualRail.annihilationRate = options.annihilationRate;
    pm.addPass(cc::ode::createODEDualRail(dualRail));
  }
  cc::ConvertODEToACROptions convert;
  convert.canonicalFallback = options.canonicalFallback;
  pm.addPass(cc::createConvertODEToACR(convert));
  if (options.merge)
    pm.addPass(cc::acr::createACRNetworkMerge());
  pm.addPass(cc::acr::createACRNetworkToCRN());
}

} // namespace

void mlir::cc::registerODEToCRNPipeline() {
  PassPipelineRegistration<ODEToCRNOptions>(
      "ode-to-crn",
      "Compile ode.system to acr.crn: normalize, dual-rail where needed, "
      "compose reactions from rules (canonical fallback), merge, lower",
      buildODEToCRNPipeline);
}
