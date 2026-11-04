//===- Passes.h - ODE to ACR conversion passes ------------------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#ifndef CHEMCOMP_CONVERSION_ODETOACR_PASSES_H
#define CHEMCOMP_CONVERSION_ODETOACR_PASSES_H

#include "ChemComp/Dialect/ACR/ACRDialect.h"
#include "mlir/IR/BuiltinOps.h"
#include "mlir/Pass/Pass.h"

#include <string>
#include <vector>

namespace mlir {
namespace cc {

#define GEN_PASS_DECL
#include "ChemComp/Conversion/ODEToACR/Passes.h.inc"

#define GEN_PASS_REGISTRATION
#include "ChemComp/Conversion/ODEToACR/Passes.h.inc"

/// Register the `ode-to-crn` pipeline: ode-normalize, ode-dual-rail,
/// convert-ode-to-acr, acr-network-merge, acr-network-to-crn.
void registerODEToCRNPipeline();

namespace ode {

/// One composable rule: a reaction over rule-local species.
struct ODERule {
  const char *name;
  int priority;
  /// e.g. "A + B -> C".
  const char *reaction;
  const char *summary;
  std::vector<std::string> species;
  std::vector<unsigned> reactants;
  std::vector<unsigned> products;
};

/// The fixed rule set, by descending priority.
const std::vector<ODERule> &getODERules();

} // namespace ode
} // namespace cc
} // namespace mlir

#endif // CHEMCOMP_CONVERSION_ODETOACR_PASSES_H
