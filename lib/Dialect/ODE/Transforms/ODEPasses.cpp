//===- ODEPasses.cpp - ode-normalize and ode-check-kinetic ----------------===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#include "ChemComp/Dialect/ODE/ODEPasses.h"
#include "ChemComp/Dialect/ODE/ODESystem.h"

namespace mlir::cc::ode {
#define GEN_PASS_DEF_ODENORMALIZE
#define GEN_PASS_DEF_ODECHECKKINETIC
#include "ChemComp/Dialect/ODE/ODEPasses.h.inc"
} // namespace mlir::cc::ode

using namespace mlir;
using namespace mlir::cc::ode;

namespace {

struct ODENormalizePass
    : public mlir::cc::ode::impl::ODENormalizeBase<ODENormalizePass> {
  void runOnOperation() final {
    getOperation().walk([](SystemOp op) {
      System system = System::read(op);
      system.normalize();
      system.write(op);
    });
  }
};

struct ODECheckKineticPass
    : public mlir::cc::ode::impl::ODECheckKineticBase<ODECheckKineticPass> {
  using ODECheckKineticBase::ODECheckKineticBase;

  void runOnOperation() final {
    bool failedAny = false;
    getOperation().walk([&](SystemOp op) {
      System system = System::read(op);
      std::vector<CrossEffect> effects = findCrossEffects(system);
      if (effects.empty()) {
        op.emitRemark() << "system '" << system.name
                        << "' is kinetic: every negative term of der(x) "
                           "contains x";
        return;
      }
      // Point each diagnostic at the offending term.
      SmallVector<TermOp> termOps(op.getBody().getOps<TermOp>());
      for (const CrossEffect &effect : effects) {
        Operation *at = termOps[effect.termIndex];
        InFlightDiagnostic diag = require ? at->emitError() : at->emitRemark();
        diag << "negative cross-effect: " << describe(system, effect);
      }
      failedAny |= require;
    });
    if (failedAny)
      signalPassFailure();
  }
};

} // namespace
