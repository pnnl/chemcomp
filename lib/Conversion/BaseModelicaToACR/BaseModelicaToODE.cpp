//===- BaseModelicaToODE.cpp - Base Modelica ODEs to ode.system -----------===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#include "ChemComp/Conversion/BaseModelicaToACR/Passes.h"

#include "ODEPolynomial.h"

#include "ChemComp/Dialect/ODE/ODEOps.h"
#include "ChemComp/Dialect/ODE/ODESystem.h"

#include "marco/Dialect/BaseModelica/IR/Dialect.h"

#include "mlir/IR/Builders.h"

#include <set>

namespace mlir {
namespace cc {
#define GEN_PASS_DEF_CONVERTBASEMODELICATOODE
#include "ChemComp/Conversion/BaseModelicaToACR/Passes.h.inc"
} // namespace cc
} // namespace mlir

using namespace mlir;
using namespace mlir::cc;
namespace bmodelica = mlir::bmodelica;

namespace {

struct ConvertBaseModelicaToODEPass
    : public cc::impl::ConvertBaseModelicaToODEBase<
          ConvertBaseModelicaToODEPass> {
  using ConvertBaseModelicaToODEBase::ConvertBaseModelicaToODEBase;

  void runOnOperation() final {
    bool failedAny = false;
    for (bmodelica::ModelOp model : llvm::make_early_inc_range(
             getOperation().getOps<bmodelica::ModelOp>())) {
      FailureOr<ode::ODESystem> extracted = ode::extractODESystem(model);
      if (failed(extracted) || failed(emit(model, *extracted))) {
        failedAny = true;
        continue;
      }
      if (!keepModel)
        model.erase();
    }
    if (failedAny)
      signalPassFailure();
  }

  LogicalResult emit(bmodelica::ModelOp model, const ode::ODESystem &in) {
    ode::System out;
    out.name = model.getSymName().str();
    if (keepModel)
      out.name += "_ode";

    for (const std::string &state : in.states) {
      ode::State s{state, std::nullopt, std::nullopt, std::nullopt};
      if (auto it = in.initialValues.find(state); it != in.initialValues.end())
        s.initial = it->second;
      out.states.push_back(std::move(s));
    }

    std::set<std::string> params;
    for (size_t i = 0; i < in.states.size(); ++i) {
      for (const ode::ODETerm &term : in.equations[i]) {
        ode::Term t{in.states[i], term.coefficient, term.parameters, {},
                    std::nullopt};
        for (size_t j = 0; j < in.states.size(); ++j)
          if (term.stateExponents[j])
            t.states[in.states[j]] = term.stateExponents[j];
        for (const auto &[param, exponent] : term.parameters)
          params.insert(param);
        out.terms.push_back(std::move(t));
      }
    }
    for (const std::string &param : params) {
      ode::Param p{param, std::nullopt};
      if (auto it = in.parameterValues.find(param);
          it != in.parameterValues.end()) {
        if (it->second < 0)
          return model.emitError()
                 << "parameter '" << param << "' is negative (" << it->second
                 << "); rate constants must not be negative, move the sign "
                    "into the equation";
        p.value = it->second;
      }
      out.params.push_back(std::move(p));
    }
    out.normalize();

    OpBuilder builder(model);
    auto system = ode::SystemOp::create(builder, model.getLoc(), out.name);
    out.write(system);
    return success();
  }
};

} // namespace
