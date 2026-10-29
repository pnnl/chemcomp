//===- BaseModelicaToACR.cpp - Base Modelica ODEs to ACR networks ---------===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#include "ChemComp/Conversion/BaseModelicaToACR/Passes.h"

#include "KnownCRNRules.h"
#include "ODEPolynomial.h"

#include "ChemComp/Dialect/ACR/ACROps.h"
#include "ChemComp/Dialect/ACR/ACRTypes.h"

#include "marco/Dialect/BaseModelica/IR/Dialect.h"

#include "mlir/IR/Builders.h"
#include "mlir/IR/BuiltinAttributes.h"
#include "mlir/IR/BuiltinTypes.h"

namespace mlir {
namespace cc {
#define GEN_PASS_DEF_CONVERTBASEMODELICATOACR
#include "ChemComp/Conversion/BaseModelicaToACR/Passes.h.inc"
} // namespace cc
} // namespace mlir

using namespace mlir;
using namespace mlir::cc;
using namespace mlir::cc::ode;
namespace bmodelica = mlir::bmodelica;

namespace {

/// Stoichiometry of the matched network, with species columns in model state
/// order.
struct NetworkMatrices {
  SmallVector<int32_t> reactants;
  SmallVector<int32_t> products;
  SmallVector<int32_t> net;
};

NetworkMatrices buildMatrices(const CRNMatch &match, size_t numStates) {
  const KnownCRNRule &rule = *match.rule;
  size_t numReactions = rule.reactions.size();
  NetworkMatrices m;
  m.reactants.assign(numReactions * numStates, 0);
  m.products.assign(numReactions * numStates, 0);
  m.net.assign(numReactions * numStates, 0);
  for (size_t r = 0; r < numReactions; ++r) {
    for (size_t s = 0; s < rule.speciesNames.size(); ++s) {
      size_t index = r * numStates + match.stateOfSpecies[s];
      m.reactants[index] = rule.reactions[r].reactants[s];
      m.products[index] = rule.reactions[r].products[s];
      m.net[index] = m.products[index] - m.reactants[index];
    }
  }
  return m;
}

/// Emit `acr.crn @<model>` holding the matched network.
void emitNetwork(OpBuilder &builder, bmodelica::ModelOp model, StringRef name,
                 const ODESystem &system, const CRNMatch &match) {
  MLIRContext *context = builder.getContext();
  Location loc = model.getLoc();
  auto numStates = static_cast<int64_t>(system.states.size());
  auto numReactions = static_cast<int64_t>(match.rule->reactions.size());

  auto speciesTensorType =
      RankedTensorType::get({numStates}, acr::SpeciesType::get(context));
  auto rateTensorType =
      RankedTensorType::get({numReactions}, builder.getF32Type());
  auto funcType = builder.getFunctionType({speciesTensorType, rateTensorType},
                                          {speciesTensorType});

  auto crnOp = acr::CRNOp::create(builder, loc, name, funcType,
                                  /*sym_visibility=*/StringAttr{},
                                  /*arg_attrs=*/ArrayAttr{},
                                  /*res_attrs=*/ArrayAttr{});
  Block *body = crnOp.addEntryBlock();
  OpBuilder::InsertionGuard guard(builder);
  builder.setInsertionPointToStart(body);

  auto generic = acr::CRNGenericOp::create(
      builder, loc, TypeRange{speciesTensorType},
      ValueRange{body->getArgument(0), body->getArgument(1)});

  NetworkMatrices matrices = buildMatrices(match, system.states.size());
  auto matrixType =
      RankedTensorType::get({numReactions, numStates}, builder.getI32Type());
  generic->setAttr("stoich_matrix",
                   DenseIntElementsAttr::get(matrixType, matrices.net));
  generic->setAttr("reactant_matrix",
                   DenseIntElementsAttr::get(matrixType, matrices.reactants));
  generic->setAttr("product_matrix",
                   DenseIntElementsAttr::get(matrixType, matrices.products));
  generic->setAttr("crn_rule", builder.getStringAttr(match.rule->name));

  SmallVector<Attribute> stateNames;
  for (const std::string &state : system.states)
    stateNames.push_back(builder.getStringAttr(state));
  generic->setAttr("ode_state_names", builder.getArrayAttr(stateNames));

  SmallVector<Attribute> rateNames;
  SmallVector<double> rateValues;
  bool allRatesKnown = true;
  for (const RateConstant &rate : match.rates) {
    rateNames.push_back(builder.getStringAttr(rate.str()));
    double value = rate.magnitude;
    for (const auto &[parameter, exponent] : rate.parameters) {
      auto it = system.parameterValues.find(parameter);
      if (it == system.parameterValues.end()) {
        allRatesKnown = false;
        break;
      }
      for (unsigned i = 0; i < exponent; ++i)
        value *= it->second;
    }
    rateValues.push_back(value);
  }
  generic->setAttr("rate_names", builder.getArrayAttr(rateNames));
  if (allRatesKnown)
    generic->setAttr("rate_values", builder.getDenseF64ArrayAttr(rateValues));

  acr::ReturnOp::create(builder, loc, generic->getResults());
}

struct ConvertBaseModelicaToACRPass
    : public cc::impl::ConvertBaseModelicaToACRBase<
          ConvertBaseModelicaToACRPass> {
  using ConvertBaseModelicaToACRBase::ConvertBaseModelicaToACRBase;

  void runOnOperation() final {
    ModuleOp module = getOperation();
    OpBuilder builder(module.getContext());
    bool failedAny = false;

    for (bmodelica::ModelOp model :
         llvm::make_early_inc_range(module.getOps<bmodelica::ModelOp>())) {
      FailureOr<ODESystem> system = extractODESystem(model);
      if (failed(system)) {
        failedAny |= !allowUnmatched;
        continue;
      }

      std::optional<CRNMatch> match = matchKnownCRN(*system);
      if (!match) {
        InFlightDiagnostic diag =
            allowUnmatched ? model.emitWarning() : model.emitError();
        diag << "no known chemical reaction network matches the ODE system "
                "of model '"
             << model.getSymName() << "'";
        diag.attachNote() << "ODE system: " << system->str();
        failedAny |= !allowUnmatched;
        continue;
      }

      // The network takes the model's name; a kept model holds on to it.
      std::string name = model.getSymName().str();
      if (keepModel)
        name += "_crn";
      builder.setInsertionPoint(model);
      emitNetwork(builder, model, name, *system, *match);
      if (!keepModel)
        model.erase();
    }

    if (failedAny)
      signalPassFailure();
  }
};

} // namespace
