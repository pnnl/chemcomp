//===- NetworkPasses.cpp - acr-network-merge and acr-network-to-crn -------===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#include "ChemComp/Dialect/ACR/ACRPasses.h"
#include "ChemComp/Dialect/ACR/ACRTypes.h"

#include "mlir/IR/Builders.h"
#include "mlir/IR/BuiltinAttributes.h"
#include "mlir/IR/BuiltinTypes.h"
#include "llvm/ADT/STLExtras.h"

#include <cmath>
#include <cstdio>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace mlir::cc::acr {
#define GEN_PASS_DEF_ACRNETWORKMERGE
#define GEN_PASS_DEF_ACRNETWORKTOCRN
#include "ChemComp/Dialect/ACR/ACRPasses.h.inc"
} // namespace mlir::cc::acr

using namespace mlir;
using namespace mlir::cc::acr;

namespace {

using Complex = std::map<std::string, int64_t>;

Complex readComplex(DictionaryAttr dict) {
  Complex complex;
  for (NamedAttribute entry : dict)
    complex[entry.getName().str()] =
        cast<IntegerAttr>(entry.getValue()).getInt();
  return complex;
}

DictionaryAttr writeComplex(Builder &builder, const Complex &complex) {
  SmallVector<NamedAttribute> entries;
  for (const auto &[name, count] : complex)
    if (count)
      entries.push_back(
          builder.getNamedAttr(name, builder.getI64IntegerAttr(count)));
  return builder.getDictionaryAttr(entries);
}

std::string formatRate(double value) {
  char buffer[32];
  std::snprintf(buffer, sizeof(buffer), "%g", value);
  return buffer;
}

/// A plain copy of an acr.rxn.
struct Reaction {
  Complex reactants;
  Complex products;
  std::optional<double> rate;
  std::string rateName;
  std::string rule;
  std::string from;
  Location loc;

  static Reaction read(RxnOp op) {
    Reaction r{readComplex(op.getReactants()),
               readComplex(op.getProducts()),
               std::nullopt,
               op.getRateName().value_or("").str(),
               op.getRule().value_or("").str(),
               op.getFrom().value_or("").str(),
               op.getLoc()};
    if (FloatAttr rate = op.getRateAttr())
      r.rate = rate.getValueAsDouble();
    return r;
  }

  void write(OpBuilder &builder) const {
    auto op = RxnOp::create(builder, loc, writeComplex(builder, reactants),
                            writeComplex(builder, products), nullptr, nullptr,
                            nullptr, nullptr);
    if (rate)
      op.setRateAttr(builder.getF64FloatAttr(*rate));
    if (!rateName.empty())
      op.setRateNameAttr(builder.getStringAttr(rateName));
    if (!rule.empty())
      op.setRuleAttr(builder.getStringAttr(rule));
    if (!from.empty())
      op.setFromAttr(builder.getStringAttr(from));
  }
};

std::string join(const std::string &a, const std::string &b,
                 const char *separator) {
  if (a.empty())
    return b;
  if (b.empty() || a == b)
    return a;
  return a + separator + b;
}

bool sameRate(const Reaction &a, const Reaction &b) {
  if (a.rateName != b.rateName)
    return false;
  if (a.rate && b.rate)
    return std::fabs(*a.rate - *b.rate) <=
           1e-12 * std::max(std::fabs(*a.rate), std::fabs(*b.rate));
  return !a.rate && !b.rate && !a.rateName.empty();
}

/// Merge `b` into `a` if that keeps the mass-action ODEs; see the pass
/// description for the two cases.
bool tryMerge(Reaction &a, const Reaction &b) {
  if (a.reactants != b.reactants)
    return false;

  if (a.products == b.products) {
    if (a.rate && b.rate)
      a.rate = *a.rate + *b.rate;
    else
      a.rate.reset();
    a.rateName = a.rateName.empty() || b.rateName.empty()
                     ? std::string()
                     : a.rateName + " + " + b.rateName;
    a.rule = join(a.rule, b.rule, "+");
    a.from = join(a.from, b.from, "; ");
    return true;
  }

  if (!sameRate(a, b))
    return false;
  // y -> y + v1 and y -> y + v2 at the same rate act like y -> y + v1 + v2.
  Complex merged = a.products;
  for (const auto &[name, count] : b.products)
    merged[name] += count;
  for (const auto &[name, count] : a.reactants)
    merged[name] -= count;
  for (const auto &[name, count] : merged)
    if (count < 0)
      return false;
  for (auto it = merged.begin(); it != merged.end();)
    it = it->second == 0 ? merged.erase(it) : std::next(it);
  a.products = merged;
  a.rule = join(a.rule, b.rule, "+");
  a.from = join(a.from, b.from, "; ");
  return true;
}

struct ACRNetworkMergePass
    : public mlir::cc::acr::impl::ACRNetworkMergeBase<ACRNetworkMergePass> {
  void runOnOperation() final {
    getOperation().walk([](NetworkOp network) {
      std::vector<Reaction> reactions;
      for (RxnOp rxn : network.getBody().getOps<RxnOp>())
        reactions.push_back(Reaction::read(rxn));

      bool changed = true;
      while (changed) {
        changed = false;
        llvm::erase_if(reactions, [](const Reaction &r) {
          return r.reactants == r.products;
        });
        for (size_t i = 0; i < reactions.size() && !changed; ++i) {
          for (size_t j = i + 1; j < reactions.size(); ++j) {
            if (tryMerge(reactions[i], reactions[j])) {
              reactions.erase(reactions.begin() + j);
              changed = true;
              break;
            }
          }
        }
      }

      for (RxnOp rxn :
           llvm::make_early_inc_range(network.getBody().getOps<RxnOp>()))
        rxn.erase();
      OpBuilder builder = OpBuilder::atBlockEnd(&network.getBody().front());
      for (const Reaction &reaction : reactions)
        reaction.write(builder);
    });
  }
};

struct ACRNetworkToCRNPass
    : public mlir::cc::acr::impl::ACRNetworkToCRNBase<ACRNetworkToCRNPass> {
  void runOnOperation() final {
    WalkResult result = getOperation().walk([](NetworkOp network) {
      return failed(lower(network)) ? WalkResult::interrupt()
                                    : WalkResult::advance();
    });
    if (result.wasInterrupted())
      signalPassFailure();
  }

  static LogicalResult lower(NetworkOp network) {
    SmallVector<SpeciesOp> species(network.getBody().getOps<SpeciesOp>());
    SmallVector<RxnOp> reactions(network.getBody().getOps<RxnOp>());
    if (species.empty() || reactions.empty())
      return network.emitError("acr-network-to-crn: network '")
             << network.getSymName() << "' needs at least one species and "
             << "one reaction";

    auto numSpecies = static_cast<int64_t>(species.size());
    auto numReactions = static_cast<int64_t>(reactions.size());
    std::map<std::string, int64_t> column;
    for (int64_t s = 0; s < numSpecies; ++s)
      column[species[s].getSymName().str()] = s;

    SmallVector<int32_t> net(numReactions * numSpecies, 0),
        reactantMatrix(numReactions * numSpecies, 0),
        productMatrix(numReactions * numSpecies, 0);
    SmallVector<Attribute> rateNames, rules;
    SmallVector<double> rates;
    bool allRates = true;
    OpBuilder builder(network);
    for (int64_t r = 0; r < numReactions; ++r) {
      RxnOp rxn = reactions[r];
      for (NamedAttribute entry : rxn.getReactants())
        reactantMatrix[r * numSpecies + column[entry.getName().str()]] =
            cast<IntegerAttr>(entry.getValue()).getInt();
      for (NamedAttribute entry : rxn.getProducts())
        productMatrix[r * numSpecies + column[entry.getName().str()]] =
            cast<IntegerAttr>(entry.getValue()).getInt();
      FloatAttr rate = rxn.getRateAttr();
      std::string name = rxn.getRateName().value_or("").str();
      if (name.empty() && rate)
        name = formatRate(rate.getValueAsDouble());
      rateNames.push_back(builder.getStringAttr(name));
      rules.push_back(builder.getStringAttr(rxn.getRule().value_or("")));
      if (rate)
        rates.push_back(rate.getValueAsDouble());
      else
        allRates = false;
    }
    for (int64_t i = 0; i < numReactions * numSpecies; ++i)
      net[i] = productMatrix[i] - reactantMatrix[i];

    SmallVector<Attribute> symbols, stateNames;
    SmallVector<int64_t> signs;
    SmallVector<double> initial;
    bool allInitial = true;
    for (SpeciesOp sp : species) {
      symbols.push_back(builder.getStringAttr(sp.getSymName()));
      stateNames.push_back(
          builder.getStringAttr(sp.getOdeState().value_or(sp.getSymName())));
      signs.push_back(sp.getRailSign() ? static_cast<int64_t>(*sp.getRailSign())
                                       : 1);
      if (FloatAttr value = sp.getInitialAttr())
        initial.push_back(value.getValueAsDouble());
      else
        allInitial = false;
    }

    MLIRContext *context = network.getContext();
    Location loc = network.getLoc();
    auto speciesType =
        RankedTensorType::get({numSpecies}, SpeciesType::get(context));
    auto rateType = RankedTensorType::get({numReactions}, builder.getF32Type());
    auto crn = CRNOp::create(
        builder, loc, network.getSymName(),
        builder.getFunctionType({speciesType, rateType}, {speciesType}),
        /*sym_visibility=*/StringAttr{}, /*arg_attrs=*/ArrayAttr{},
        /*res_attrs=*/ArrayAttr{});
    Block *body = crn.addEntryBlock();
    OpBuilder::InsertionGuard guard(builder);
    builder.setInsertionPointToStart(body);
    auto generic = CRNGenericOp::create(
        builder, loc, TypeRange{speciesType},
        ValueRange{body->getArgument(0), body->getArgument(1)});

    auto matrixType =
        RankedTensorType::get({numReactions, numSpecies}, builder.getI32Type());
    generic->setAttr("stoich_matrix",
                     DenseIntElementsAttr::get(matrixType, net));
    generic->setAttr("reactant_matrix",
                     DenseIntElementsAttr::get(matrixType, reactantMatrix));
    generic->setAttr("product_matrix",
                     DenseIntElementsAttr::get(matrixType, productMatrix));
    generic->setAttr("species_symbols", builder.getArrayAttr(symbols));
    generic->setAttr("ode_state_names", builder.getArrayAttr(stateNames));
    generic->setAttr("ode_state_signs", builder.getDenseI64ArrayAttr(signs));
    generic->setAttr("rate_names", builder.getArrayAttr(rateNames));
    generic->setAttr("rules", builder.getArrayAttr(rules));
    if (allRates)
      generic->setAttr("rate_values", builder.getDenseF64ArrayAttr(rates));
    if (allInitial)
      generic->setAttr("initial_values", builder.getDenseF64ArrayAttr(initial));
    ReturnOp::create(builder, loc, generic->getResults());

    network.erase();
    return success();
  }
};

} // namespace
