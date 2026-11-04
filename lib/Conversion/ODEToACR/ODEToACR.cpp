//===- ODEToACR.cpp - Compose reaction networks from ODE terms ------------===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// Under mass action, a reaction y -> y' with rate k adds
//
//   k * (y'_s - y_s) * x^y
//
// to der(s) for every species s. Composing a network for an ODE system means
// covering its terms with such contributions:
//
//   - rules (ODERules.td) are reactions over rule-local species; a rule
//     matches when some assignment of its species to states and some k make
//     all its contributions appear among the unclaimed terms;
//   - the canonical rule (Hars & Toth) covers any single term of a kinetic
//     system: c*x^y in der(x_i) is the contribution of y -> y + sign(c)*X_i.
//
//===----------------------------------------------------------------------===//

#include "ChemComp/Conversion/ODEToACR/Passes.h"

#include "ChemComp/Dialect/ACR/ACROps.h"
#include "ChemComp/Dialect/ODE/ODEOps.h"
#include "ChemComp/Dialect/ODE/ODESystem.h"

#include "mlir/IR/Builders.h"

#include <cmath>

namespace mlir {
namespace cc {
#define GEN_PASS_DEF_CONVERTODETOACR
#include "ChemComp/Conversion/ODEToACR/Passes.h.inc"
} // namespace cc
} // namespace mlir

using namespace mlir;
using namespace mlir::cc;
using namespace mlir::cc::ode;

const std::vector<ODERule> &mlir::cc::ode::getODERules() {
  static const std::vector<ODERule> rules = {
#include "ChemComp/Conversion/ODEToACR/ODERules.inc"
  };
  return rules;
}

namespace {

/// A reaction to emit, over state names.
struct Reaction {
  Factors reactants;
  Factors products;
  /// Rate constant: magnitude * product of parameters.
  double magnitude;
  Factors params;
  std::string rule;
  std::string from;
};

bool approxEqual(double a, double b) {
  return std::fabs(a - b) <= 1e-9 * std::max(std::fabs(a), std::fabs(b));
}

/// "b", "0.5*k", "2".
std::string rateName(double magnitude, const Factors &params) {
  std::string result;
  if (magnitude != 1.0 || params.empty())
    result = formatNumber(magnitude);
  std::string product = factorsString(params);
  if (!result.empty() && !product.empty())
    result += "*";
  return result + product;
}

class Composer {
public:
  explicit Composer(const System &system)
      : system(system), claimed(system.terms.size(), false) {}

  /// Apply `rule` as many times as it matches.
  void applyRule(const ODERule &rule) {
    std::vector<int> assignment;
    while (findMatch(rule, assignment)) {
    }
  }

  /// One reaction per unclaimed term (canonical Hungarian-Lemma rule).
  void applyCanonical() {
    for (size_t i = 0; i < system.terms.size(); ++i) {
      if (claimed[i])
        continue;
      const Term &term = system.terms[i];
      Factors products = term.states;
      if (term.coeff > 0) {
        ++products[term.state];
      } else if (--products[term.state] == 0) {
        products.erase(term.state);
      }
      claimed[i] = true;
      reactions.push_back({term.states, products, std::fabs(term.coeff),
                           term.params, "canonical",
                           "der(" + term.state + "): " + term.str()});
    }
  }

  std::vector<size_t> unclaimed() const {
    std::vector<size_t> result;
    for (size_t i = 0; i < claimed.size(); ++i)
      if (!claimed[i])
        result.push_back(i);
    return result;
  }

  std::vector<Reaction> reactions;

private:
  /// Try every injective assignment of rule species to states, in
  /// lexicographic order; claim and emit the first match.
  bool findMatch(const ODERule &rule, std::vector<int> &assignment) {
    if (assignment.size() == rule.species.size())
      return tryAssignment(rule, assignment);
    for (int s = 0; s < static_cast<int>(system.states.size()); ++s) {
      if (llvm::is_contained(assignment, s))
        continue;
      assignment.push_back(s);
      bool found = findMatch(rule, assignment);
      assignment.pop_back();
      if (found)
        return true;
    }
    return false;
  }

  /// An unclaimed term of der(state) with the given monomial and, optionally,
  /// parameters and coefficient.
  std::optional<size_t> findTerm(const std::string &state,
                                 const Factors &monomial, const Factors *params,
                                 double coeff, size_t from = 0) const {
    for (size_t i = from; i < system.terms.size(); ++i) {
      const Term &term = system.terms[i];
      if (claimed[i] || term.state != state || term.states != monomial)
        continue;
      if (params && (term.params != *params || !approxEqual(term.coeff, coeff)))
        continue;
      return i;
    }
    return std::nullopt;
  }

  bool tryAssignment(const ODERule &rule, const std::vector<int> &assignment) {
    auto stateOf = [&](size_t s) -> const std::string & {
      return system.states[assignment[s]].name;
    };
    Factors monomial, products;
    std::vector<int> change(rule.species.size());
    for (size_t s = 0; s < rule.species.size(); ++s) {
      if (rule.reactants[s])
        monomial[stateOf(s)] = rule.reactants[s];
      if (rule.products[s])
        products[stateOf(s)] = rule.products[s];
      change[s] = static_cast<int>(rule.products[s]) -
                  static_cast<int>(rule.reactants[s]);
    }

    // The first changing species fixes the rate constant k and its
    // parameters; every other changing species must then have its term.
    size_t first = 0;
    while (change[first] == 0)
      ++first;
    for (std::optional<size_t> lead =
             findTerm(stateOf(first), monomial, nullptr, 0);
         lead;
         lead = findTerm(stateOf(first), monomial, nullptr, 0, *lead + 1)) {
      const Term &leadTerm = system.terms[*lead];
      double k = leadTerm.coeff / change[first];
      if (k <= 0)
        continue;

      std::vector<size_t> terms = {*lead};
      bool complete = true;
      for (size_t s = first + 1; s < rule.species.size() && complete; ++s) {
        if (change[s] == 0)
          continue;
        std::optional<size_t> other =
            findTerm(stateOf(s), monomial, &leadTerm.params, change[s] * k);
        if (other)
          terms.push_back(*other);
        else
          complete = false;
      }
      if (!complete)
        continue;

      std::string from;
      for (size_t i : terms) {
        claimed[i] = true;
        const Term &term = system.terms[i];
        if (!from.empty())
          from += "; ";
        from += "der(" + term.state + "): " + term.str();
      }
      reactions.push_back(
          {monomial, products, k, leadTerm.params, rule.name, from});
      return true;
    }
    return false;
  }

  const System &system;
  std::vector<bool> claimed;
};

DictionaryAttr toComplex(Builder &builder, const Factors &factors) {
  SmallVector<NamedAttribute> entries;
  for (const auto &[name, count] : factors)
    entries.push_back(builder.getNamedAttr(
        name, builder.getI64IntegerAttr(static_cast<int64_t>(count))));
  return builder.getDictionaryAttr(entries);
}

struct ConvertODEToACRPass
    : public cc::impl::ConvertODEToACRBase<ConvertODEToACRPass> {
  using ConvertODEToACRBase::ConvertODEToACRBase;

  void runOnOperation() final {
    bool failedAny = false;
    for (SystemOp op : llvm::make_early_inc_range(
             getOperation().getBody()->getOps<SystemOp>()))
      failedAny |= failed(convert(op));
    if (failedAny)
      signalPassFailure();
  }

  LogicalResult convert(SystemOp op) {
    System system = System::read(op);
    system.normalize();

    // Hungarian Lemma: without this, no network can have these ODEs.
    std::vector<CrossEffect> effects = findCrossEffects(system);
    if (!effects.empty()) {
      InFlightDiagnostic diag = op.emitError()
                                << "system '" << system.name
                                << "' has no mass-action network: ";
      llvm::interleaveComma(effects, diag, [&](const CrossEffect &effect) {
        diag << describe(system, effect);
      });
      diag.attachNote() << "run ode-dual-rail first (the ode-to-crn pipeline "
                           "does), or rewrite the model";
      return failure();
    }

    Composer composer(system);
    for (const ODERule &rule : getODERules())
      composer.applyRule(rule);

    if (!canonicalFallback) {
      std::vector<size_t> leftover = composer.unclaimed();
      if (!leftover.empty()) {
        InFlightDiagnostic diag =
            op.emitError() << "no composable rule covers " << leftover.size()
                           << " term(s) of system '" << system.name << "': ";
        llvm::interleave(
            leftover, diag,
            [&](size_t i) {
              const Term &term = system.terms[i];
              diag << "der(" << term.state << "): " << term.str();
            },
            "; ");
        diag.attachNote() << "enable canonical-fallback to cover them with "
                             "one reaction per term";
        return failure();
      }
    }
    composer.applyCanonical();

    OpBuilder builder(op);
    auto network = acr::NetworkOp::create(builder, op.getLoc(), system.name);
    builder.setInsertionPointToStart(&network.getBody().front());
    for (const State &state : system.states) {
      auto species = acr::SpeciesOp::create(builder, op.getLoc(), state.name,
                                            nullptr, nullptr, nullptr);
      if (state.initial)
        species.setInitialAttr(builder.getF64FloatAttr(*state.initial));
      if (state.railOf) {
        species.setOdeStateAttr(builder.getStringAttr(*state.railOf));
        species.setRailSignAttr(builder.getI64IntegerAttr(state.sign()));
      }
    }
    for (const Reaction &reaction : composer.reactions) {
      auto rxn = acr::RxnOp::create(builder, op.getLoc(),
                                    toComplex(builder, reaction.reactants),
                                    toComplex(builder, reaction.products),
                                    nullptr, nullptr, nullptr, nullptr);
      std::optional<double> rate = reaction.magnitude;
      for (const auto &[param, exponent] : reaction.params) {
        std::optional<double> value = system.paramValue(param);
        if (!value) {
          rate.reset();
          break;
        }
        *rate *= std::pow(*value, exponent);
      }
      if (rate)
        rxn.setRateAttr(builder.getF64FloatAttr(*rate));
      rxn.setRateNameAttr(
          builder.getStringAttr(rateName(reaction.magnitude, reaction.params)));
      rxn.setRuleAttr(builder.getStringAttr(reaction.rule));
      rxn.setFromAttr(builder.getStringAttr(reaction.from));
    }
    op.erase();
    return success();
  }
};

} // namespace
