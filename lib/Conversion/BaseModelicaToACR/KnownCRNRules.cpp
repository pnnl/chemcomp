//===- KnownCRNRules.cpp - ODE systems of known reaction networks ---------===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#include "KnownCRNRules.h"

#include <algorithm>
#include <cmath>
#include <numeric>

using namespace mlir::cc::ode;

//===----------------------------------------------------------------------===//
// Rule table
//===----------------------------------------------------------------------===//

llvm::ArrayRef<KnownCRNRule> mlir::cc::ode::getKnownCRNRules() {
  // Stoichiometry vectors are indexed by speciesNames.
  static const std::vector<KnownCRNRule> rules = {
      // der(x) = -k*x
      {"decay", {"X"}, {{/*X ->*/ {1}, /*0*/ {0}}}},

      // der(x) = k1 - k2*x
      {"birth-death", {"X"}, {{/*0 ->*/ {0}, /*X*/ {1}}, {{1}, {0}}}},

      // der(a) = -k1*a + k2*b, der(b) = k1*a - k2*b
      {"reversible-isomerization",
       {"A", "B"},
       {{/*A ->*/ {1, 0}, /*B*/ {0, 1}}, {/*B ->*/ {0, 1}, /*A*/ {1, 0}}}},

      // der(x) = k1*x - k2*x*y, der(y) = k2*x*y - k3*y
      {"lotka-volterra",
       {"X", "Y"},
       {{/*X ->*/ {1, 0}, /*2X*/ {2, 0}},
        {/*X + Y ->*/ {1, 1}, /*2Y*/ {0, 2}},
        {/*Y ->*/ {0, 1}, /*0*/ {0, 0}}}},
  };
  return rules;
}

//===----------------------------------------------------------------------===//
// Matching
//===----------------------------------------------------------------------===//

std::string RateConstant::str() const {
  std::string result;
  if (magnitude != 1.0 || parameters.empty())
    result = formatNumber(magnitude);
  for (const auto &[parameter, exponent] : parameters) {
    if (!result.empty())
      result += "*";
    result += parameter;
    if (exponent != 1)
      result += "^" + std::to_string(exponent);
  }
  return result;
}

namespace {

/// A term of a rule's mass-action ODE: net * k_reaction * prod x^reactants.
struct RuleTerm {
  std::vector<unsigned> exponents;
  int net;
  size_t reaction;
};

} // namespace

/// Mass-action ODEs of a rule, one term list per species. Returns nullopt if
/// two reactions contribute the same monomial to one species, which the
/// term-by-term matching below cannot disambiguate.
static std::optional<std::vector<std::vector<RuleTerm>>>
getMassActionODEs(const KnownCRNRule &rule) {
  size_t numSpecies = rule.speciesNames.size();
  std::vector<std::vector<RuleTerm>> odes(numSpecies);
  for (size_t r = 0; r < rule.reactions.size(); ++r) {
    const RuleReaction &reaction = rule.reactions[r];
    for (size_t s = 0; s < numSpecies; ++s) {
      int net = static_cast<int>(reaction.products[s]) -
                static_cast<int>(reaction.reactants[s]);
      if (net == 0)
        continue;
      for (const RuleTerm &other : odes[s])
        if (other.exponents == reaction.reactants)
          return std::nullopt;
      odes[s].push_back({reaction.reactants, net, r});
    }
  }
  return odes;
}

static bool nearlyEqual(double a, double b) {
  return std::fabs(a - b) <= 1e-9 * std::max(std::fabs(a), std::fabs(b));
}

/// Try to match `rule` against `system` with species s mapped to model state
/// stateOfSpecies[s].
static std::optional<std::vector<RateConstant>>
matchWithMapping(const std::vector<std::vector<RuleTerm>> &odes,
                 size_t numReactions, const ODESystem &system,
                 const std::vector<size_t> &stateOfSpecies) {
  std::vector<std::optional<RateConstant>> rates(numReactions);

  for (size_t s = 0; s < odes.size(); ++s) {
    const std::vector<ODETerm> &modelTerms =
        system.equations[stateOfSpecies[s]];
    if (modelTerms.size() != odes[s].size())
      return std::nullopt;

    for (const RuleTerm &ruleTerm : odes[s]) {
      // Rule exponents, reindexed to model states.
      std::vector<unsigned> exponents(system.states.size(), 0);
      for (size_t t = 0; t < ruleTerm.exponents.size(); ++t)
        exponents[stateOfSpecies[t]] = ruleTerm.exponents[t];

      const ODETerm *modelTerm = nullptr;
      for (const ODETerm &candidate : modelTerms) {
        if (candidate.stateExponents != exponents)
          continue;
        if (modelTerm) // Several model terms share this monomial.
          return std::nullopt;
        modelTerm = &candidate;
      }
      if (!modelTerm)
        return std::nullopt;

      // The sign comes from the stoichiometry; the rate constant is positive.
      if ((modelTerm->coefficient < 0) != (ruleTerm.net < 0))
        return std::nullopt;
      RateConstant rate{std::fabs(modelTerm->coefficient) /
                            std::abs(ruleTerm.net),
                        modelTerm->parameters};

      std::optional<RateConstant> &bound = rates[ruleTerm.reaction];
      if (!bound) {
        bound = rate;
      } else if (bound->parameters != rate.parameters ||
                 !nearlyEqual(bound->magnitude, rate.magnitude)) {
        return std::nullopt;
      }
    }
  }

  std::vector<RateConstant> result;
  for (std::optional<RateConstant> &rate : rates) {
    if (!rate)
      return std::nullopt;
    result.push_back(*rate);
  }
  return result;
}

std::optional<CRNMatch> mlir::cc::ode::matchKnownCRN(const ODESystem &system) {
  for (const KnownCRNRule &rule : getKnownCRNRules()) {
    size_t numSpecies = rule.speciesNames.size();
    if (numSpecies != system.states.size())
      continue;
    auto odes = getMassActionODEs(rule);
    if (!odes)
      continue;

    // Try every assignment of rule species to model states; the systems in
    // the table are small.
    std::vector<size_t> stateOfSpecies(numSpecies);
    std::iota(stateOfSpecies.begin(), stateOfSpecies.end(), 0);
    do {
      if (auto rates = matchWithMapping(*odes, rule.reactions.size(), system,
                                        stateOfSpecies))
        return CRNMatch{&rule, stateOfSpecies, std::move(*rates)};
    } while (
        std::next_permutation(stateOfSpecies.begin(), stateOfSpecies.end()));
  }
  return std::nullopt;
}
