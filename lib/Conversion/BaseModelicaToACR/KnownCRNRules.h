//===- KnownCRNRules.h - ODE systems of known reaction networks -*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// A fixed table of chemical reaction networks. Each rule is given as reactions
// only; its ODE system is derived by mass-action kinetics, so a rule and the
// network it emits cannot disagree.
//
// Matching is a full match: an ODE system matches a rule when some renaming of
// the states and some assignment of rate constants turn the rule's mass-action
// ODEs into exactly the model's ODEs.
//
//===----------------------------------------------------------------------===//

#ifndef CHEMCOMP_LIB_CONVERSION_BASEMODELICATOACR_KNOWNCRNRULES_H
#define CHEMCOMP_LIB_CONVERSION_BASEMODELICATOACR_KNOWNCRNRULES_H

#include "ODEPolynomial.h"

#include "llvm/ADT/ArrayRef.h"

#include <optional>
#include <string>
#include <vector>

namespace mlir {
namespace cc {
namespace ode {

/// A mass-action reaction with stoichiometry indexed by rule species.
struct RuleReaction {
  std::vector<unsigned> reactants;
  std::vector<unsigned> products;
};

/// A known chemical reaction network.
struct KnownCRNRule {
  std::string name;
  std::vector<std::string> speciesNames;
  std::vector<RuleReaction> reactions;
};

/// The rate constant of one reaction: magnitude * product of parameters.
struct RateConstant {
  double magnitude;
  Monomial parameters;

  std::string str() const;
};

/// A rule matched against an ODE system.
struct CRNMatch {
  const KnownCRNRule *rule;
  /// Model state for each rule species.
  std::vector<size_t> stateOfSpecies;
  /// Rate constant of each rule reaction.
  std::vector<RateConstant> rates;
};

/// The table of known networks.
llvm::ArrayRef<KnownCRNRule> getKnownCRNRules();

/// Find the first rule whose mass-action ODEs equal `system`.
std::optional<CRNMatch> matchKnownCRN(const ODESystem &system);

} // namespace ode
} // namespace cc
} // namespace mlir

#endif // CHEMCOMP_LIB_CONVERSION_BASEMODELICATOACR_KNOWNCRNRULES_H
