//===- DualRail.cpp - ode-dual-rail ---------------------------------------===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// Dual-rail encoding x = x_p - x_m (e.g. Oishi & Klavins 2011). After the
// substitution every term is routed to the rail where it is a production or a
// self-degradation, so the result has no negative cross-effects and therefore
// is the mass-action system of some network (Hungarian Lemma).
//
//===----------------------------------------------------------------------===//

#include "ChemComp/Dialect/ODE/ODEPasses.h"
#include "ChemComp/Dialect/ODE/ODESystem.h"

#include "mlir/IR/Builders.h"

#include <set>

namespace mlir::cc::ode {
#define GEN_PASS_DEF_ODEDUALRAIL
#include "ChemComp/Dialect/ODE/ODEPasses.h.inc"
} // namespace mlir::cc::ode

using namespace mlir;
using namespace mlir::cc::ode;

namespace {

struct Rails {
  std::string plus;
  std::string minus;
};

/// A name not used by any state or parameter of `system`.
std::string freshName(const System &system, std::string name) {
  auto taken = [&](const std::string &candidate) {
    if (system.stateIndex(candidate) >= 0)
      return true;
    for (const Param &param : system.params)
      if (param.name == candidate)
        return true;
    return false;
  };
  while (taken(name))
    name += "_";
  return name;
}

double binomial(unsigned n, unsigned k) {
  double result = 1;
  for (unsigned i = 1; i <= k; ++i)
    result = result * (n - k + i) / i;
  return result;
}

/// `system` with the states in `split` replaced by their two rails.
System splitStates(const System &system, const std::set<std::string> &split,
                   double eta) {
  std::map<std::string, Rails> rails;
  for (const State &state : system.states)
    if (split.count(state.name))
      rails[state.name] = {freshName(system, state.name + "_p"),
                           freshName(system, state.name + "_m")};

  System result;
  result.name = system.name;
  result.params = system.params;
  for (const State &state : system.states) {
    auto it = rails.find(state.name);
    if (it == rails.end()) {
      result.states.push_back(state);
      continue;
    }
    State plus{it->second.plus, std::nullopt, state.name, 1};
    State minus{it->second.minus, std::nullopt, state.name, -1};
    if (state.initial) {
      plus.initial = *state.initial > 0 ? *state.initial : 0.0;
      minus.initial = *state.initial < 0 ? -*state.initial : 0.0;
    }
    result.states.push_back(plus);
    result.states.push_back(minus);
  }

  for (const Term &term : system.terms) {
    // Expand prod x_j^e with x_j = x_jp - x_jm for the split states.
    std::vector<std::pair<double, Factors>> expanded = {{term.coeff, {}}};
    for (const auto &[name, exponent] : term.states) {
      auto it = rails.find(name);
      std::vector<std::pair<double, Factors>> next;
      for (const auto &[coeff, factors] : expanded) {
        if (it == rails.end()) {
          Factors f = factors;
          f[name] += exponent;
          next.push_back({coeff, f});
          continue;
        }
        for (unsigned k = 0; k <= exponent; ++k) {
          Factors f = factors;
          if (exponent - k)
            f[it->second.plus] += exponent - k;
          if (k)
            f[it->second.minus] += k;
          double sign = (k % 2) ? -1.0 : 1.0;
          next.push_back({coeff * sign * binomial(exponent, k), f});
        }
      }
      expanded = std::move(next);
    }

    auto target = rails.find(term.state);
    for (const auto &[coeff, factors] : expanded) {
      if (target == rails.end()) {
        result.terms.push_back(
            {term.state, coeff, term.params, factors, term.loc});
        continue;
      }
      const Rails &r = target->second;
      // der(x) = der(x_p) - der(x_m): a term can be added to x_p as is, or
      // to x_m with the opposite sign. Pick the rail where it is either a
      // production or a degradation of that rail.
      bool toMinus = coeff > 0 ? factors.count(r.minus) > 0
                               : factors.count(r.plus) == 0;
      result.terms.push_back({toMinus ? r.minus : r.plus,
                              toMinus ? -coeff : coeff, term.params, factors,
                              term.loc});
    }
  }

  // Annihilation x_p + x_m -> 0 removes the same amount from both rails.
  for (const auto &[name, r] : rails) {
    Factors both{{r.plus, 1}, {r.minus, 1}};
    result.terms.push_back({r.plus, -eta, {}, both, std::nullopt});
    result.terms.push_back({r.minus, -eta, {}, both, std::nullopt});
  }

  result.normalize();
  return result;
}

struct ODEDualRailPass
    : public mlir::cc::ode::impl::ODEDualRailBase<ODEDualRailPass> {
  using ODEDualRailBase::ODEDualRailBase;

  void runOnOperation() final {
    if (mode != "auto" && mode != "all") {
      getOperation().emitError("ode-dual-rail: unknown mode '")
          << mode << "'; expected 'auto' or 'all'";
      return signalPassFailure();
    }
    if (annihilationRate <= 0) {
      getOperation().emitError(
          "ode-dual-rail: annihilation-rate must be positive");
      return signalPassFailure();
    }

    WalkResult result = getOperation().walk([&](SystemOp op) {
      return failed(runOnSystem(op)) ? WalkResult::interrupt()
                                     : WalkResult::advance();
    });
    if (result.wasInterrupted())
      signalPassFailure();
  }

  LogicalResult runOnSystem(SystemOp op) {
    if (op->hasAttr("dual_rail"))
      return success(); // Already encoded.

    System system = System::read(op);
    std::set<std::string> split;
    for (const std::string &name : extraStates) {
      if (system.stateIndex(name) < 0)
        return op.emitError("ode-dual-rail: '")
               << name << "' is not a state of system '" << system.name << "'";
      split.insert(name);
    }
    if (mode == "all")
      for (const State &state : system.states)
        split.insert(state.name);

    // Split the targets of cross-effects until none is left. Each round adds
    // at least one state (split rails never have cross-effects), so this ends
    // after at most one round per state.
    System encoded = system;
    while (true) {
      if (!split.empty())
        encoded = splitStates(system, split, annihilationRate);
      std::vector<CrossEffect> effects = findCrossEffects(encoded);
      if (effects.empty())
        break;
      size_t before = split.size();
      for (const CrossEffect &effect : effects) {
        const Term &term = encoded.terms[effect.termIndex];
        split.insert(
            encoded.states[encoded.stateIndex(term.state)].encodedState().str());
      }
      if (split.size() == before)
        return op.emitError("ode-dual-rail: could not remove the negative "
                            "cross-effects of system '")
               << system.name << "'";
    }
    if (split.empty())
      return success();

    encoded.write(op);
    Builder builder(op.getContext());
    SmallVector<Attribute> names;
    for (const State &state : system.states)
      if (split.count(state.name))
        names.push_back(builder.getStringAttr(state.name));
    op->setAttr("dual_rail", builder.getArrayAttr(names));
    return success();
  }
};

} // namespace
