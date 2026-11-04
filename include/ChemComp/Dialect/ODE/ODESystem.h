//===- ODESystem.h - In-memory view of an ode.system ------------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// Passes on ODE systems mostly rewrite whole sets of terms (normalizing,
// substituting x = x_p - x_m, claiming terms for reactions). They read an
// `ode.system` into this plain C++ form, transform it, and write it back.
//
// The Hungarian Lemma (Hars & Toth 1981) lives here too: a polynomial ODE
// system is the mass-action system of some reaction network iff no term has a
// negative cross-effect, i.e. every negative term in der(x) contains x.
//
//===----------------------------------------------------------------------===//

#ifndef CHEMCOMP_DIALECT_ODE_ODESYSTEM_H
#define CHEMCOMP_DIALECT_ODE_ODESYSTEM_H

#include "mlir/IR/OpImplementation.h"
#include "mlir/Support/LLVM.h"
#include "llvm/ADT/StringRef.h"

#include <map>
#include <optional>
#include <string>
#include <vector>

namespace mlir {
namespace cc {
namespace ode {

class SystemOp;

/// Format a number compactly, as printf "%g" does ("2", "0.5", "1e-06").
std::string formatNumber(double value);

/// Shortest decimal spelling that reads back as exactly `value`.
std::string formatExact(double value);

/// Print `value` with formatExact; parse it back (integers are accepted).
void printNumber(OpAsmPrinter &printer, double value);
ParseResult parseNumber(OpAsmParser &parser, double &value);

/// Product of symbols raised to positive powers.
using Factors = std::map<std::string, unsigned>;

/// "b*x*y^2", or "" for an empty product.
std::string factorsString(const Factors &factors);

/// One monomial `coeff * params * states` of der(state).
struct Term {
  std::string state;
  double coeff = 0.0;
  Factors params;
  Factors states;
  /// Where the term came from, for diagnostics; unset for derived terms.
  std::optional<Location> loc;

  /// The monomial alone, e.g. "-b*x*y", "2", "-x".
  std::string str() const;
};

struct State {
  std::string name;
  std::optional<double> initial;
  /// Set on the rails of a dual-rail pair: the encoded state and the sign.
  std::optional<std::string> railOf;
  std::optional<int64_t> railSign;

  /// The original state this one stands for.
  llvm::StringRef encodedState() const {
    return railOf ? llvm::StringRef(*railOf) : llvm::StringRef(name);
  }
  int64_t sign() const { return railSign.value_or(1); }
};

struct Param {
  std::string name;
  std::optional<double> value;
};

/// A plain copy of an `ode.system`.
struct System {
  std::string name;
  std::vector<State> states;
  std::vector<Param> params;
  std::vector<Term> terms;

  static System read(SystemOp op);
  /// Replace the body of `op` by this system.
  void write(SystemOp op) const;

  /// Index of a state in `states`, or -1.
  int stateIndex(llvm::StringRef name) const;
  /// Value of a parameter, if it has one.
  std::optional<double> paramValue(llvm::StringRef name) const;

  /// Sum like terms, drop zero terms, and order terms by state, then by
  /// monomial. Equal systems normalize to equal term lists.
  void normalize();

  /// "der(x) = -k*x; der(y) = ..."
  std::string str() const;
};

/// A negative term of der(x) that does not contain x. Such a term would need
/// a reaction consuming x at a rate independent of x, which mass action
/// cannot express (Hungarian Lemma).
struct CrossEffect {
  size_t termIndex;
};

/// All negative cross-effects of `system`; empty iff it is kinetic.
std::vector<CrossEffect> findCrossEffects(const System &system);

/// "der(y) has -x, which does not contain y".
std::string describe(const System &system, const CrossEffect &effect);

} // namespace ode
} // namespace cc
} // namespace mlir

#endif // CHEMCOMP_DIALECT_ODE_ODESYSTEM_H
