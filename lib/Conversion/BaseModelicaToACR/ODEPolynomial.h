//===- ODEPolynomial.h - Polynomial ODE systems from Base Modelica -*- C++ -*-//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// Normalizes the equations of a bmodelica.model into a polynomial ODE system
//
//   x_i' = sum_t c_t * prod_p p^e_p * prod_j x_j^e_j
//
// where x are the states, p the parameters and c numeric constants. Matching
// against known networks works on this canonical form, so it does not depend
// on operand order or on how the expression was written.
//
//===----------------------------------------------------------------------===//

#ifndef CHEMCOMP_LIB_CONVERSION_BASEMODELICATOACR_ODEPOLYNOMIAL_H
#define CHEMCOMP_LIB_CONVERSION_BASEMODELICATOACR_ODEPOLYNOMIAL_H

#include "mlir/IR/Location.h"
#include "mlir/Support/LLVM.h"
#include "llvm/ADT/StringRef.h"

#include <map>
#include <optional>
#include <string>
#include <vector>

namespace mlir {
namespace bmodelica {
class ModelOp;
} // namespace bmodelica

namespace cc {
namespace ode {

/// Product of symbols raised to positive powers. Symbols are variable names;
/// derivatives are named "der(<state>)".
using Monomial = std::map<std::string, unsigned>;

/// Sum of monomials with numeric coefficients. Terms with a zero coefficient
/// are never stored.
class Polynomial {
public:
  static Polynomial constant(double value);
  static Polynomial symbol(llvm::StringRef name);

  Polynomial operator+(const Polynomial &other) const;
  Polynomial operator-(const Polynomial &other) const;
  Polynomial operator*(const Polynomial &other) const;
  Polynomial operator-() const;
  Polynomial scaled(double factor) const;

  /// The value of a polynomial without symbols.
  std::optional<double> getConstantValue() const;

  const std::map<Monomial, double> &getTerms() const { return terms; }

  std::string str() const;

private:
  void addTerm(const Monomial &monomial, double coefficient);

  std::map<Monomial, double> terms;
};

/// A monomial of a state's right-hand side, split into its state part and its
/// rate-constant part.
struct ODETerm {
  /// Exponent of each state, indexed like ODESystem::states.
  std::vector<unsigned> stateExponents;
  /// Numeric coefficient (sign included).
  double coefficient;
  /// Parameters (and their exponents) multiplying the term.
  Monomial parameters;
};

/// An explicit polynomial ODE system extracted from a model.
struct ODESystem {
  /// States, in the order their variables are declared in the model.
  std::vector<std::string> states;
  /// Right-hand side of each state's derivative.
  std::vector<std::vector<ODETerm>> equations;
  /// Values of the parameters with a constant binding or start value.
  std::map<std::string, double> parameterValues;

  /// Human-readable form, e.g. "der(x) = -k*x".
  std::string str() const;
};

/// Format a number compactly, as printf "%g" does ("2", "0.5", "1e-06").
std::string formatNumber(double value);

/// Extract the ODE system of a model, or emit an error on the model and fail
/// if its equations are not an explicit scalar polynomial ODE system.
FailureOr<ODESystem> extractODESystem(bmodelica::ModelOp model);

} // namespace ode
} // namespace cc
} // namespace mlir

#endif // CHEMCOMP_LIB_CONVERSION_BASEMODELICATOACR_ODEPOLYNOMIAL_H
