//===- ODEPolynomial.cpp - Polynomial ODE systems from Base Modelica ------===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#include "ODEPolynomial.h"

#include "marco/Dialect/BaseModelica/IR/Dialect.h"

#include "mlir/IR/Diagnostics.h"
#include "llvm/ADT/DenseMap.h"
#include "llvm/ADT/StringMap.h"
#include "llvm/Support/Format.h"
#include "llvm/Support/raw_ostream.h"

#include <cmath>

using namespace mlir;
using namespace mlir::cc::ode;
namespace bmodelica = mlir::bmodelica;

//===----------------------------------------------------------------------===//
// Polynomial
//===----------------------------------------------------------------------===//

Polynomial Polynomial::constant(double value) {
  Polynomial result;
  result.addTerm({}, value);
  return result;
}

Polynomial Polynomial::symbol(StringRef name) {
  Polynomial result;
  result.addTerm({{name.str(), 1}}, 1.0);
  return result;
}

void Polynomial::addTerm(const Monomial &monomial, double coefficient) {
  double &value = terms[monomial];
  value += coefficient;
  if (value == 0.0)
    terms.erase(monomial);
}

Polynomial Polynomial::operator+(const Polynomial &other) const {
  Polynomial result = *this;
  for (const auto &[monomial, coefficient] : other.terms)
    result.addTerm(monomial, coefficient);
  return result;
}

Polynomial Polynomial::operator-(const Polynomial &other) const {
  return *this + (-other);
}

Polynomial Polynomial::operator-() const { return scaled(-1.0); }

Polynomial Polynomial::scaled(double factor) const {
  Polynomial result;
  for (const auto &[monomial, coefficient] : terms)
    result.addTerm(monomial, coefficient * factor);
  return result;
}

Polynomial Polynomial::operator*(const Polynomial &other) const {
  Polynomial result;
  for (const auto &[lhsMonomial, lhsCoefficient] : terms) {
    for (const auto &[rhsMonomial, rhsCoefficient] : other.terms) {
      Monomial product = lhsMonomial;
      for (const auto &[symbol, exponent] : rhsMonomial)
        product[symbol] += exponent;
      result.addTerm(product, lhsCoefficient * rhsCoefficient);
    }
  }
  return result;
}

std::optional<double> Polynomial::getConstantValue() const {
  if (terms.empty())
    return 0.0;
  if (terms.size() == 1 && terms.begin()->first.empty())
    return terms.begin()->second;
  return std::nullopt;
}

std::string mlir::cc::ode::formatNumber(double value) {
  std::string result;
  llvm::raw_string_ostream os(result);
  os << llvm::format("%g", value);
  return result;
}

static std::string formatMonomial(const Monomial &monomial) {
  std::string result;
  for (const auto &[symbol, exponent] : monomial) {
    if (!result.empty())
      result += "*";
    result += symbol;
    if (exponent != 1)
      result += "^" + std::to_string(exponent);
  }
  return result;
}

/// Format `coefficient * monomial` as a signed term, e.g. "- 2*k*x".
static std::string formatTerm(double coefficient, const Monomial &monomial,
                              bool first) {
  std::string sign =
      coefficient < 0 ? (first ? "-" : " - ") : (first ? "" : " + ");
  double magnitude = std::fabs(coefficient);
  if (monomial.empty())
    return sign + formatNumber(magnitude);
  if (magnitude == 1.0)
    return sign + formatMonomial(monomial);
  return sign + formatNumber(magnitude) + "*" + formatMonomial(monomial);
}

std::string Polynomial::str() const {
  if (terms.empty())
    return "0";
  std::string result;
  bool first = true;
  for (const auto &[monomial, coefficient] : terms) {
    result += formatTerm(coefficient, monomial, first);
    first = false;
  }
  return result;
}

std::string ODESystem::str() const {
  std::string result;
  for (size_t i = 0; i < states.size(); ++i) {
    Polynomial rhs;
    for (const ODETerm &term : equations[i]) {
      Polynomial product = Polynomial::constant(term.coefficient);
      auto multiply = [&](const std::string &symbol, unsigned exponent) {
        for (unsigned k = 0; k < exponent; ++k)
          product = product * Polynomial::symbol(symbol);
      };
      for (const auto &[parameter, exponent] : term.parameters)
        multiply(parameter, exponent);
      for (size_t j = 0; j < states.size(); ++j)
        multiply(states[j], term.stateExponents[j]);
      rhs = rhs + product;
    }
    if (i != 0)
      result += "; ";
    result += "der(" + states[i] + ") = " + rhs.str();
  }
  return result;
}

//===----------------------------------------------------------------------===//
// Expression folding
//===----------------------------------------------------------------------===//

namespace {

/// Folds bmodelica expressions into polynomials.
class ExpressionFolder {
public:
  /// `derivativeOf` maps derivative variables (der_x) to their state (x).
  explicit ExpressionFolder(const llvm::StringMap<std::string> &derivativeOf)
      : derivativeOf(derivativeOf) {}

  FailureOr<Polynomial> fold(Value value);

private:
  FailureOr<Polynomial> foldOp(Operation *op);
  FailureOr<Polynomial> unsupported(Operation *op, const Twine &reason);

  const llvm::StringMap<std::string> &derivativeOf;
  llvm::DenseMap<Value, Polynomial> cache;
};

} // namespace

static std::string derivativeSymbol(StringRef state) {
  return ("der(" + state + ")").str();
}

static std::optional<double> getNumericValue(Attribute attr) {
  if (auto real = dyn_cast<bmodelica::RealAttr>(attr))
    return real.getValue().convertToDouble();
  if (auto integer = dyn_cast<bmodelica::IntegerAttr>(attr))
    return static_cast<double>(integer.getValue().getSExtValue());
  if (auto real = dyn_cast<FloatAttr>(attr))
    return real.getValueAsDouble();
  if (auto integer = dyn_cast<IntegerAttr>(attr))
    return static_cast<double>(integer.getValue().getSExtValue());
  return std::nullopt;
}

FailureOr<Polynomial> ExpressionFolder::unsupported(Operation *op,
                                                    const Twine &reason) {
  op->emitError() << reason;
  return failure();
}

FailureOr<Polynomial> ExpressionFolder::fold(Value value) {
  auto it = cache.find(value);
  if (it != cache.end())
    return it->second;

  Operation *op = value.getDefiningOp();
  if (!op) {
    emitError(value.getLoc())
        << "expected the ODE to be built from operations, found a block "
           "argument (equations with inductions are not supported)";
    return failure();
  }

  FailureOr<Polynomial> result = foldOp(op);
  if (succeeded(result))
    cache[value] = *result;
  return result;
}

FailureOr<Polynomial> ExpressionFolder::foldOp(Operation *op) {
  if (auto constant = dyn_cast<bmodelica::ConstantOp>(op)) {
    if (std::optional<double> value = getNumericValue(constant.getValue()))
      return Polynomial::constant(*value);
    return unsupported(op, "unsupported constant in an ODE expression");
  }

  if (auto get = dyn_cast<bmodelica::VariableGetOp>(op)) {
    StringRef name = get.getVariable();
    auto it = derivativeOf.find(name);
    if (it != derivativeOf.end())
      return Polynomial::symbol(derivativeSymbol(it->second));
    return Polynomial::symbol(name);
  }

  if (auto der = dyn_cast<bmodelica::DerOp>(op)) {
    auto get = der.getOperand().getDefiningOp<bmodelica::VariableGetOp>();
    if (!get)
      return unsupported(op, "der() is only supported on a variable");
    return Polynomial::symbol(derivativeSymbol(get.getVariable()));
  }

  if (auto neg = dyn_cast<bmodelica::NegateOp>(op)) {
    FailureOr<Polynomial> operand = fold(neg.getOperand());
    if (failed(operand))
      return failure();
    return -*operand;
  }

  // Binary arithmetic.
  auto foldBinary =
      [&](Value lhs,
          Value rhs) -> FailureOr<std::pair<Polynomial, Polynomial>> {
    FailureOr<Polynomial> lhsPoly = fold(lhs);
    if (failed(lhsPoly))
      return failure();
    FailureOr<Polynomial> rhsPoly = fold(rhs);
    if (failed(rhsPoly))
      return failure();
    return std::make_pair(*lhsPoly, *rhsPoly);
  };

  if (auto add = dyn_cast<bmodelica::AddOp>(op)) {
    auto operands = foldBinary(add.getLhs(), add.getRhs());
    if (failed(operands))
      return failure();
    return operands->first + operands->second;
  }

  if (auto sub = dyn_cast<bmodelica::SubOp>(op)) {
    auto operands = foldBinary(sub.getLhs(), sub.getRhs());
    if (failed(operands))
      return failure();
    return operands->first - operands->second;
  }

  if (auto mul = dyn_cast<bmodelica::MulOp>(op)) {
    auto operands = foldBinary(mul.getLhs(), mul.getRhs());
    if (failed(operands))
      return failure();
    return operands->first * operands->second;
  }

  if (auto div = dyn_cast<bmodelica::DivOp>(op)) {
    auto operands = foldBinary(div.getLhs(), div.getRhs());
    if (failed(operands))
      return failure();
    std::optional<double> divisor = operands->second.getConstantValue();
    if (!divisor || *divisor == 0.0)
      return unsupported(op, "only division by a non-zero numeric constant "
                             "is supported in an ODE expression");
    return operands->first.scaled(1.0 / *divisor);
  }

  if (auto pow = dyn_cast<bmodelica::PowOp>(op)) {
    auto operands = foldBinary(pow.getBase(), pow.getExponent());
    if (failed(operands))
      return failure();
    std::optional<double> exponent = operands->second.getConstantValue();
    if (!exponent || *exponent < 0 || std::floor(*exponent) != *exponent ||
        *exponent > 16)
      return unsupported(op, "only small non-negative integer exponents are "
                             "supported in an ODE expression");
    Polynomial result = Polynomial::constant(1.0);
    for (int i = 0; i < static_cast<int>(*exponent); ++i)
      result = result * operands->first;
    return result;
  }

  return unsupported(op, "unsupported operation '" +
                             op->getName().getStringRef() +
                             "' in an ODE expression; only constants, "
                             "variables, der, +, -, *, / by a constant and "
                             "integer powers are supported");
}

//===----------------------------------------------------------------------===//
// Model extraction
//===----------------------------------------------------------------------===//

/// Value of a parameter given by a constant binding equation or start value.
static std::optional<double> getRegionConstant(Region &region) {
  if (region.empty())
    return std::nullopt;
  Operation *terminator = region.front().getTerminator();
  if (!terminator || terminator->getNumOperands() != 1)
    return std::nullopt;
  llvm::StringMap<std::string> noDerivatives;
  ExpressionFolder folder(noDerivatives);
  // Folding may fail on non-constant expressions: diagnose silently.
  ScopedDiagnosticHandler silence(region.getContext(),
                                  [](Diagnostic &) { return success(); });
  FailureOr<Polynomial> value = folder.fold(terminator->getOperand(0));
  if (failed(value))
    return std::nullopt;
  return value->getConstantValue();
}

FailureOr<ODESystem> mlir::cc::ode::extractODESystem(bmodelica::ModelOp model) {
  ODESystem system;

  // Classify the variables.
  llvm::StringMap<std::string> derivativeOf;
  llvm::StringMap<bool> isParameter;
  bmodelica::DerivativesMap derivativesMap = model.getDerivativesMap();
  SmallVector<bmodelica::VariableOp> variables;
  model.collectVariables(variables);

  for (bmodelica::VariableOp variable : variables) {
    if (!variable.getType().isScalar())
      return model.emitError()
             << "array variable '" << variable.getSymName()
             << "' is not supported; only scalar ODE models can be mapped "
                "to chemical reaction networks";
    if (variable.getType().isReadOnly())
      isParameter[variable.getSymName()] = true;
    auto symbol = FlatSymbolRefAttr::get(variable.getSymNameAttr());
    if (auto derived = derivativesMap.getDerivedVariable(symbol))
      derivativeOf[variable.getSymName()] =
          derived->getRootReference().getValue().str();
  }

  // Parameter values.
  for (auto binding : model.getOps<bmodelica::BindingEquationOp>())
    if (isParameter.count(binding.getVariable()))
      if (auto value = getRegionConstant(binding->getRegion(0)))
        system.parameterValues[binding.getVariable().str()] = *value;
  for (auto start : model.getOps<bmodelica::StartOp>()) {
    StringRef name = start.getVariable().getRootReference().getValue();
    if (isParameter.count(name) && !system.parameterValues.count(name.str()))
      if (auto value = getRegionConstant(start->getRegion(0)))
        system.parameterValues[name.str()] = *value;
  }

  // Fold each dynamic equation into `der(x) = rhs`.
  SmallVector<bmodelica::EquationInstanceOp> equations;
  model.collectDynamicEquations(equations);
  if (equations.empty())
    return model.emitError() << "model has no dynamic equations";

  llvm::StringMap<Polynomial> rhsOf;
  ExpressionFolder folder(derivativeOf);

  for (bmodelica::EquationInstanceOp instance : equations) {
    bmodelica::EquationTemplateOp eqTemplate = instance.getTemplate();
    if (!eqTemplate.getInductionVariables().empty())
      return instance.emitError()
             << "equations with induction variables are not supported";

    auto sides = dyn_cast<bmodelica::EquationSidesOp>(
        eqTemplate->getRegion(0).front().getTerminator());
    if (!sides || sides.getLhsValues().size() != 1 ||
        sides.getRhsValues().size() != 1)
      return instance.emitError() << "expected a scalar equation";

    FailureOr<Polynomial> lhs = folder.fold(sides.getLhsValues()[0]);
    if (failed(lhs))
      return failure();
    FailureOr<Polynomial> rhs = folder.fold(sides.getRhsValues()[0]);
    if (failed(rhs))
      return failure();
    Polynomial residual = *lhs - *rhs;

    // The residual must be `c * der(x) + rest`, with rest free of derivatives.
    std::optional<std::string> state;
    double derivativeCoefficient = 0.0;
    Polynomial rest;
    for (const auto &[monomial, coefficient] : residual.getTerms()) {
      bool hasDerivative = llvm::any_of(monomial, [](const auto &entry) {
        return StringRef(entry.first).starts_with("der(");
      });
      if (!hasDerivative) {
        Polynomial term = Polynomial::constant(coefficient);
        for (const auto &[symbol, exponent] : monomial)
          for (unsigned i = 0; i < exponent; ++i)
            term = term * Polynomial::symbol(symbol);
        rest = rest + term;
        continue;
      }
      if (monomial.size() != 1 || monomial.begin()->second != 1 || state)
        return instance.emitError()
               << "expected an explicit ODE equation 'der(x) = f(...)' that "
                  "is linear in a single derivative";
      StringRef symbol = monomial.begin()->first;
      state = symbol.drop_front(4).drop_back(1).str();
      derivativeCoefficient = coefficient;
    }
    if (!state)
      return instance.emitError()
             << "algebraic equations are not supported; every equation must "
                "define the derivative of a state";
    if (rhsOf.count(*state))
      return instance.emitError()
             << "state '" << *state << "' has more than one ODE equation";

    rhsOf[*state] = rest.scaled(-1.0 / derivativeCoefficient);
  }

  // States in declaration order.
  for (bmodelica::VariableOp variable : variables)
    if (rhsOf.count(variable.getSymName()))
      system.states.push_back(variable.getSymName().str());

  llvm::StringMap<size_t> stateIndex;
  for (size_t i = 0; i < system.states.size(); ++i)
    stateIndex[system.states[i]] = i;

  // Split each monomial into its state and rate-constant parts.
  for (const std::string &state : system.states) {
    std::vector<ODETerm> terms;
    for (const auto &[monomial, coefficient] : rhsOf[state].getTerms()) {
      ODETerm term;
      term.stateExponents.assign(system.states.size(), 0);
      term.coefficient = coefficient;
      for (const auto &[symbol, exponent] : monomial) {
        auto it = stateIndex.find(symbol);
        if (it != stateIndex.end()) {
          term.stateExponents[it->second] = exponent;
          continue;
        }
        if (isParameter.count(symbol)) {
          term.parameters[symbol] = exponent;
          continue;
        }
        return model.emitError()
               << "the ODE of '" << state << "' uses '" << symbol
               << "', which is neither a state with its own ODE nor a "
                  "parameter";
      }
      terms.push_back(std::move(term));
    }
    system.equations.push_back(std::move(terms));
  }

  return system;
}
