//===- ODESystem.cpp - In-memory view of an ode.system --------------------===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#include "ChemComp/Dialect/ODE/ODESystem.h"
#include "ChemComp/Dialect/ODE/ODEOps.h"

#include "mlir/IR/Builders.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <tuple>

using namespace mlir;
using namespace mlir::cc::ode;

//===----------------------------------------------------------------------===//
// Numbers
//===----------------------------------------------------------------------===//

std::string mlir::cc::ode::formatNumber(double value) {
  char buffer[32];
  std::snprintf(buffer, sizeof(buffer), "%g", value);
  return buffer;
}

std::string mlir::cc::ode::formatExact(double value) {
  // The shortest "%.Ng" spelling that reads back exactly ("10" rather than
  // "1e+01").
  std::string best;
  char buffer[32];
  for (int precision = 1; precision <= 17; ++precision) {
    std::snprintf(buffer, sizeof(buffer), "%.*g", precision, value);
    if (std::strtod(buffer, nullptr) != value)
      continue;
    if (best.empty() || std::strlen(buffer) < best.size())
      best = buffer;
  }
  return best;
}

void mlir::cc::ode::printNumber(OpAsmPrinter &printer, double value) {
  // MLIR float literals need a '.': "2" -> "2.0", "1e-06" -> "1.0e-06".
  std::string text = formatExact(value);
  if (text.find_first_of(".ni") == std::string::npos) {
    size_t exponent = text.find('e');
    text.insert(exponent == std::string::npos ? text.size() : exponent, ".0");
  }
  printer << text;
}

ParseResult mlir::cc::ode::parseNumber(OpAsmParser &parser, double &value) {
  return parser.parseFloat(value);
}

//===----------------------------------------------------------------------===//
// Terms
//===----------------------------------------------------------------------===//

std::string mlir::cc::ode::factorsString(const Factors &factors) {
  std::string result;
  for (const auto &[name, exponent] : factors) {
    if (!result.empty())
      result += "*";
    result += name;
    if (exponent != 1)
      result += "^" + std::to_string(exponent);
  }
  return result;
}

std::string Term::str() const {
  std::string product = factorsString(params);
  std::string monomial = factorsString(states);
  if (!product.empty() && !monomial.empty())
    product += "*";
  product += monomial;

  double magnitude = std::fabs(coeff);
  std::string result = coeff < 0 ? "-" : "";
  if (product.empty())
    return result + formatNumber(magnitude);
  if (magnitude != 1.0)
    result += formatNumber(magnitude) + "*";
  return result + product;
}

//===----------------------------------------------------------------------===//
// System
//===----------------------------------------------------------------------===//

static Factors readFactors(ArrayAttr list) {
  Factors factors;
  for (Attribute attr : list)
    ++factors[cast<FlatSymbolRefAttr>(attr).getValue().str()];
  return factors;
}

static ArrayAttr writeFactors(Builder &builder, const Factors &factors) {
  SmallVector<Attribute> list;
  for (const auto &[name, exponent] : factors)
    for (unsigned i = 0; i < exponent; ++i)
      list.push_back(FlatSymbolRefAttr::get(builder.getContext(), name));
  return builder.getArrayAttr(list);
}

System System::read(SystemOp op) {
  System system;
  system.name = op.getSymName().str();
  for (Operation &nested : op.getBody().front()) {
    if (auto state = dyn_cast<StateOp>(nested)) {
      State s;
      s.name = state.getSymName().str();
      if (FloatAttr initial = state.getInitialAttr())
        s.initial = initial.getValueAsDouble();
      if (std::optional<StringRef> railOf = state.getRailOf())
        s.railOf = railOf->str();
      if (std::optional<uint64_t> sign = state.getRailSign())
        s.railSign = static_cast<int64_t>(*sign);
      system.states.push_back(std::move(s));
    } else if (auto param = dyn_cast<ParamOp>(nested)) {
      Param p{param.getSymName().str(), std::nullopt};
      if (FloatAttr value = param.getValueAttr())
        p.value = value.getValueAsDouble();
      system.params.push_back(std::move(p));
    } else if (auto term = dyn_cast<TermOp>(nested)) {
      system.terms.push_back({term.getState().str(),
                              term.getCoeff().convertToDouble(),
                              readFactors(term.getParams()),
                              readFactors(term.getStates()), term.getLoc()});
    }
  }
  return system;
}

void System::write(SystemOp op) const {
  Block &body = op.getBody().front();
  while (!body.empty())
    body.back().erase();

  OpBuilder builder = OpBuilder::atBlockEnd(&body);
  Location loc = op.getLoc();
  for (const State &state : states) {
    auto stateOp = StateOp::create(builder, loc, state.name,
                                   /*initial=*/nullptr, /*rail_of=*/nullptr,
                                   /*rail_sign=*/nullptr);
    if (state.initial)
      stateOp.setInitialAttr(builder.getF64FloatAttr(*state.initial));
    if (state.railOf)
      stateOp.setRailOfAttr(builder.getStringAttr(*state.railOf));
    if (state.railSign)
      stateOp.setRailSignAttr(builder.getI64IntegerAttr(*state.railSign));
  }
  for (const Param &param : params) {
    auto paramOp = ParamOp::create(builder, loc, param.name, nullptr);
    if (param.value)
      paramOp.setValueAttr(builder.getF64FloatAttr(*param.value));
  }
  for (const Term &term : terms) {
    TermOp::create(builder, term.loc.value_or(loc),
                   FlatSymbolRefAttr::get(builder.getContext(), term.state),
                   builder.getF64FloatAttr(term.coeff),
                   writeFactors(builder, term.params),
                   writeFactors(builder, term.states));
  }
}

int System::stateIndex(StringRef name) const {
  for (size_t i = 0; i < states.size(); ++i)
    if (states[i].name == name)
      return static_cast<int>(i);
  return -1;
}

std::optional<double> System::paramValue(StringRef name) const {
  for (const Param &param : params)
    if (param.name == name)
      return param.value;
  return std::nullopt;
}

void System::normalize() {
  using Key = std::tuple<int, Factors, Factors>;
  struct Sum {
    double value = 0;
    // Sum of the magnitudes, to tell cancellation from a small coefficient.
    double scale = 0;
    std::optional<Location> loc;
  };
  std::map<Key, Sum> sums;
  for (const Term &term : terms) {
    Sum &sum = sums[Key{stateIndex(term.state), term.states, term.params}];
    sum.value += term.coeff;
    sum.scale += std::fabs(term.coeff);
    if (!sum.loc)
      sum.loc = term.loc;
  }

  std::vector<Term> result;
  for (const auto &[key, sum] : sums) {
    if (std::fabs(sum.value) <= 1e-12 * sum.scale)
      continue;
    const auto &[index, monomial, product] = key;
    result.push_back({states[index].name, sum.value, product, monomial,
                      sum.loc});
  }
  terms = std::move(result);
}

std::string System::str() const {
  std::string result;
  for (const State &state : states) {
    std::string rhs;
    for (const Term &term : terms) {
      if (term.state != state.name)
        continue;
      std::string text = term.str();
      if (rhs.empty())
        rhs = text;
      else if (text.front() == '-')
        rhs += " - " + text.substr(1);
      else
        rhs += " + " + text;
    }
    if (!result.empty())
      result += "; ";
    result += "der(" + state.name + ") = " + (rhs.empty() ? "0" : rhs);
  }
  return result;
}

//===----------------------------------------------------------------------===//
// Hungarian Lemma
//===----------------------------------------------------------------------===//

std::vector<CrossEffect> mlir::cc::ode::findCrossEffects(const System &system) {
  // Parameters are non-negative (ode.param verifies it), so the sign of a
  // term is the sign of its numeric coefficient.
  std::vector<CrossEffect> effects;
  for (size_t i = 0; i < system.terms.size(); ++i) {
    const Term &term = system.terms[i];
    if (term.coeff < 0 && !term.states.count(term.state))
      effects.push_back({i});
  }
  return effects;
}

std::string mlir::cc::ode::describe(const System &system,
                                    const CrossEffect &effect) {
  const Term &term = system.terms[effect.termIndex];
  return "der(" + term.state + ") has " + term.str() +
         ", which does not contain " + term.state;
}
