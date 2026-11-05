//===- TranslateODEToJSON.cpp - Export ode.system to JSON -----------------===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// Writes an ode.system as the reference for checking a compiled network:
//
//   {"name": "Decay",
//    "states": [{"name": "x", "initial": 1}],
//    "terms": [{"state": "x", "coeff": -0.5, "monomial": {"x": 1}}]}
//
// Parameters are folded into the coefficients, so every parameter used by a
// term needs a value.
//
//===----------------------------------------------------------------------===//

#include "ChemComp/Dialect/ODE/ODEOps.h"
#include "ChemComp/Dialect/ODE/ODESystem.h"
#include "ChemComp/Target/JSON/JSONEmitter.h"

#include "llvm/Support/Format.h"
#include "llvm/Support/JSON.h"

#include <cmath>

using namespace mlir;
using namespace mlir::cc::ode;

static void emitDouble(llvm::json::OStream &J, double value) {
  SmallString<32> buffer;
  llvm::raw_svector_ostream os(buffer);
  os << llvm::format("%.17g", value);
  J.rawValue(buffer);
}

LogicalResult mlir::cc::ode::translateODEToJSON(Operation *op,
                                                llvm::raw_ostream &os) {
  SmallVector<SystemOp> systems;
  op->walk([&](SystemOp system) { systems.push_back(system); });
  if (systems.size() != 1)
    return op->emitError("ode-to-json expects exactly one ode.system, found ")
           << systems.size();

  System system = System::read(systems.front());
  for (const Term &term : system.terms)
    for (const auto &[param, exponent] : term.params)
      if (!system.paramValue(param))
        return systems.front().emitError("ode-to-json: parameter '")
               << param << "' has no value";

  llvm::json::OStream J(os, /*IndentSize=*/2);
  J.object([&] {
    J.attribute("name", system.name);
    J.attributeArray("states", [&] {
      for (const State &state : system.states) {
        J.object([&] {
          J.attribute("name", state.name);
          J.attributeBegin("initial");
          emitDouble(J, state.initial.value_or(0.0));
          J.attributeEnd();
        });
      }
    });
    J.attributeArray("terms", [&] {
      for (const Term &term : system.terms) {
        double coeff = term.coeff;
        for (const auto &[param, exponent] : term.params)
          coeff *= std::pow(*system.paramValue(param), exponent);
        J.object([&] {
          J.attribute("state", term.state);
          J.attributeBegin("coeff");
          emitDouble(J, coeff);
          J.attributeEnd();
          J.attributeObject("monomial", [&] {
            for (const auto &[state, exponent] : term.states)
              J.attribute(state, static_cast<int64_t>(exponent));
          });
        });
      }
    });
  });
  os << "\n";
  return success();
}
