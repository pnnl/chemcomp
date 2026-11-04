//===- TranslateToJSON.cpp - Export ACR CRN to JSON --------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#include "ChemComp/Dialect/ACR/ACRDialect.h"
#include "ChemComp/Dialect/ACR/ACROps.h"
#include "ChemComp/Dialect/ACR/ACRTypes.h"
#include "ChemComp/Target/JSON/JSONEmitter.h"
#include "mlir/Dialect/Arith/IR/Arith.h"
#include "mlir/IR/BuiltinOps.h"
#include "mlir/IR/Operation.h"
#include "llvm/Support/Format.h"
#include "llvm/Support/JSON.h"
#include "llvm/Support/raw_ostream.h"

using namespace mlir;

/// Try to extract a constant float value from a Value by walking the def-use
/// chain looking for an arith.constant. Returns the APFloat to preserve
/// exact f32 semantics.
static std::optional<APFloat> extractConstantRate(Value val) {
  if (auto constOp = val.getDefiningOp<arith::ConstantOp>()) {
    if (auto floatAttr = dyn_cast<FloatAttr>(constOp.getValue()))
      return floatAttr.getValue();
  }
  return std::nullopt;
}

/// Emit a floating point value to the JSON stream. Uses shortest
/// representation that round-trips through f32.
static void emitRate(llvm::json::OStream &J, const APFloat &val) {
  // Convert to double and format with enough precision for f32 round-trip
  // (6 significant digits suffice for float32).
  double d = val.convertToDouble();
  SmallString<32> buf;
  llvm::raw_svector_ostream os(buf);
  // Use %g-style formatting with 6 significant digits — this gives
  // clean output like "0.1" for the f32 nearest to 0.1.
  os << llvm::format("%g", d);
  J.rawValue(buf);
}

/// Emit a double compactly ("2", "0.5", "1.1") with up to 15 significant
/// digits, so values recorded as f64 reach the emulator unchanged.
static void emitDouble(llvm::json::OStream &J, double d) {
  SmallString<32> buf;
  llvm::raw_svector_ostream os(buf);
  os << llvm::format("%.15g", d);
  J.rawValue(buf);
}

/// Emit a single CRNGenericOp as JSON.
static LogicalResult emitCRNGenericOp(cc::acr::CRNGenericOp op,
                                      llvm::json::OStream &J) {
  auto stoichAttr = op->getAttrOfType<DenseIntElementsAttr>("stoich_matrix");
  if (!stoichAttr)
    return op.emitError("missing stoich_matrix attribute");

  auto stoichType = stoichAttr.getType();
  int64_t R = stoichType.getDimSize(0); // reactions
  int64_t S = stoichType.getDimSize(1); // species

  // Read stoich matrix values into a flat vector, then index as [r*S + s].
  SmallVector<int32_t> matrix;
  for (auto val : stoichAttr.getValues<int32_t>())
    matrix.push_back(val);

  // Reactants and products. The net matrix loses catalysts (X -> 2X has net
  // +1, which reads back as 0 -> X), so prefer the full matrices when the
  // producer recorded them; otherwise split the net matrix by sign.
  SmallVector<int32_t> reactants(R * S), products(R * S);
  auto reactantAttr = op->getAttrOfType<DenseIntElementsAttr>("reactant_matrix");
  auto productAttr = op->getAttrOfType<DenseIntElementsAttr>("product_matrix");
  if (reactantAttr && productAttr &&
      reactantAttr.getType().getShape() == stoichType.getShape() &&
      productAttr.getType().getShape() == stoichType.getShape()) {
    llvm::copy(reactantAttr.getValues<int32_t>(), reactants.begin());
    llvm::copy(productAttr.getValues<int32_t>(), products.begin());
  } else {
    for (int64_t i = 0; i < R * S; ++i) {
      reactants[i] = matrix[i] < 0 ? -matrix[i] : 0;
      products[i] = matrix[i] > 0 ? matrix[i] : 0;
    }
  }

  // Prefer the concrete names attached by the crn-db-match pass; fall back to
  // generic S0, S1, ... where the CRN was never matched or a species was left
  // unbound. Without this the matched names would stop at the IR and never
  // reach the emulator.
  auto namesAttr = op->getAttrOfType<ArrayAttr>("species_names");
  SmallVector<std::string> speciesNames;
  for (int64_t s = 0; s < S; ++s) {
    if (namesAttr && s < (int64_t)namesAttr.size()) {
      if (auto name = dyn_cast<StringAttr>(namesAttr[s])) {
        if (!name.getValue().empty()) {
          speciesNames.push_back(name.getValue().str());
          continue;
        }
      }
    }
    speciesNames.push_back(("S" + Twine(s)).str());
  }

  // Extract rate constants.
  // Detect tensor vs scalar form.
  auto inputs = op.getInputs();
  APFloat zeroF32(0.0f);
  SmallVector<APFloat> rates(R, zeroF32);
  // Rates from `rate_values` keep their double precision.
  std::optional<SmallVector<double>> exactRates;

  bool isTensorForm =
      !inputs.empty() && isa<RankedTensorType>(inputs[0].getType());

  // Rates recorded on the op (e.g. by convert-bmodelica-to-acr) take the place
  // of a constant rate tensor.
  auto rateValues = op->getAttrOfType<DenseF64ArrayAttr>("rate_values");

  if (isTensorForm) {
    // Tensor form: inputs[1] is tensor<Rxf32>.
    // Try to trace it to an arith.constant.
    if (inputs.size() >= 2) {
      if (rateValues && static_cast<int64_t>(rateValues.size()) == R) {
        exactRates.emplace(rateValues.asArrayRef());
      } else if (auto constOp = inputs[1].getDefiningOp<arith::ConstantOp>()) {
        if (auto denseAttr =
                dyn_cast<DenseFPElementsAttr>(constOp.getValue())) {
          unsigned idx = 0;
          for (auto val : denseAttr.getValues<APFloat>()) {
            if (idx < (unsigned)R)
              rates[idx] = val;
            ++idx;
          }
        }
      } else {
        op.emitWarning("could not resolve rate tensor to constant; "
                       "using 0.0 for all rates");
      }
    }
  } else {
    // Scalar form: inputs[S..S+R-1] are f32 rate values.
    for (int64_t r = 0; r < R; ++r) {
      int64_t idx = S + r;
      if (idx < (int64_t)inputs.size()) {
        if (auto rate = extractConstantRate(inputs[idx]))
          rates[r] = *rate;
        else
          op.emitWarning("could not resolve rate constant for reaction ")
              << r << "; using 0.0";
      }
    }
  }

  // Emit JSON.
  J.object([&] {
    J.attributeArray("reactions", [&] {
      for (int64_t r = 0; r < R; ++r) {
        J.object([&] {
          J.attributeArray("species", [&] {
            for (int64_t s = 0; s < S; ++s)
              J.value(speciesNames[s]);
          });
          J.attributeArray("reactants", [&] {
            for (int64_t s = 0; s < S; ++s)
              J.value(reactants[r * S + s]);
          });
          J.attributeArray("products", [&] {
            for (int64_t s = 0; s < S; ++s)
              J.value(products[r * S + s]);
          });
          J.attributeBegin("rate");
          if (exactRates)
            emitDouble(J, (*exactRates)[r]);
          else
            emitRate(J, rates[r]);
          J.attributeEnd();
          J.attribute("name", ("R" + Twine(r + 1)).str());
        });
      }
    });
    // Initial amounts, when the producer recorded them (e.g. ODE start
    // values); 1 otherwise.
    auto initialValues = op->getAttrOfType<DenseF64ArrayAttr>("initial_values");
    bool haveInitial =
        initialValues && static_cast<int64_t>(initialValues.size()) == S;
    J.attributeObject("populations", [&] {
      for (int64_t s = 0; s < S; ++s) {
        if (!haveInitial) {
          J.attribute(speciesNames[s], 1);
          continue;
        }
        J.attributeBegin(speciesNames[s]);
        emitDouble(J, initialValues[s]);
        J.attributeEnd();
      }
    });
    J.attributeArray("tracked_species", [&] {});

    // Which ODE state each species stands for, and with which sign: a
    // dual-rail pair x = x+ - x- maps both rails to x, with signs +1 and -1.
    auto stateNames = op->getAttrOfType<ArrayAttr>("ode_state_names");
    auto stateSigns = op->getAttrOfType<DenseI64ArrayAttr>("ode_state_signs");
    if (stateNames && static_cast<int64_t>(stateNames.size()) == S) {
      J.attributeObject("species_map", [&] {
        for (int64_t s = 0; s < S; ++s) {
          auto state = dyn_cast<StringAttr>(stateNames[s]);
          if (!state)
            continue;
          int64_t sign = 1;
          if (stateSigns && static_cast<int64_t>(stateSigns.size()) == S)
            sign = stateSigns[s];
          J.attributeObject(speciesNames[s], [&] {
            J.attribute("state", state.getValue());
            J.attribute("sign", sign);
          });
        }
      });
    }
  });

  return success();
}

LogicalResult mlir::cc::acr::translateToJSON(Operation *op, llvm::raw_ostream &os) {
  bool foundAny = false;
  LogicalResult result = success();

  op->walk([&](cc::acr::CRNGenericOp genericOp) {
    if (failed(result))
      return;
    if (foundAny)
      os << "\n";
    llvm::json::OStream J(os, /*IndentSize=*/4);
    if (failed(emitCRNGenericOp(genericOp, J)))
      result = failure();
    os << "\n";
    foundAny = true;
  });

  if (!foundAny)
    return op->emitError("no acr.crn_generic ops found");

  return result;
}
