//===- TranslateFromJSON.cpp - Import CRN JSON to MLIR -----*- C++ -*-===//
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
#include "mlir/IR/Builders.h"
#include "mlir/IR/BuiltinOps.h"
#include "mlir/IR/BuiltinTypes.h"
#include "mlir/IR/Location.h"
#include "llvm/Support/JSON.h"
#include "llvm/Support/SourceMgr.h"

using namespace mlir;

OwningOpRef<Operation *>
mlir::cc::acr::translateFromJSON(llvm::SourceMgr &sourceMgr, MLIRContext *context) {
  auto loc = FileLineColLoc::get(
      StringAttr::get(context,
                      sourceMgr.getMemoryBuffer(1)->getBufferIdentifier()),
      0, 0);

  // Parse the JSON input.
  StringRef jsonStr = sourceMgr.getMemoryBuffer(1)->getBuffer();
  auto parsed = llvm::json::parse(jsonStr);
  if (!parsed) {
    emitError(loc) << "failed to parse JSON: "
                   << llvm::toString(parsed.takeError());
    return {};
  }

  auto *root = parsed->getAsObject();
  if (!root) {
    emitError(loc) << "expected JSON object at top level";
    return {};
  }

  // Extract reactions array.
  auto *reactions = root->getArray("reactions");
  if (!reactions) {
    emitError(loc) << "expected 'reactions' array in JSON";
    return {};
  }

  int64_t R = reactions->size();
  if (R == 0) {
    emitError(loc) << "reactions array is empty";
    return {};
  }

  // Determine number of species from first reaction.
  auto *firstReaction = (*reactions)[0].getAsObject();
  if (!firstReaction) {
    emitError(loc) << "expected reaction object";
    return {};
  }
  auto *speciesArr = firstReaction->getArray("species");
  if (!speciesArr) {
    emitError(loc) << "expected 'species' array in reaction";
    return {};
  }
  int64_t S = speciesArr->size();

  // Build stoichiometric matrix and extract rates.
  SmallVector<int32_t> matrixValues;
  SmallVector<double> rates;

  for (int64_t r = 0; r < R; ++r) {
    auto *rxn = (*reactions)[r].getAsObject();
    if (!rxn) {
      emitError(loc) << "expected reaction object at index " << r;
      return {};
    }

    auto *reactants = rxn->getArray("reactants");
    auto *products = rxn->getArray("products");
    if (!reactants || !products) {
      emitError(loc) << "reaction " << r
                     << " missing 'reactants' or 'products'";
      return {};
    }
    if ((int64_t)reactants->size() != S || (int64_t)products->size() != S) {
      emitError(loc) << "reaction " << r << " has mismatched species count";
      return {};
    }

    for (int64_t s = 0; s < S; ++s) {
      auto reactant = (*reactants)[s].getAsInteger();
      auto product = (*products)[s].getAsInteger();
      if (!reactant || !product) {
        emitError(loc) << "non-integer value in reaction " << r;
        return {};
      }
      // stoich[r][s] = products[s] - reactants[s]
      matrixValues.push_back(static_cast<int32_t>(*product - *reactant));
    }

    // Extract rate.
    auto rate = rxn->getNumber("rate");
    if (!rate) {
      emitError(loc) << "reaction " << r << " missing 'rate'";
      return {};
    }
    rates.push_back(*rate);
  }

  // Build the MLIR module.
  context->loadDialect<cc::acr::ACRDialect, arith::ArithDialect>();
  OpBuilder builder(context);

  auto module = ModuleOp::create(loc);
  builder.setInsertionPointToEnd(module.getBody());

  // Types.
  auto speciesType = cc::acr::SpeciesType::get(context);
  auto speciesTensorType = RankedTensorType::get({S}, speciesType);
  auto rateTensorType = RankedTensorType::get({R}, builder.getF32Type());
  auto funcType = builder.getFunctionType({speciesTensorType, rateTensorType},
                                          {speciesTensorType});

  // Create acr.crn op.
  auto crnOp = cc::acr::CRNOp::create(builder, loc, "crn", funcType,
                                      /*sym_visibility=*/StringAttr{},
                                      /*arg_attrs=*/ArrayAttr{},
                                      /*res_attrs=*/ArrayAttr{});

  // Create the body block with arguments.
  Block *body = crnOp.addEntryBlock();
  builder.setInsertionPointToStart(body);

  Value speciesArg = body->getArgument(0);
  Value rateArg = body->getArgument(1);

  // Build stoich_matrix attribute.
  auto matrixTensorType = RankedTensorType::get({R, S}, builder.getI32Type());
  auto stoichAttr = DenseIntElementsAttr::get(matrixTensorType, matrixValues);

  // Create acr.crn_generic op (tensor form).
  auto genericOp =
      cc::acr::CRNGenericOp::create(builder, loc, TypeRange{speciesTensorType},
                                    ValueRange{speciesArg, rateArg});
  genericOp->setAttr("stoich_matrix", stoichAttr);

  // Create acr.return.
  cc::acr::ReturnOp::create(builder, loc, genericOp->getResults());

  return OwningOpRef<Operation *>(module.getOperation());
}
