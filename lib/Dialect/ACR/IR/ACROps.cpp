//===- ACROps.cpp - ACR dialect ops ---------------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#include "ChemComp/Dialect/ACR/ACROps.h"
#include "ChemComp/Dialect/ACR/ACRDialect.h"
#include "ChemComp/Dialect/ACR/ACRTypes.h"
#include "mlir/IR/Builders.h"
#include "mlir/IR/BuiltinTypes.h"
#include "mlir/Interfaces/FunctionImplementation.h"

using namespace mlir;
using namespace mlir::acr;

//===----------------------------------------------------------------------===//
// ACR_ReactionOp — custom parser/printer (function-like)
//===----------------------------------------------------------------------===//

ParseResult ReactionOp::parse(OpAsmParser &parser, OperationState &result) {
  auto buildFuncType =
      [](Builder &builder, ArrayRef<Type> argTypes, ArrayRef<Type> results,
         function_interface_impl::VariadicFlag,
         std::string &) { return builder.getFunctionType(argTypes, results); };

  return function_interface_impl::parseFunctionOp(
      parser, result, /*allowVariadic=*/false,
      getFunctionTypeAttrName(result.name), buildFuncType,
      getArgAttrsAttrName(result.name), getResAttrsAttrName(result.name));
}

void ReactionOp::print(OpAsmPrinter &p) {
  function_interface_impl::printFunctionOp(
      p, *this, /*isVariadic=*/false, getFunctionTypeAttrName(),
      getArgAttrsAttrName(), getResAttrsAttrName());
}

//===----------------------------------------------------------------------===//
// ACR_CRNOp — custom parser/printer (function-like)
//===----------------------------------------------------------------------===//

ParseResult CRNOp::parse(OpAsmParser &parser, OperationState &result) {
  auto buildFuncType =
      [](Builder &builder, ArrayRef<Type> argTypes, ArrayRef<Type> results,
         function_interface_impl::VariadicFlag,
         std::string &) { return builder.getFunctionType(argTypes, results); };

  return function_interface_impl::parseFunctionOp(
      parser, result, /*allowVariadic=*/false,
      getFunctionTypeAttrName(result.name), buildFuncType,
      getArgAttrsAttrName(result.name), getResAttrsAttrName(result.name));
}

void CRNOp::print(OpAsmPrinter &p) {
  function_interface_impl::printFunctionOp(
      p, *this, /*isVariadic=*/false, getFunctionTypeAttrName(),
      getArgAttrsAttrName(), getResAttrsAttrName());
}

//===----------------------------------------------------------------------===//
// ACR_ReactOp — verifier
//===----------------------------------------------------------------------===//

LogicalResult ReactOp::verify() {
  auto inputs = getInputs().getTypes();
  auto results = getResults().getTypes();

  if (inputs.empty())
    return emitOpError("expected at least one input");

  // Detect tensor vs scalar form by checking the first input.
  if (auto speciesTensor = dyn_cast<RankedTensorType>(inputs[0])) {
    // --- Tensor form ---
    // Expected: (tensor<Nx!acr.species>, tensor<Nxf32>, f32)
    //        -> tensor<Nx!acr.species>
    if (!isa<acr::SpeciesType>(speciesTensor.getElementType()))
      return emitOpError("tensor form requires first input to be "
                         "tensor<Nx!acr.species>");
    if (speciesTensor.getRank() != 1)
      return emitOpError("species tensor must be 1-D");

    int64_t N = speciesTensor.getDimSize(0);

    if (inputs.size() != 3)
      return emitOpError("tensor form expects exactly 3 inputs: "
                         "tensor<Nx!acr.species>, tensor<Nxf32>, f32");

    // Check stoichiometry tensor.
    auto stoichTensor = dyn_cast<RankedTensorType>(inputs[1]);
    if (!stoichTensor || stoichTensor.getRank() != 1 ||
        !stoichTensor.getElementType().isF32())
      return emitOpError("tensor form expects second input to be "
                         "tensor<Nxf32>");
    if (stoichTensor.getDimSize(0) != N)
      return emitOpError("stoichiometry tensor dimension (")
             << stoichTensor.getDimSize(0)
             << ") must match species tensor dimension (" << N << ")";

    // Check rate constant.
    if (!inputs[2].isF32())
      return emitOpError("tensor form expects third input to be f32 rate");

    // Check results.
    if (results.size() != 1)
      return emitOpError("tensor form expects exactly 1 result");
    auto resultTensor = dyn_cast<RankedTensorType>(results[0]);
    if (!resultTensor || resultTensor.getRank() != 1 ||
        !isa<acr::SpeciesType>(resultTensor.getElementType()))
      return emitOpError("tensor form result must be tensor<Nx!acr.species>");
    if (resultTensor.getDimSize(0) != N)
      return emitOpError("result tensor dimension must match species "
                         "tensor dimension");
  } else {
    // --- Scalar form ---
    // Expected: (species..., f32 stoich..., f32 rate) -> (species...)
    // At least 1 species, 1 stoich, 1 rate = minimum 3 inputs.
    if (inputs.size() < 3)
      return emitOpError("scalar form expects at least 3 inputs");

    // Count leading species types.
    unsigned numSpecies = 0;
    for (auto ty : inputs) {
      if (isa<acr::SpeciesType>(ty))
        ++numSpecies;
      else
        break;
    }
    if (numSpecies == 0)
      return emitOpError("expected first input(s) to be !acr.species");

    // After species: numSpecies f32 stoich values + 1 f32 rate.
    unsigned expectedInputs = 2 * numSpecies + 1;
    if (inputs.size() != expectedInputs)
      return emitOpError("scalar form with ")
             << numSpecies << " species expects " << expectedInputs
             << " inputs (species + stoich + rate), got " << inputs.size();

    // Verify stoich + rate are f32.
    for (unsigned i = numSpecies; i < inputs.size(); ++i)
      if (!inputs[i].isF32())
        return emitOpError("expected f32 for input #") << i;

    // Verify results are all species and count matches.
    if (results.size() != numSpecies)
      return emitOpError("scalar form expects ")
             << numSpecies << " results, got " << results.size();
    for (auto [i, ty] : llvm::enumerate(results))
      if (!isa<acr::SpeciesType>(ty))
        return emitOpError("expected !acr.species result at index ") << i;
  }

  return success();
}

//===----------------------------------------------------------------------===//
// ACR_CRNGenericOp — verifier
//===----------------------------------------------------------------------===//

LogicalResult CRNGenericOp::verify() {
  // The stoich_matrix attribute is required.
  auto stoichAttr =
      (*this)->getAttrOfType<DenseIntElementsAttr>("stoich_matrix");
  if (!stoichAttr)
    return emitOpError("requires a 'stoich_matrix' dense integer attribute");

  auto stoichType = stoichAttr.getType();
  if (stoichType.getRank() != 2)
    return emitOpError("'stoich_matrix' must be a 2-D tensor");

  int64_t R = stoichType.getDimSize(0); // reactions
  int64_t S = stoichType.getDimSize(1); // species

  auto inputs = getInputs().getTypes();
  auto results = getResults().getTypes();

  if (inputs.empty())
    return emitOpError("expected at least one input");

  // Detect tensor vs scalar form.
  if (auto speciesTensor = dyn_cast<RankedTensorType>(inputs[0])) {
    // --- Tensor form ---
    // Expected: (tensor<Sx!acr.species>, tensor<Rxf32>) ->
    // tensor<Sx!acr.species>
    if (!isa<acr::SpeciesType>(speciesTensor.getElementType()) ||
        speciesTensor.getRank() != 1)
      return emitOpError("tensor form requires first input to be "
                         "tensor<Sx!acr.species>");
    if (speciesTensor.getDimSize(0) != S)
      return emitOpError("species tensor dimension (")
             << speciesTensor.getDimSize(0)
             << ") must match stoich_matrix columns (" << S << ")";

    if (inputs.size() != 2)
      return emitOpError("tensor form expects exactly 2 inputs");

    auto rateTensor = dyn_cast<RankedTensorType>(inputs[1]);
    if (!rateTensor || rateTensor.getRank() != 1 ||
        !rateTensor.getElementType().isF32())
      return emitOpError("tensor form expects second input to be "
                         "tensor<Rxf32>");
    if (rateTensor.getDimSize(0) != R)
      return emitOpError("rate tensor dimension (")
             << rateTensor.getDimSize(0) << ") must match stoich_matrix rows ("
             << R << ")";

    if (results.size() != 1)
      return emitOpError("tensor form expects exactly 1 result");
    auto resultTensor = dyn_cast<RankedTensorType>(results[0]);
    if (!resultTensor || resultTensor.getRank() != 1 ||
        !isa<acr::SpeciesType>(resultTensor.getElementType()))
      return emitOpError("tensor form result must be "
                         "tensor<Sx!acr.species>");
    if (resultTensor.getDimSize(0) != S)
      return emitOpError("result tensor dimension must match "
                         "stoich_matrix columns");
  } else {
    // --- Scalar form ---
    // Expected: (S x !acr.species, R x f32) -> (S x !acr.species)
    unsigned expectedInputs = S + R;
    if ((int64_t)inputs.size() != expectedInputs)
      return emitOpError("scalar form expects ")
             << expectedInputs << " inputs (" << S << " species + " << R
             << " rates), got " << inputs.size();

    for (int64_t i = 0; i < S; ++i)
      if (!isa<acr::SpeciesType>(inputs[i]))
        return emitOpError("expected !acr.species for input #") << i;
    for (int64_t i = S; i < S + R; ++i)
      if (!inputs[i].isF32())
        return emitOpError("expected f32 rate for input #") << i;

    if ((int64_t)results.size() != S)
      return emitOpError("scalar form expects ")
             << S << " species results, got " << results.size();
    for (int64_t i = 0; i < S; ++i)
      if (!isa<acr::SpeciesType>(results[i]))
        return emitOpError("expected !acr.species result at index ") << i;
  }

  return success();
}

//===----------------------------------------------------------------------===//
// Generated op implementations
//===----------------------------------------------------------------------===//

#define GET_OP_CLASSES
#include "ChemComp/Dialect/ACR/ACROps.cpp.inc"
