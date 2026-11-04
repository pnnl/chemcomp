//===- ODEOps.cpp - ODE dialect operations --------------------------------===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#include "ChemComp/Dialect/ODE/ODEOps.h"
#include "ChemComp/Dialect/ODE/ODESystem.h"

#include "mlir/IR/Builders.h"
#include "mlir/IR/OpImplementation.h"

using namespace mlir;
using namespace mlir::cc::ode;

//===----------------------------------------------------------------------===//
// ode.system
//===----------------------------------------------------------------------===//

void SystemOp::build(OpBuilder &builder, OperationState &state,
                     StringRef name) {
  state.addAttribute(getSymNameAttrName(state.name),
                     builder.getStringAttr(name));
  Region *body = state.addRegion();
  body->emplaceBlock();
}

LogicalResult SystemOp::verifyRegions() {
  for (Operation &op : getBody().front()) {
    if (!isa<StateOp, ParamOp, TermOp>(op))
      return op.emitOpError("is not allowed in an ode.system; expected "
                            "ode.state, ode.param or ode.term");
  }
  return success();
}

//===----------------------------------------------------------------------===//
// ode.state:  ode.state @x [initial <number>] attr-dict
//===----------------------------------------------------------------------===//

ParseResult StateOp::parse(OpAsmParser &parser, OperationState &result) {
  StringAttr name;
  if (parser.parseSymbolName(name, getSymNameAttrName(result.name),
                             result.attributes))
    return failure();
  if (succeeded(parser.parseOptionalKeyword("initial"))) {
    double value;
    if (parseNumber(parser, value))
      return failure();
    result.addAttribute(getInitialAttrName(result.name),
                        parser.getBuilder().getF64FloatAttr(value));
  }
  return parser.parseOptionalAttrDict(result.attributes);
}

void StateOp::print(OpAsmPrinter &p) {
  p << ' ';
  p.printSymbolName(getSymName());
  if (FloatAttr initial = getInitialAttr()) {
    p << " initial ";
    printNumber(p, initial.getValueAsDouble());
  }
  p.printOptionalAttrDict((*this)->getAttrs(),
                          {getSymNameAttrName(), getInitialAttrName()});
}

//===----------------------------------------------------------------------===//
// ode.param:  ode.param @k [= <number>] attr-dict
//===----------------------------------------------------------------------===//

ParseResult ParamOp::parse(OpAsmParser &parser, OperationState &result) {
  StringAttr name;
  if (parser.parseSymbolName(name, getSymNameAttrName(result.name),
                             result.attributes))
    return failure();
  if (succeeded(parser.parseOptionalEqual())) {
    double value;
    if (parseNumber(parser, value))
      return failure();
    result.addAttribute(getValueAttrName(result.name),
                        parser.getBuilder().getF64FloatAttr(value));
  }
  return parser.parseOptionalAttrDict(result.attributes);
}

void ParamOp::print(OpAsmPrinter &p) {
  p << ' ';
  p.printSymbolName(getSymName());
  if (FloatAttr value = getValueAttr()) {
    p << " = ";
    printNumber(p, value.getValueAsDouble());
  }
  p.printOptionalAttrDict((*this)->getAttrs(),
                          {getSymNameAttrName(), getValueAttrName()});
}

LogicalResult ParamOp::verify() {
  if (FloatAttr value = getValueAttr(); value && value.getValueAsDouble() < 0)
    return emitOpError("parameter values are rate constants and must not be "
                       "negative; move the sign into the terms");
  return success();
}

//===----------------------------------------------------------------------===//
// ode.term:  ode.term @x += <number> [params [...]] [states [...]] attr-dict
//===----------------------------------------------------------------------===//

ParseResult TermOp::parse(OpAsmParser &parser, OperationState &result) {
  Builder &builder = parser.getBuilder();
  FlatSymbolRefAttr state;
  double coeff;
  if (parser.parseAttribute(state, getStateAttrName(result.name),
                            result.attributes) ||
      parser.parsePlus() || parser.parseEqual() || parseNumber(parser, coeff))
    return failure();
  result.addAttribute(getCoeffAttrName(result.name),
                      builder.getF64FloatAttr(coeff));

  auto parseList = [&](StringRef keyword, StringAttr attrName) {
    ArrayAttr list = builder.getArrayAttr({});
    if (succeeded(parser.parseOptionalKeyword(keyword)) &&
        parser.parseAttribute(list))
      return failure();
    result.addAttribute(attrName, list);
    return success();
  };
  if (failed(parseList("params", getParamsAttrName(result.name))) ||
      failed(parseList("states", getStatesAttrName(result.name))))
    return failure();
  return parser.parseOptionalAttrDict(result.attributes);
}

void TermOp::print(OpAsmPrinter &p) {
  p << ' ';
  p.printAttributeWithoutType(getStateAttr());
  p << " += ";
  printNumber(p, getCoeff().convertToDouble());
  if (!getParams().empty())
    p << " params " << getParams();
  if (!getStates().empty())
    p << " states " << getStates();
  p.printOptionalAttrDict((*this)->getAttrs(),
                          {getStateAttrName(), getCoeffAttrName(),
                           getParamsAttrName(), getStatesAttrName()});
}

LogicalResult TermOp::verify() {
  Operation *system = (*this)->getParentOp();
  auto lookup = [&](FlatSymbolRefAttr ref) {
    return SymbolTable::lookupSymbolIn(system, ref.getAttr());
  };
  if (!isa_and_nonnull<StateOp>(lookup(getStateAttr())))
    return emitOpError("derivative of '")
           << getState() << "', which is not an ode.state of this system";
  for (Attribute attr : getParams()) {
    auto ref = dyn_cast<FlatSymbolRefAttr>(attr);
    if (!ref || !isa_and_nonnull<ParamOp>(lookup(ref)))
      return emitOpError("'") << attr << "' is not an ode.param of this system";
  }
  for (Attribute attr : getStates()) {
    auto ref = dyn_cast<FlatSymbolRefAttr>(attr);
    if (!ref || !isa_and_nonnull<StateOp>(lookup(ref)))
      return emitOpError("'") << attr << "' is not an ode.state of this system";
  }
  return success();
}

#define GET_OP_CLASSES
#include "ChemComp/Dialect/ODE/ODEOps.cpp.inc"
