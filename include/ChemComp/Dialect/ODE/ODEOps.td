//===- ODEOps.td - ODE dialect operations -------------------*- tablegen -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#ifndef ODE_OPS
#define ODE_OPS

include "ChemComp/Dialect/ODE/ODEDialect.td"
include "mlir/IR/SymbolInterfaces.td"
include "mlir/IR/RegionKindInterface.td"

//===----------------------------------------------------------------------===//
// ode.system
//===----------------------------------------------------------------------===//

def ODE_SystemOp : ODE_Op<"system", [
    Symbol, SymbolTable, IsolatedFromAbove, NoTerminator, SingleBlock,
    NoRegionArguments, GraphRegionNoTerminator
]> {
  let summary = "A polynomial ODE system";
  let description = [{
    Holds the states, parameters and terms of one ODE system. The order of the
    `ode.state` operations is the order of the states everywhere downstream
    (species columns, JSON output).

    ```mlir
    // der(x) = -k*x
    ode.system @Decay {
      ode.state @x initial 1
      ode.param @k = 0.5
      ode.term @x += -1 params [@k] states [@x]
    }
    ```
  }];

  let arguments = (ins SymbolNameAttr:$sym_name);
  let regions = (region SizedRegion<1>:$body);
  let assemblyFormat = "$sym_name attr-dict-with-keyword $body";

  let skipDefaultBuilders = 1;
  let builders = [OpBuilder<(ins "::llvm::StringRef":$name)>];
  let hasRegionVerifier = 1;
}

//===----------------------------------------------------------------------===//
// ode.state
//===----------------------------------------------------------------------===//

def ODE_StateOp : ODE_Op<"state", [Symbol, HasParent<"SystemOp">]> {
  let summary = "A state variable of an ODE system";
  let description = [{
    Declares a state `x`, optionally with its initial value. After
    `ode-dual-rail`, the two rails of a state record the state they encode and
    their sign: `x = x_p - x_m`.

    ```mlir
    ode.state @x initial 2
    ode.state @x_p initial 2 {rail_of = "x", rail_sign = 1 : i64}
    ```
  }];

  let arguments = (ins SymbolNameAttr:$sym_name,
                       OptionalAttr<F64Attr>:$initial,
                       OptionalAttr<StrAttr>:$rail_of,
                       OptionalAttr<I64Attr>:$rail_sign);
  let hasCustomAssemblyFormat = 1;
}

//===----------------------------------------------------------------------===//
// ode.param
//===----------------------------------------------------------------------===//

def ODE_ParamOp : ODE_Op<"param", [Symbol, HasParent<"SystemOp">]> {
  let summary = "A parameter (rate constant) of an ODE system";
  let description = [{
    Declares a parameter, optionally with its value. Parameters are rate
    constants and must not be negative.

    ```mlir
    ode.param @k = 0.5
    ode.param @kd
    ```
  }];

  let arguments = (ins SymbolNameAttr:$sym_name, OptionalAttr<F64Attr>:$value);
  let hasCustomAssemblyFormat = 1;
  let hasVerifier = 1;
}

//===----------------------------------------------------------------------===//
// ode.term
//===----------------------------------------------------------------------===//

def ODE_TermOp : ODE_Op<"term", [HasParent<"SystemOp">]> {
  let summary = "One monomial of a state's derivative";
  let description = [{
    Adds `coeff * prod(params) * prod(states)` to the derivative of `state`.
    A symbol listed n times is raised to the power n, so `states [@x, @x, @y]`
    is `x^2*y`. An empty `states` list is a constant term.

    ```mlir
    // der(x) += -b*x*y
    ode.term @x += -1 params [@b] states [@x, @y]
    // der(x) += 2
    ode.term @x += 2
    ```
  }];

  let arguments = (ins FlatSymbolRefAttr:$state,
                       F64Attr:$coeff,
                       SymbolRefArrayAttr:$params,
                       SymbolRefArrayAttr:$states);
  let hasCustomAssemblyFormat = 1;
  let hasVerifier = 1;
}

#endif // ODE_OPS
