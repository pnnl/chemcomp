//===- ODEDialect.td - ODE dialect definition ---------------*- tablegen -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#ifndef ODE_DIALECT
#define ODE_DIALECT

include "mlir/IR/OpBase.td"

def ODE_Dialect : Dialect {
  let name = "ode";
  let summary = "Polynomial systems of ordinary differential equations";
  let description = [{
    The `ode` dialect holds explicit polynomial ODE systems

        x_i' = sum_t c_t * prod_p p^a_p * prod_j x_j^b_j

    with one `ode.term` operation per monomial. It is the input of the
    composable ODE-to-reaction-network flow: passes rewrite the system term by
    term (normalization, dual-rail encoding), and `convert-ode-to-acr` turns
    the terms into mass-action reactions.

    The dialect has no MARCO dependency. Models written in Modelica reach it
    through `convert-bmodelica-to-ode` when ChemComp is built with MARCO.
  }];
  let cppNamespace = "::mlir::cc::ode";
}

class ODE_Op<string mnemonic, list<Trait> traits = []> :
    Op<ODE_Dialect, mnemonic, traits>;

#endif // ODE_DIALECT
