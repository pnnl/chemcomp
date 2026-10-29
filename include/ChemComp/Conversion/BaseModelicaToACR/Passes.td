//===- Passes.td - Base Modelica to ACR conversion passes --*- tablegen -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// These passes depend on the MARCO Base Modelica dialect, so they live in their
// own TableGen file and are only built with CHEMCOMP_ENABLE_MARCO=ON.
//
//===----------------------------------------------------------------------===//

#ifndef CHEMCOMP_CONVERSION_BASEMODELICATOACR_PASSES
#define CHEMCOMP_CONVERSION_BASEMODELICATOACR_PASSES

include "mlir/Pass/PassBase.td"

def ConvertBaseModelicaToACR
    : Pass<"convert-bmodelica-to-acr", "::mlir::ModuleOp"> {
  let summary = "Map Base Modelica ODE models onto known chemical reaction "
                "networks";
  let description = [{
    Converts each `bmodelica.model` whose dynamic equations form an explicit,
    scalar, polynomial ODE system into an `acr.crn` holding an
    `acr.crn_generic` with the network's net stoichiometry.

    Every equation of the model must define the derivative of one state
    variable (`der(x) = f(...)`, either before or after
    `derivatives-materialization`). The right-hand sides are normalized into
    polynomials over the states, with parameters and numeric constants as
    coefficients. The polynomial system is then compared, up to renaming of
    states and rate constants, against a fixed table of known networks whose
    mass-action kinetics produce it (decay, birth-death, reversible
    isomerization, Lotka-Volterra).

    The emitted `acr.crn_generic` carries, besides `stoich_matrix`:
      - `crn_rule`: name of the matched rule;
      - `ode_state_names`: Base Modelica state for each species column;
      - `reactant_matrix` / `product_matrix`: full stoichiometry, which the
        net `stoich_matrix` loses for catalysts (e.g. X -> 2X);
      - `rate_names`: rate constant expression of each reaction;
      - `rate_values`: numeric rate constants, when every parameter has a
        constant binding or start value.

    The network is named after the model; with `keep-model`, the model keeps
    its name and the network is named `<model>_crn`. Models that match no
    rule are an error, unless `allow-unmatched` is set.
  }];

  let dependentDialects = ["::mlir::cc::acr::ACRDialect"];

  let options = [
    Option<"keepModel", "keep-model", "bool", /*default=*/"false",
           "Keep the converted bmodelica.model operations">,
    Option<"allowUnmatched", "allow-unmatched", "bool", /*default=*/"false",
           "Warn instead of failing on models that match no known network">,
  ];
}

#endif // CHEMCOMP_CONVERSION_BASEMODELICATOACR_PASSES
