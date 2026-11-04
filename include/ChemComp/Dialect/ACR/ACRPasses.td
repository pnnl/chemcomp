//===- ACRPsss.td - ACR dialect passes -------*- tablegen -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#ifndef ACR_PASS
#define ACR_PASS

include "mlir/Pass/PassBase.td"

def ACRSwitchBarFoo: Pass<"acr-switch-bar-foo", "::mlir::ModuleOp"> {
  let summary = "Switches the name of a FuncOp named `bar` to `foo` and folds.";
  let description = [{
    Switches the name of a FuncOp named `bar` to `foo` and folds.
    ```
    func.func @bar() {
      return
    }
    // Gets transformed to:
    func.func @foo() {
      return
    }
    ```
  }];
}

def ACRNetworkMerge : Pass<"acr-network-merge", "::mlir::ModuleOp"> {
  let summary = "Merge reactions of an acr.network without changing its ODEs";
  let description = [{
    Rewrites the reactions of every `acr.network` into fewer reactions with
    the same mass-action ODEs. Only the net reaction vector of each reactant
    complex matters (`sum_k rate_k * (products_k - reactants)`), so:

      - identical reactions are summed into one (rates add);
      - two reactions with the same reactants and the same rate (equal value
        and equal `rate_name`) become one reaction whose change is the sum of
        both, provided no product count becomes negative. For example
        `X + Y -> Y` and `X + Y -> X + 2Y`, both at rate b, become
        `X + Y -> 2Y`;
      - reactions whose products equal their reactants are dropped.

    It mostly tidies the output of the canonical fallback of
    `convert-ode-to-acr`; reactions emitted by composable rules usually have
    nothing left to merge.
  }];
}

def ACRNetworkToCRN : Pass<"acr-network-to-crn", "::mlir::ModuleOp"> {
  let summary = "Lower acr.network to acr.crn with an acr.crn_generic";
  let description = [{
    Replaces every `acr.network @N` with `acr.crn @N` holding one
    `acr.crn_generic`, the form `crn-db-match` and `mlir-to-json` consume.
    Species columns and reaction rows keep the network's order. Besides
    `stoich_matrix` (net), the op carries `reactant_matrix`, `product_matrix`,
    `species_symbols`, `ode_state_names`, `ode_state_signs`, `rate_names`,
    `rules`, and, when every value is known, `rate_values` and
    `initial_values`.
  }];
}

#endif // ACR_PASS
