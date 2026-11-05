//===- Passes.td - ODE to ACR conversion passes -------------*- tablegen -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#ifndef CHEMCOMP_CONVERSION_ODETOACR_PASSES
#define CHEMCOMP_CONVERSION_ODETOACR_PASSES

include "mlir/Pass/PassBase.td"

def ConvertODEToACR : Pass<"convert-ode-to-acr", "::mlir::ModuleOp"> {
  let summary = "Compose a reaction network from the terms of an ODE system";
  let description = [{
    Turns each `ode.system` into an `acr.network` whose mass-action kinetics
    are exactly the system's equations.

    The system must be kinetic (no negative cross-effects; run
    `ode-dual-rail` first, as the `ode-to-crn` pipeline does). Then:

    1. **Composable rules.** A fixed set of reaction motifs (ODERules.td:
       binding, dissociation, predation, conversion, degradation, ...) is
       tried in priority order. A rule matches when, for some assignment of
       its species to states and some rate constant k, every term its
       reaction contributes under mass action is present with exactly that
       coefficient and parameters. The rule then claims those terms and emits
       its reaction. Rules are applied greedily until none matches.
    2. **Canonical fallback** (Hars & Toth). Each term no rule claimed becomes
       one reaction: `+c*x^y` in der(x_i) gives `y -> y + X_i`, `-c*x^y` gives
       `y -> y - X_i`, both with rate c. Disable it with
       `canonical-fallback=false` to compose with rules only; unclaimed terms
       are then an error.

    Every reaction records its `rule`, its rate expression (`rate_name`), its
    numeric `rate` when all parameters have values, and the terms it accounts
    for (`from`).
  }];

  let dependentDialects = ["::mlir::cc::acr::ACRDialect"];

  let options = [
    Option<"canonicalFallback", "canonical-fallback", "bool",
           /*default=*/"true",
           "Cover terms that no rule claims with the canonical Hungarian-Lemma "
           "reaction (otherwise they are an error)">,
  ];
}

#endif // CHEMCOMP_CONVERSION_ODETOACR_PASSES
