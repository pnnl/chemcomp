//===- ODEPasses.td - ODE dialect passes --------------------*- tablegen -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#ifndef ODE_PASSES
#define ODE_PASSES

include "mlir/Pass/PassBase.td"

def ODENormalize : Pass<"ode-normalize", "::mlir::ModuleOp"> {
  let summary = "Sum like terms and put ODE systems in canonical order";
  let description = [{
    In every `ode.system`, terms with the same state, monomial and parameters
    are summed into one; terms that cancel are dropped; the remaining terms
    are ordered by state (declaration order), then by monomial. Two systems
    that are equal as polynomials normalize to the same IR.
  }];
}

def ODECheckKinetic : Pass<"ode-check-kinetic", "::mlir::ModuleOp"> {
  let summary = "Report negative cross-effects (Hungarian Lemma)";
  let description = [{
    A polynomial ODE system is the mass-action system of some reaction network
    if and only if every negative term in `der(x)` contains `x` (Hars & Toth
    1981). This pass reports each term that breaks the condition, e.g.
    `der(y) has -x, which does not contain y`, as a remark, or as an error with
    `require`. Systems without such terms get a remark saying they are
    kinetic.
  }];
  let options = [
    Option<"require", "require", "bool", /*default=*/"false",
           "Report negative cross-effects as errors">,
  ];
}

def ODEDualRail : Pass<"ode-dual-rail", "::mlir::ModuleOp"> {
  let summary = "Encode states as differences x = x_p - x_m of two rails";
  let description = [{
    Rewrites a state `x` as the difference of two non-negative states,
    `x = x_p - x_m`. Substituting into the equations and expanding gives only
    terms of the form `c * monomial`; each is routed to one rail so that every
    negative term of a rail contains that rail:

      - a positive term goes to `der(x_p)`, or, as a negative term, to
        `der(x_m)` when it contains `x_m` (it then degrades the minus rail);
      - a negative term stays in `der(x_p)` when it contains `x_p`, otherwise
        it goes to `der(x_m)` with the opposite sign.

    Both rails also get `-eta * x_p * x_m` (the reaction `x_p + x_m -> 0`),
    which cancels in `der(x_p) - der(x_m)`. So `x_p - x_m` follows the original
    equation exactly, for any `eta`; `eta` only keeps the rails small.

    With `mode=auto` (the default) only states that need it are split: the
    targets of negative cross-effects, repeated until the system is kinetic
    (splitting one state can create cross-effects in others). `mode=all`
    splits every state. `states` adds states to split in either mode.

    The rails are named `<x>_p` and `<x>_m` and carry `rail_of` and
    `rail_sign`; the system records the split states in `dual_rail`.
  }];
  let options = [
    Option<"mode", "mode", "std::string", /*default=*/"\"auto\"",
           "Which states to split: 'auto' (only where needed) or 'all'">,
    ListOption<"extraStates", "states", "std::string",
               "Additional states to split">,
    Option<"annihilationRate", "annihilation-rate", "double",
           /*default=*/"1.0",
           "Rate eta of the annihilation x_p + x_m -> 0">,
  ];
}

#endif // ODE_PASSES
