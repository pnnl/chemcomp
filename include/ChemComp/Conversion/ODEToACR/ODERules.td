//===- ODERules.td - Composable ODE-to-reaction rules -------*- tablegen -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// The fixed rule set of convert-ode-to-acr. `chemcomp-tblgen -gen-ode-rules`
// turns it into a C++ table; there is no runtime rule loading.
//
// A rule is a single reaction motif, written as reactants and products over
// rule-local species names; a name listed n times has count n. The ODE terms
// a rule claims are derived from the motif by mass action: the reaction
// y -> y' with rate k adds k * (y'_s - y_s) * x^y to der(s) for each species
// s. So a rule cannot disagree with the reaction it emits.
//
// Rules are tried by descending priority (greedy). Larger motifs, which claim
// more terms, come first.
//
//===----------------------------------------------------------------------===//

#ifndef CHEMCOMP_CONVERSION_ODETOACR_ODERULES
#define CHEMCOMP_CONVERSION_ODETOACR_ODERULES

class ODERule<string name, int priority, list<string> reactants,
              list<string> products, string summary> {
  string Name = name;
  int Priority = priority;
  list<string> Reactants = reactants;
  list<string> Products = products;
  string Summary = summary;
}

// Three-term motifs.
def Binding : ODERule<"binding", 100, ["A", "B"], ["C"],
    "two species bind into a complex">;
def Dissociation : ODERule<"dissociation", 95, ["C"], ["A", "B"],
    "a complex falls apart into two species">;

// Two-term motifs.
def Predation : ODERule<"predation", 90, ["X", "Y"], ["Y", "Y"],
    "Y converts X into another Y (predation, infection)">;
def CubicAutocatalysis : ODERule<"cubic-autocatalysis", 85,
    ["X", "X", "Y"], ["X", "X", "X"],
    "two X convert Y into X (Brusselator step)">;
def Conversion : ODERule<"conversion", 80, ["A"], ["B"],
    "A turns into B (isomerization, recovery)">;
def Annihilation : ODERule<"annihilation", 75, ["X", "Y"], [],
    "X and Y destroy each other (dual-rail pairs)">;

// One-term motifs.
def CatalyticProduction : ODERule<"catalytic-production", 60, ["Y"],
    ["Y", "X"], "Y produces X without being consumed">;
def Autocatalysis : ODERule<"autocatalysis", 55, ["X"], ["X", "X"],
    "X replicates">;
def SecondOrderDecay : ODERule<"second-order-decay", 50, ["X", "X"], ["X"],
    "two X meet and one disappears">;
def Degradation : ODERule<"degradation", 45, ["X"], [],
    "X decays">;
def Production : ODERule<"production", 40, [], ["X"],
    "X is produced at a constant rate">;

#endif // CHEMCOMP_CONVERSION_ODETOACR_ODERULES
