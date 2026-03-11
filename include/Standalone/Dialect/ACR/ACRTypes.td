//===- ACRTypes.td - ACR dialect types -------*- tablegen -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#ifndef ACR_TYPES
#define ACR_TYPES

include "mlir/IR/AttrTypeBase.td"
include "Standalone/Dialect/ACR/ACRDialect.td"

//===----------------------------------------------------------------------===//
// ACR type definitions
//===----------------------------------------------------------------------===//

class ACR_Type<string name, string typeMnemonic, list<Trait> traits = []>
    : TypeDef<ACR_Dialect, name, traits> {
  let mnemonic = typeMnemonic;
}

// The species type represents the concentration of a chemical species.
// It is a zero-parameter type that prints as `!acr.species`.
def ACR_SpeciesType : ACR_Type<"Species", "species"> {
  let summary = "Type representing the concentration of a chemical species";
  let description = [{
    The `!acr.species` type represents the dynamic concentration of a chemical
    species in an abstract chemical reaction. Concentrations can be passed as
    individual values or grouped into vectors for bulk operations.

    Example:
    ```mlir
    // Single species argument
    acr.reaction @r(%h2: !acr.species, %o2: !acr.species) -> !acr.species
    // Vector of species
    acr.reaction @r(%sv: vector<3x!acr.species>) -> vector<3x!acr.species>
    ```
  }];
  // Zero parameters — prints as `!acr.species` with no angle brackets.
}

#endif // ACR_TYPES
