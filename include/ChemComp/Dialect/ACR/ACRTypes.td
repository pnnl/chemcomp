//===- ACRTypes.td - ACR dialect types -------*- tablegen -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#ifndef ACR_TYPES
#define ACR_TYPES

include "mlir/IR/AttrTypeBase.td"
include "mlir/IR/BuiltinTypeInterfaces.td"
include "ChemComp/Dialect/ACR/ACRDialect.td"

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
    individual values or grouped into tensors for bulk operations.

    Example:
    ```mlir
    // Single species argument
    acr.reaction @r(%h2: !acr.species, %o2: !acr.species) -> !acr.species
    // Tensor of species
    acr.reaction @r(%sv: tensor<3x!acr.species>) -> tensor<3x!acr.species>
    ```
  }];
  // Zero parameters — prints as `!acr.species` with no angle brackets.
}

// A 1-D ranked tensor whose element type is !acr.species.
def ACR_SpeciesTensor : 1DTensorOf<[ACR_SpeciesType]>;

// Matches either a scalar !acr.species or a 1-D tensor of species.
def ACR_SpeciesLike : AnyTypeOf<[ACR_SpeciesType, ACR_SpeciesTensor],
                                 "species scalar or tensor">;

#endif // ACR_TYPES
