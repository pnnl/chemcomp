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

def ACR_CustomType : ACR_Type<"Custom", "custom"> {
    let summary = "ACR custom type";
    let description = "Custom type in acr dialect";
    let parameters = (ins StringRefParameter<"the custom value">:$value);
    let assemblyFormat = "`<` $value `>`";
}

def ACR_CRType : ACR_Type<"CR", "cr"> {
    let summary = "ACR chemical reaction type";
    let description = "Type representing the number of species in a chemical reaction";
    let parameters = (ins "uint64_t":$numSpecies);
    let assemblyFormat = "`<` $numSpecies `>`";
}

#endif // ACR_TYPES
