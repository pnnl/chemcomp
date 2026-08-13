//===- ACRDialect.td - ACR dialect -----------*- tablegen -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#ifndef ACR_DIALECT
#define ACR_DIALECT

include "mlir/IR/OpBase.td"

//===----------------------------------------------------------------------===//
// ACR dialect definition.
//===----------------------------------------------------------------------===//

def ACR_Dialect : Dialect {
    let name = "acr";
    let summary = "A dialect for Absract Chemical Reaction Networks.";
    let description = [{
        TODO
    }];
    let cppNamespace = "::mlir::acr";

    let useDefaultTypePrinterParser = 1;
    let extraClassDeclaration = [{
        void registerTypes();
    }];
}

//===----------------------------------------------------------------------===//
// Base acr operation definition.
//===----------------------------------------------------------------------===//

class ACR_Op<string mnemonic, list<Trait> traits = []> :
        Op<ACR_Dialect, mnemonic, traits>;

#endif // ACR_DIALECT
