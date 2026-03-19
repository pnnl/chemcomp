//===- ACROps.td - ACR dialect ops -----------*- tablegen -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#ifndef ACR_OPS
#define ACR_OPS

include "Standalone/Dialect/ACR/ACRTypes.td"
include "mlir/Interfaces/FunctionInterfaces.td"
include "mlir/Interfaces/CallInterfaces.td"
include "mlir/Interfaces/ControlFlowInterfaces.td"
include "mlir/Interfaces/SideEffectInterfaces.td"
include "mlir/IR/SymbolInterfaces.td"

//===----------------------------------------------------------------------===//
// acr.stoich_coeff
//
// Produces a stoichiometric coefficient as an f32 value.
// Negative coefficients indicate reactants, positive indicate products.
//
// Example:
//   %h2_stoich = acr.stoich_coeff -2 : f32
//===----------------------------------------------------------------------===//

def ACR_StoichCoeffOp : ACR_Op<"stoich_coeff", [Pure]> {
  let summary = "Stoichiometric coefficient constant";
  let description = [{
    Produces a stoichiometric coefficient as an f32 scalar. Negative values
    denote reactants (consumed), positive values denote products (produced).

    Example:
    ```mlir
    %h2_stoich  = acr.stoich_coeff -2 : f32  // 2 H2 consumed
    %h2o_stoich = acr.stoich_coeff  2 : f32  // 2 H2O produced
    ```
  }];

  let arguments = (ins SI64Attr:$value);
  let results = (outs F32:$result);

  // Custom format: `acr.stoich_coeff -2 : f32`
  // The SI64Attr prints without a type suffix for i64 integers.
  let assemblyFormat = "$value `:` type($result) attr-dict";
}

//===----------------------------------------------------------------------===//
// acr.react
//
// Applies a reaction to a set of species given stoichiometry and a rate.
// Two calling conventions are supported (distinguished by operand types):
//
//  Scalar  – individual !acr.species + individual f32 stoich values + f32 rate
//  Tensor  – tensor<Nx!acr.species>  + tensor<Nxf32> stoich        + f32 rate
//
// Example (scalar):
//   %h2_o, %o2_o = acr.react(%h2, %o2, %hs, %os, %k)
//       : (!acr.species, !acr.species, f32, f32, f32) -> (!acr.species, !acr.species)
//
// Example (tensor):
//   %sv_out = acr.react(%sv, %stoich, %k)
//       : (tensor<3x!acr.species>, tensor<3xf32>, f32) -> tensor<3x!acr.species>
//===----------------------------------------------------------------------===//

def ACR_ReactOp : ACR_Op<"react"> {
  let summary = "Apply a chemical reaction to species concentrations";
  let description = [{
    Applies a single reaction step to a set of species, updating their
    concentrations according to the stoichiometric coefficients and rate.

    The operation accepts either individual `!acr.species` SSA values together
    with matching individual `f32` stoichiometric coefficients and a `f32`
    rate constant (scalar form), or a `tensor<Nx!acr.species>` species tensor
    together with a `tensor<Nxf32>` stoichiometry tensor and a `f32` rate
    (tensor form).
  }];

  let arguments = (ins Variadic<AnyType>:$inputs);
  let results = (outs Variadic<AnyType>:$results);

  let assemblyFormat = [{
    `(` $inputs `)` `:` functional-type($inputs, results) attr-dict
  }];

  let hasVerifier = 1;
}

//===----------------------------------------------------------------------===//
// acr.reaction
//
// A module-level, function-like operation that represents a single abstract
// chemical reaction. Arguments are the species concentrations and rate
// constant; the body contains stoichiometry via acr.stoich_coeff and the
// acr.react call, and terminates with a func.return.
//
// Supports an optional `private` visibility keyword (same as func.func).
//
// Example (scalar definition):
//   acr.reaction @r1(%h2: !acr.species, %o2: !acr.species, %k: f32)
//       -> (!acr.species, !acr.species) { ... }
//
// Example (tensor definition):
//   acr.reaction @r1(%sv: tensor<3x!acr.species>, %k: f32)
//       -> tensor<3x!acr.species> { ... }
//
// Example (external declaration):
//   acr.reaction private @r1(%h2: !acr.species, %o2: !acr.species, %k: f32)
//       -> (!acr.species, !acr.species)
//===----------------------------------------------------------------------===//

def ACR_ReactionOp : ACR_Op<"reaction", [
    AutomaticAllocationScope,
    FunctionOpInterface,
    IsolatedFromAbove
]> {
  let summary = "Abstract chemical reaction (function-like)";
  let description = [{
    An `acr.reaction` operation defines (or declares) a single abstract
    chemical reaction. Its arguments are the concentrations of the participating
    species and the rate constant; its results are the updated concentrations.
    The body contains `acr.stoich_coeff` and `acr.react` operations and is
    terminated by a `func.return`.

    An external declaration (no body) uses the `private` visibility keyword
    and omits the `{...}` block.
  }];

  let arguments = (ins
    SymbolNameAttr:$sym_name,
    TypeAttrOf<FunctionType>:$function_type,
    OptionalAttr<StrAttr>:$sym_visibility,
    OptionalAttr<DictArrayAttr>:$arg_attrs,
    OptionalAttr<DictArrayAttr>:$res_attrs
  );
  let regions = (region AnyRegion:$body);

  let extraClassDeclaration = [{
    //===------------------------------------------------------------------===//
    // FunctionOpInterface / CallableOpInterface Methods
    //===------------------------------------------------------------------===//
    ArrayRef<Type> getArgumentTypes() { return getFunctionType().getInputs(); }
    ArrayRef<Type> getResultTypes()   { return getFunctionType().getResults(); }
    ::mlir::Region *getCallableRegion() {
      return isExternal() ? nullptr : &getBody();
    }

    //===------------------------------------------------------------------===//
    // SymbolOpInterface Methods
    //===------------------------------------------------------------------===//
    bool isDeclaration() { return isExternal(); }

    /// Returns true when the reaction has no body (external declaration).
    bool isExternal() { return getBody().empty(); }
  }];

  let hasCustomAssemblyFormat = 1;
}

//===----------------------------------------------------------------------===//
// acr.reaction_instance
//
// Calls an acr.reaction (or acr.crn) by symbol reference.
// Mirrors func.call but lives inside acr.crn bodies.
//
// Example:
//   %h2_o, %o2_o, %h2o_o =
//       acr.reaction_instance @r1(%h2, %o2, %h2o, %k)
//           : (!acr.species, !acr.species, !acr.species, f32)
//          -> (!acr.species, !acr.species, !acr.species)
//===----------------------------------------------------------------------===//

def ACR_ReactionInstanceOp : ACR_Op<"reaction_instance", [
    CallOpInterface
]> {
  let summary = "Instantiate (call) a named acr.reaction or acr.crn";
  let description = [{
    `acr.reaction_instance` calls an `acr.reaction` or `acr.crn` operation
    identified by a symbol reference. It passes the species concentrations
    and rate constant as arguments and receives the updated concentrations
    as results.
  }];

  let arguments = (ins FlatSymbolRefAttr:$callee,
                       Variadic<AnyType>:$operands);
  let results = (outs Variadic<AnyType>:$results);

  let extraClassDeclaration = [{
    operand_range getArgOperands() { return getOperands(); }
    MutableOperandRange getArgOperandsMutable() {
      return getOperandsMutable();
    }
    CallInterfaceCallable getCallableForCallee() {
      return (*this)->getAttrOfType<SymbolRefAttr>("callee");
    }
    void setCalleeFromCallable(CallInterfaceCallable callee) {
      (*this)->setAttr("callee", callee.get<SymbolRefAttr>());
    }
  }];

  let assemblyFormat = [{
    $callee `(` $operands `)` attr-dict `:` functional-type($operands, results)
  }];
}

//===----------------------------------------------------------------------===//
// acr.crn
//
// A module-level, function-like operation that defines a Chemical Reaction
// Network (CRN) by composing multiple acr.reaction_instance calls.
//
// Example (scalar):
//   acr.crn @net(%h2: !acr.species, %o2: !acr.species, %k: f32)
//       -> (!acr.species, !acr.species) {
//     %h2_o, %o2_o = acr.reaction_instance @r1(%h2, %o2, %k)
//         : (!acr.species, !acr.species, f32) -> (!acr.species, !acr.species)
//     acr.return %h2_o, %o2_o : !acr.species, !acr.species
//   }
//
// Example (tensor):
//   acr.crn @net(%sv: tensor<3x!acr.species>, %k: f32)
//       -> tensor<3x!acr.species> {
//     ...
//     acr.return %sv_out : tensor<3x!acr.species>
//   }
//===----------------------------------------------------------------------===//

def ACR_CRNOp : ACR_Op<"crn", [
    AutomaticAllocationScope,
    FunctionOpInterface,
    IsolatedFromAbove
]> {
  let summary = "Chemical reaction network (function-like)";
  let description = [{
    An `acr.crn` operation defines a Chemical Reaction Network by composing
    one or more `acr.reaction_instance` calls. It has the same function-like
    structure as `acr.reaction`.
  }];

  let arguments = (ins
    SymbolNameAttr:$sym_name,
    TypeAttrOf<FunctionType>:$function_type,
    OptionalAttr<StrAttr>:$sym_visibility,
    OptionalAttr<DictArrayAttr>:$arg_attrs,
    OptionalAttr<DictArrayAttr>:$res_attrs
  );
  let regions = (region AnyRegion:$body);

  let extraClassDeclaration = [{
    ArrayRef<Type> getArgumentTypes() { return getFunctionType().getInputs(); }
    ArrayRef<Type> getResultTypes()   { return getFunctionType().getResults(); }
    ::mlir::Region *getCallableRegion() {
      return isExternal() ? nullptr : &getBody();
    }
    bool isDeclaration() { return isExternal(); }
    bool isExternal() { return getBody().empty(); }
  }];

  let hasCustomAssemblyFormat = 1;
}

//===----------------------------------------------------------------------===//
// acr.crn_generic
//
// A compact representation of a CRN using a stoichiometric matrix attribute.
// Intended as the canonicalized form of acr.crn after lowering.
//
// The stoich_matrix is a 2-D integer tensor where each row is a reaction and
// each column is a species (reactions × species, i.e. ChemComp transposed
// convention).
//
// Example (tensor form, inside acr.crn body):
//   %out = acr.crn_generic
//       {stoich_matrix = dense<[[-2,-1,2,0,0],[0,-1,0,-1,1]]>
//           : tensor<2x5xi32>}
//       (%species_vec, %rate_vec) :
//       (tensor<5x!acr.species>, tensor<2xf32>) -> tensor<5x!acr.species>
//
// Example (scalar form):
//   %out0, %out1, %out2 = acr.crn_generic
//       {stoich_matrix = dense<[[-2, -1, 2], [0, -1, 1]]> : tensor<2x3xi32>}
//       (%s0, %s1, %s2, %r0, %r1)
//           : (!acr.species, !acr.species, !acr.species, f32, f32)
//          -> (!acr.species, !acr.species, !acr.species)
//===----------------------------------------------------------------------===//

def ACR_CRNGenericOp : ACR_Op<"crn_generic"> {
  let summary = "Compact CRN using a stoichiometric matrix attribute";
  let description = [{
    `acr.crn_generic` encodes an entire CRN via a 2-D `stoich_matrix`
    tensor attribute (reactions × species). It takes species concentrations
    and rate constants as runtime arguments (either as tensors or individual
    scalars), and returns the updated species concentrations.

    This operation is the target of canonicalization from `acr.crn` /
    `acr.reaction_instance` chains.
  }];

  let arguments = (ins Variadic<AnyType>:$inputs);
  let results = (outs Variadic<AnyType>:$results);

  // Format: attr-dict first (contains stoich_matrix), then operands and types.
  // Printed as:
  //   {stoich_matrix = dense<...> : tensor<2x5xi32>}
  //       (%species_vec, %rate_vec) :
  //       (tensor<5x!acr.species>, tensor<2xf32>) -> tensor<5x!acr.species>
  let assemblyFormat = [{
    attr-dict `(` $inputs `)` `:` functional-type($inputs, results)
  }];

  let hasVerifier = 1;
}

//===----------------------------------------------------------------------===//
// acr.return
//
// Terminator for acr.reaction and acr.crn region bodies.
//
// Example (scalar):
//   acr.return %h2_o, %o2_o : !acr.species, !acr.species
//
// Example (tensor):
//   acr.return %sv_out : tensor<3x!acr.species>
//===----------------------------------------------------------------------===//

def ACR_ReturnOp : ACR_Op<"return", [Pure, ReturnLike, Terminator]> {
  let summary = "Return from an acr.reaction or acr.crn body";
  let description = [{
    `acr.return` terminates the body of an `acr.reaction` or `acr.crn`
    operation. The types of its operands must match the result types of the
    enclosing function-like op.
  }];

  let arguments = (ins Variadic<AnyType>:$operands);
  let assemblyFormat = "($operands^ `:` type($operands))? attr-dict";
}

#endif // ACR_OPS
