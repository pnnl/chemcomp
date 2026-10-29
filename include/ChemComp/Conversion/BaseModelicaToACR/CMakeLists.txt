set(LLVM_TARGET_DEFINITIONS Passes.td)
mlir_tablegen(Passes.h.inc -gen-pass-decls -name BaseModelicaToACR)
add_public_tablegen_target(ChemCompBaseModelicaToACRPassIncGen)
add_mlir_doc(Passes BaseModelicaToACRPasses ./ -gen-pass-doc)
