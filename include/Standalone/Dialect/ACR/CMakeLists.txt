add_mlir_dialect(ACROps acr)
add_mlir_doc(ACRDialect ACRDialect ACR/ -gen-dialect-doc)
add_mlir_doc(ACROps ACROps ACR/ -gen-op-doc)

set(LLVM_TARGET_DEFINITIONS ACRPasses.td)
mlir_tablegen(ACRPasses.h.inc --gen-pass-decls)
add_public_tablegen_target(MLIRACRPassesIncGen)
