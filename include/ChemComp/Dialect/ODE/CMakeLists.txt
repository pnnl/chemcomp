add_mlir_dialect(ODEOps ode)
add_mlir_doc(ODEDialect ODEDialect ODE/ -gen-dialect-doc)
add_mlir_doc(ODEOps ODEOps ODE/ -gen-op-doc)

set(LLVM_TARGET_DEFINITIONS ODEPasses.td)
mlir_tablegen(ODEPasses.h.inc --gen-pass-decls)
add_public_tablegen_target(MLIRODEPassesIncGen)
add_mlir_doc(ODEPasses ODEPasses ODE/ -gen-pass-doc)
