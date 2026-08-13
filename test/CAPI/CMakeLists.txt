# Here we create a single aggregate shared library with the parts of the CAPI
# that we want to bundle together. Then we link a simple C executable
# against it to demonstrate that it does have the fully self contained
# core MLIR library and our own chemcomp dialect.
add_mlir_aggregate(ChemCompCAPITestLib
  SHARED
  EMBED_LIBS
  MLIRCAPIIR
  # TODO: Remove this in favor of showing fine grained dialect registration
  # (once available).
  MLIRCAPIRegisterEverything
  ChemCompCAPI
)

add_llvm_executable(chemcomp-capi-test
  chemcomp-capi-test.c
)
llvm_update_compile_flags(chemcomp-capi-test)
target_link_libraries(chemcomp-capi-test
  PRIVATE ChemCompCAPITestLib)
