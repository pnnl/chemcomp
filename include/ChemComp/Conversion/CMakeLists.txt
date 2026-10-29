add_subdirectory(CRNToDB)

# Conversion pass declarations are centralized here: every conversion's pass
# def lives in Passes.td and shares a single IncGen target.
set(LLVM_TARGET_DEFINITIONS Passes.td)
mlir_tablegen(Passes.h.inc -gen-pass-decls -name Conversion)
add_public_tablegen_target(ChemCompConversionPassIncGen)
add_mlir_doc(Passes ConversionPasses ./ -gen-pass-doc)

if(CHEMCOMP_ENABLE_MARCO)
  add_subdirectory(BaseModelicaToACR)
endif()
