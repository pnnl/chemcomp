add_subdirectory(CRNToDB)

# Conversion pass declarations are centralized here: every conversion's pass
# def lives in Passes.td and shares a single IncGen target.
#
# The group name is ChemCompConversion, not Conversion: -name generates
# mlir::register<name>Passes(), and upstream MLIR already defines
# mlir::registerConversionPasses() in the same namespace.
set(LLVM_TARGET_DEFINITIONS Passes.td)
mlir_tablegen(Passes.h.inc -gen-pass-decls -name ChemCompConversion)
add_public_tablegen_target(ChemCompConversionPassIncGen)
add_mlir_doc(Passes ConversionPasses ./ -gen-pass-doc)
