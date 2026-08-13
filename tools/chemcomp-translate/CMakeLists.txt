set(LLVM_LINK_COMPONENTS
  Support
  )

get_property(dialect_libs GLOBAL PROPERTY MLIR_DIALECT_LIBS)
get_property(translation_libs GLOBAL PROPERTY MLIR_TRANSLATION_LIBS)

add_llvm_executable(chemcomp-translate
  chemcomp-translate.cpp
  )
llvm_update_compile_flags(chemcomp-translate)
target_link_libraries(chemcomp-translate
  PRIVATE
  ${dialect_libs}
  ${translation_libs}
  MLIRACR
  MLIRIR
  MLIRParser
  MLIRPass
  MLIRSPIRVDialect
  MLIRTargetJSON
  MLIRTranslateLib
  MLIRSupport
  )

mlir_check_link_libraries(chemcomp-translate)
