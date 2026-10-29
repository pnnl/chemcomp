get_property(dialect_libs GLOBAL PROPERTY MLIR_DIALECT_LIBS)
get_property(conversion_libs GLOBAL PROPERTY MLIR_CONVERSION_LIBS)
set(LIBS
        ${dialect_libs}
        ${conversion_libs}
        MLIRArithDialect
        MLIROptLib
        MLIRRegisterAllPasses
        MLIRACR
        )

# The MARCO dialects and the conversions using them are only linked when
# CHEMCOMP_ENABLE_MARCO is ON; otherwise registerMarco() is a no-op.
if(CHEMCOMP_ENABLE_MARCO)
  set(REGISTER_MARCO_SOURCE RegisterMarco.cpp)
  list(APPEND LIBS
        MLIRBaseModelicaToACR
        marco::MLIRBaseModelica
        marco::MLIRBaseModelicaTransforms
        marco::MLIRModeling
        )
else()
  set(REGISTER_MARCO_SOURCE RegisterMarcoStub.cpp)
endif()

set(LLVM_OPTIONAL_SOURCES RegisterMarco.cpp RegisterMarcoStub.cpp)
add_llvm_executable(chemcomp-opt chemcomp-opt.cpp ${REGISTER_MARCO_SOURCE})

llvm_update_compile_flags(chemcomp-opt)
target_link_libraries(chemcomp-opt PRIVATE ${LIBS})

mlir_check_all_link_libraries(chemcomp-opt)
