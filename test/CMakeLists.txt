llvm_canonicalize_cmake_booleans(
        MLIR_ENABLE_BINDINGS_PYTHON
        CHEMCOMP_ENABLE_MARCO
)

configure_lit_site_cfg(
        ${CMAKE_CURRENT_SOURCE_DIR}/lit.site.cfg.py.in
        ${CMAKE_CURRENT_BINARY_DIR}/lit.site.cfg.py
        MAIN_CONFIG
        ${CMAKE_CURRENT_SOURCE_DIR}/lit.cfg.py
)

set(CHEMCOMP_TEST_DEPENDS
        FileCheck count not
        chemcomp-tblgen
        chemcomp-capi-test
        chemcomp-opt
        chemcomp-translate
        ChemCompPlugin
        )
if(MLIR_ENABLE_BINDINGS_PYTHON)
  list(APPEND CHEMCOMP_TEST_DEPENDS ChemCompPythonModules)
endif()

add_lit_testsuite(check-chemcomp "Running the chemcomp regression tests"
        ${CMAKE_CURRENT_BINARY_DIR}
        DEPENDS ${CHEMCOMP_TEST_DEPENDS}
        )
set_target_properties(check-chemcomp PROPERTIES FOLDER "Tests")

add_lit_testsuites(CHEMCOMP ${CMAKE_CURRENT_SOURCE_DIR} DEPENDS ${CHEMCOMP_TEST_DEPENDS})

add_subdirectory(CAPI)
