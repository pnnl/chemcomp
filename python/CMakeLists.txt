include(AddMLIRPython)

# Specifies that all MLIR packages are co-located under the `mlir_chemcomp`
# top level package (the API has been embedded in a relocatable way).
# TODO: Add an upstream cmake param for this vs having a global here.
add_compile_definitions("MLIR_PYTHON_PACKAGE_PREFIX=mlir_chemcomp.")


################################################################################
# Sources
################################################################################

declare_mlir_python_sources(ChemCompPythonSources)

declare_mlir_dialect_python_bindings(
  ADD_TO_PARENT ChemCompPythonSources
  ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/mlir_chemcomp"
  TD_FILE dialects/ACROps.td
  SOURCES
    dialects/acr.py
  DIALECT_NAME acr)

declare_mlir_python_extension(ChemCompPythonSources.Extension
  MODULE_NAME _chemcompDialects
  ADD_TO_PARENT ChemCompPythonSources
  SOURCES
    ChemCompExtension.cpp
  EMBED_CAPI_LINK_LIBS
    ChemCompCAPI
)

################################################################################
# Common CAPI
################################################################################

add_mlir_python_common_capi_library(ChemCompPythonCAPI
  INSTALL_COMPONENT ChemCompPythonModules
  INSTALL_DESTINATION python_packages/chemcomp/mlir_chemcomp/_mlir_libs
  OUTPUT_DIRECTORY "${MLIR_BINARY_DIR}/python_packages/chemcomp/mlir_chemcomp/_mlir_libs"
  RELATIVE_INSTALL_ROOT "../../../.."
  DECLARED_SOURCES
    ChemCompPythonSources
    # TODO: Remove this in favor of showing fine grained registration once
    # available.
    MLIRPythonExtension.RegisterEverything
    MLIRPythonSources.Core
)

################################################################################
# Instantiation of all Python modules
################################################################################

add_mlir_python_modules(ChemCompPythonModules
  ROOT_PREFIX "${MLIR_BINARY_DIR}/python_packages/chemcomp/mlir_chemcomp"
  INSTALL_PREFIX "python_packages/chemcomp/mlir_chemcomp"
  DECLARED_SOURCES
    ChemCompPythonSources
    # TODO: Remove this in favor of showing fine grained registration once
    # available.
    MLIRPythonExtension.RegisterEverything
    MLIRPythonSources
  COMMON_CAPI_LINK_LIBS
    ChemCompPythonCAPI
  )
