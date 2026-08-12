# A conversion library rather than a plain one: standalone-opt already links
# everything in the global MLIR_CONVERSION_LIBS property, so registering here
# means the pass reaches the tool without editing its CMakeLists.
add_mlir_conversion_library(MLIRCRN2DB
    CRNDBMatchPass.cpp
    MotifDatabase.cpp
    MotifDatabaseFixture.cpp
    MotifDatabaseJSON.cpp

    ADDITIONAL_HEADER_DIRS
    ${PROJECT_SOURCE_DIR}/include/Standalone/CRN2DB

    DEPENDS
    ChemCompMotifsIncGen
    ChemCompCRN2DBPassIncGen

    LINK_LIBS PUBLIC
    MLIRACR
    MLIRIR
    MLIRPass
    MLIRSupport
)
