# A conversion library rather than a plain one: chemcomp-opt already links
# everything in the global MLIR_CONVERSION_LIBS property, so registering here
# means the pass reaches the tool without editing its CMakeLists.
add_mlir_conversion_library(MLIRCRNToDB
    CRNDBMatchPass.cpp
    MotifDatabase.cpp
    FrontierDP.cpp
    Preprocessing.cpp
    Scoring.cpp
    MotifDatabaseFixture.cpp
    MotifDatabaseJSON.cpp

    ADDITIONAL_HEADER_DIRS
    ${PROJECT_SOURCE_DIR}/include/Standalone/Conversion/CRNToDB

    DEPENDS
    ChemCompMotifsIncGen
    ChemCompConversionPassIncGen

    LINK_LIBS PUBLIC
    MLIRACR
    MLIRIR
    MLIRPass
    MLIRSupport
)
