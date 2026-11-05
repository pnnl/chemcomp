add_mlir_conversion_library(MLIRODEToACR
    ODEToACR.cpp
    Pipeline.cpp

    ADDITIONAL_HEADER_DIRS
    ${PROJECT_SOURCE_DIR}/include/ChemComp/Conversion/ODEToACR

    DEPENDS
    ChemCompODEToACRPassIncGen
    ChemCompODERulesIncGen

    LINK_LIBS PUBLIC
    MLIRACR
    MLIRODE
    MLIRIR
    MLIRPass
    MLIRSupport
)
