add_mlir_translation_library(MLIRTargetJSON
    TranslateToJSON.cpp
    TranslateFromJSON.cpp
    TranslateODEToJSON.cpp
    TranslateRegistration.cpp

    ADDITIONAL_HEADER_DIRS
    ${PROJECT_SOURCE_DIR}/include/ChemComp/Target/JSON

    LINK_LIBS PUBLIC
    MLIRACR
    MLIRODE
    MLIRArithDialect
    MLIRIR
    MLIRSupport
    MLIRTranslateLib
)
