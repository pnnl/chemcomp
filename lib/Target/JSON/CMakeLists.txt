add_mlir_translation_library(MLIRTargetJSON
    TranslateToJSON.cpp
    TranslateFromJSON.cpp
    TranslateRegistration.cpp

    ADDITIONAL_HEADER_DIRS
    ${PROJECT_SOURCE_DIR}/include/Standalone/Target/JSON

    LINK_LIBS PUBLIC
    MLIRACR
    MLIRArithDialect
    MLIRIR
    MLIRSupport
    MLIRTranslateLib
)
