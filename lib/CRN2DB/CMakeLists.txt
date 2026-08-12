add_mlir_library(MLIRCRN2DB
    MotifDatabase.cpp

    ADDITIONAL_HEADER_DIRS
    ${PROJECT_SOURCE_DIR}/include/Standalone/CRN2DB

    LINK_LIBS PUBLIC
    MLIRSupport
)
