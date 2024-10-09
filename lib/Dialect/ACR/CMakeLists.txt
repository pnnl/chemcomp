add_mlir_dialect_library(MLIRACR
        IR/ACRTypes.cpp
        IR/ACRDialect.cpp
        IR/ACROps.cpp
        Transforms/ACRPasses.cpp

        ADDITIONAL_HEADER_DIRS
        ${PROJECT_SOURCE_DIR}/Standalone/Dialect/ACR

        DEPENDS
        MLIRACROpsIncGen
        MLIRACRPassesIncGen

        LINK_LIBS PUBLIC
        MLIRIR
        MLIRInferTypeOpInterface
        MLIRFuncDialect
        )
