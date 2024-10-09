add_mlir_dialect_library(MLIRStandalone
        IR/StandaloneTypes.cpp
        IR/StandaloneDialect.cpp
        IR/StandaloneOps.cpp
        Transforms/StandalonePasses.cpp

        ADDITIONAL_HEADER_DIRS
        ${PROJECT_SOURCE_DIR}/include/Standalone

        DEPENDS
        MLIRStandaloneOpsIncGen
        MLIRStandalonePassesIncGen

        LINK_LIBS PUBLIC
        MLIRIR
        MLIRInferTypeOpInterface
        MLIRFuncDialect
        )
