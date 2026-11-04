add_mlir_dialect_library(MLIRODE
        IR/ODEDialect.cpp
        IR/ODEOps.cpp
        Utils/ODESystem.cpp
        Transforms/ODEPasses.cpp
        Transforms/DualRail.cpp

        ADDITIONAL_HEADER_DIRS
        ${PROJECT_SOURCE_DIR}/include/ChemComp/Dialect/ODE

        DEPENDS
        MLIRODEOpsIncGen
        MLIRODEPassesIncGen

        LINK_LIBS PUBLIC
        MLIRIR
        MLIRPass
        MLIRSupport
        )
