add_llvm_library(ChemCompPlugin
        # BUILDTREE_ONLY is only for testing purposes
        MODULE BUILDTREE_ONLY
        chemcomp-plugin.cpp

        DEPENDS
        MLIRACR
        PLUGIN_TOOL
        mlir-opt

        LINK_LIBS
        MLIRACR
        )
target_include_directories(
        ChemCompPlugin
        PRIVATE
        "${STANDALONE_BINARY_DIR}/include"
)
