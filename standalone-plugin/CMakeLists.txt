add_llvm_library(StandalonePlugin
        # BUILDTREE_ONLY is only for testing purposes
        MODULE BUILDTREE_ONLY
        standalone-plugin.cpp

        DEPENDS
        MLIRACR
        PLUGIN_TOOL
        mlir-opt

        LINK_LIBS
        MLIRACR
        )
target_include_directories(
        StandalonePlugin
        PRIVATE
        "${STANDALONE_BINARY_DIR}/include"
)
