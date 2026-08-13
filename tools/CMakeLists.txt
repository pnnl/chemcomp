# chemcomp-tblgen is deliberately absent here: the top-level CMakeLists adds it
# before include/, whose custom commands invoke it. Everything below must come
# after lib/ instead, so the conversion libraries have already registered
# themselves in the global MLIR_CONVERSION_LIBS property that chemcomp-opt reads.
add_subdirectory(chemcomp-opt)
add_subdirectory(chemcomp-plugin)
add_subdirectory(chemcomp-translate)
