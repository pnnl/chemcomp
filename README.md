# ChemComp - Chemical Compiler Tool

ChemComp is set to implement a compilation flow that converts mathematical operations to chemical/biological operations of well-known reaction behavior, providing translation of relevant classes of ODEs.

# Legacy build instructions from An out-of-tree MLIR dialect

This is an example of an out-of-tree [MLIR](https://mlir.llvm.org/) dialect along with a standalone `opt`-like tool to operate on that dialect.

## Building - Component Build

This setup assumes that you have built LLVM and MLIR in `$BUILD_DIR` and installed them to `$PREFIX`. To build and launch the tests, run
```sh
mkdir build && cd build
cmake -G Ninja .. -DMLIR_DIR=$PREFIX/lib/cmake/mlir -DLLVM_EXTERNAL_LIT=$BUILD_DIR/bin/llvm-lit
cmake --build . --target check-standalone
```
To build the documentation from the TableGen description of the dialect operations, run
```sh
cmake --build . --target mlir-doc
```
**Note**: Make sure to pass `-DLLVM_INSTALL_UTILS=ON` when building LLVM with CMake in order to install `FileCheck` to the chosen installation prefix.

## Building - Monolithic Build

This setup assumes that you build the project as part of a monolithic LLVM build via the `LLVM_EXTERNAL_PROJECTS` mechanism.
To build LLVM, MLIR, the example and launch the tests run
```sh
mkdir build && cd build
cmake -G Ninja `$LLVM_SRC_DIR/llvm` \
    -DCMAKE_BUILD_TYPE=Release \
    -DLLVM_TARGETS_TO_BUILD=host
    -DLLVM_ENABLE_PROJECTS=mlir \
    -DLLVM_EXTERNAL_PROJECTS=standalone-dialect -DLLVM_EXTERNAL_STANDALONE_DIALECT_SOURCE_DIR=../
cmake --build . --target check-standalone
```
Here, `$LLVM_SRC_DIR` needs to point to the root of the monorepo.

# License

This project is made available under the Apache License v2.0 with LLVM Exceptions. See the LICENSE.txt file for more details.

## Software from third parties included in this project

The ChemComp project contains third party software which is under different license terms. All such code will be identified clearly using at least one of two mechanisms:

It will be in a separate directory tree with its own LICENSE.txt or LICENSE file at the top containing the specific license and restrictions which apply to that software, or
It will contain specific license and restriction terms at the top of every file.