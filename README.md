# ChemComp - Chemical Compiler Tool

ChemComp is set to implement a compilation flow that converts mathematical operations to chemical/biological operations of well-known reaction behavior, providing translation of relevant classes of ODEs.

# Build instructions

ChemComp is built out-of-tree against a specific MLIR revision. Use this commit for llvm-project:
26a1d6601d727a96f4301d0d8647b5a42760ae0c (tag: llvmorg-18.1.2)

ChemComp is compiled and tested against that exact MLIR version. Other revisions may not build or may behave differently.

This repository is published in isolation. The steps below are self-contained.

## Build LLVM/MLIR

```sh
LLVM_SRC_DIR=/path/to/llvm-project
LLVM_BUILD_DIR=/path/to/llvm-build
LLVM_INSTALL_DIR=/path/to/llvm-install

cmake -G Ninja -S "$LLVM_SRC_DIR/llvm" -B "$LLVM_BUILD_DIR" \
    -DCMAKE_INSTALL_PREFIX="$LLVM_INSTALL_DIR" \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DLLVM_INSTALL_UTILS=ON \
    -DLLVM_ENABLE_PROJECTS=mlir \
    -DLLVM_TARGETS_TO_BUILD=host \
    -DLLVM_ENABLE_ASSERTIONS=ON \
    -DMLIR_ENABLE_BINDINGS_PYTHON=ON

cmake --build "$LLVM_BUILD_DIR" --target mlir-opt mlir-translate mlir-cpu-runner install
```

## Build ChemComp (out-of-tree)

```sh
CHEMCOMP_SRC_DIR=/path/to/chemcomp
CHEMCOMP_BUILD_DIR=/path/to/chemcomp-build

cmake -G Ninja -S "$CHEMCOMP_SRC_DIR" -B "$CHEMCOMP_BUILD_DIR" \
    -DCMAKE_BUILD_TYPE=Debug \
    -DMLIR_DIR="$LLVM_INSTALL_DIR/lib/cmake/mlir" \
    -DLLVM_DIR="$LLVM_INSTALL_DIR/lib/cmake/llvm" \
    -DLLVM_EXTERNAL_LIT="$LLVM_BUILD_DIR/bin/llvm-lit"

cmake --build "$CHEMCOMP_BUILD_DIR" --target standalone-opt
```

To run the test suite:
```sh
cmake --build "$CHEMCOMP_BUILD_DIR" --target check-standalone
```

# License

This project is made available under the Apache License v2.0 with LLVM Exceptions. See the [LICENSE.txt](LICENSE.txt) file for more details.

## Software from third parties included in this project

The ChemComp project contains third party software which is under different license terms. All such code will be identified clearly using at least one of two mechanisms:

It will be in a separate directory tree with its own LICENSE.txt or LICENSE file at the top containing the specific license and restrictions which apply to that software, or
It will contain specific license and restriction terms at the top of every file.

# Citations

If you wish to cite this work or learn more, refer to [CITATIONS.md](CITATIONS.md).

# Disclaimer Notice

This material was prepared as an account of work sponsored by an agency of the United States Government.  Neither the United States Government nor the United States Department of Energy, nor Battelle, nor any of their employees, nor any jurisdiction or organization that has cooperated in the development of these materials, makes any warranty, express or implied, or assumes any legal liability or responsibility for the accuracy, completeness, or usefulness or any information, apparatus, product, software, or process disclosed, or represents that its use would not infringe privately owned rights.

Reference herein to any specific commercial product, process, or service by trade name, trademark, manufacturer, or otherwise does not necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or any agency thereof, or Battelle Memorial Institute. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or any agency thereof.

<div align=center>
<pre style="align-text:center">
PACIFIC NORTHWEST NATIONAL LABORATORY
operated by
BATTELLE
for the
UNITED STATES DEPARTMENT OF ENERGY
under Contract DE-AC05-76RL01830
</pre>
</div>