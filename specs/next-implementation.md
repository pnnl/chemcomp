# Next Implementation Steps for ChemComp

I would like to implement the infrastructure for chemcomp, which will be a compiler that translates ODEs into chemical reactions. This will involve creating a new MLIR dialect for chemical reactions, as well as implementing the necessary passes to perform the translation.

For now we focus on the chemical reaction dialect.

This will be a MLIR based implementation.

## Structure of the project and build instructions

We are on a docker development environment with everything needed to build our compiler.

The project structure is as follows.

```
/workspaces/chemical-compiler/chemcomp
├── CMakeLists.txt
├── LICENSE.txt
├── README.md
├── docs
│   ├── asp-dac2024
│   │   ├── mlir.json
│   │   └── reactions.mlir
│   ├── cr_explanations.md
│   ├── crn_explanation.md
│   ├── mass_action_knetics.md
│   ├── project_goals.md
│   └── types_of_chem_reactions.md
├── include
│   ├── CMakeLists.txt
│   ├── Standalone
│   │   ├── CMakeLists.txt
│   │   └── Dialect
│   │       ├── ACR
│   │       │   ├── ACRDialect.h
│   │       │   ├── ACRDialect.td
│   │       │   ├── ACROps.h
│   │       │   ├── ACROps.td
│   │       │   ├── ACRPasses.h
│   │       │   ├── ACRPasses.td
│   │       │   ├── ACRTypes.h
│   │       │   ├── ACRTypes.td
│   │       │   └── CMakeLists.txt
│   │       ├── CMakeLists.txt
│   │       └── Standalone
│   │           ├── CMakeLists.txt
│   │           ├── StandaloneDialect.h
│   │           ├── StandaloneDialect.td
│   │           ├── StandaloneOps.h
│   │           ├── StandaloneOps.td
│   │           ├── StandalonePasses.h
│   │           ├── StandalonePasses.td
│   │           ├── StandaloneTypes.h
│   │           └── StandaloneTypes.td
│   └── Standalone-c
│       └── Dialects.h
├── lib
│   ├── CAPI
│   │   ├── CMakeLists.txt
│   │   └── Dialects.cpp
│   ├── CMakeLists.txt
│   └── Dialect
│       ├── ACR
│       │   ├── CMakeLists.txt
│       │   ├── IR
│       │   │   ├── ACRDialect.cpp
│       │   │   ├── ACROps.cpp
│       │   │   └── ACRTypes.cpp
│       │   └── Transforms
│       │       └── ACRPasses.cpp
│       ├── CMakeLists.txt
│       └── Standalone
│           ├── CMakeLists.txt
│           ├── IR
│           │   ├── StandaloneDialect.cpp
│           │   ├── StandaloneOps.cpp
│           │   └── StandaloneTypes.cpp
│           └── Transforms
│               └── StandalonePasses.cpp
├── python
│   ├── CMakeLists.txt
│   ├── StandaloneExtension.cpp
│   └── mlir_standalone
│       └── dialects
│           ├── ACROps.td
│           ├── StandaloneOps.td
│           ├── acr.py
│           └── standalone.py
├── specs
│   └── next-implementation.md
├── standalone-opt
│   ├── CMakeLists.txt
│   └── standalone-opt.cpp
├── standalone-plugin
│   ├── CMakeLists.txt
│   └── standalone-plugin.cpp
├── standalone-translate
│   ├── CMakeLists.txt
│   └── standalone-translate.cpp
└── test
    ├── CAPI
    │   ├── CMakeLists.txt
    │   ├── lit.local.cfg
    │   └── standalone-capi-test.c
    ├── CMakeLists.txt
    ├── Dialects
    │   ├── ACR
    │   │   ├── roundtrip.mlir
    │   │   └── types.mlir
    │   ├── ACRN
    │   │   └── roundtrip.mlir
    │   ├── CR
    │   ├── ODE
    │   │   └── roundtrip.mlir
    │   ├── Transform
    │   │   ├── sequence.mlir
    │   │   └── transform-parse.mlir
    │   └── end-to-end
    │       ├── odeint-pendulum.py
    │       └── odeint-simple.py
    ├── Standalone
    │   ├── dummy.mlir
    │   ├── standalone-opt.mlir
    │   ├── standalone-pass-plugin.mlir
    │   ├── standalone-plugin.mlir
    │   └── standalone-translate.mlir
    ├── lit.cfg.py
    ├── lit.site.cfg.py.in
    ├── python
    │   ├── lit.local.cfg
    │   └── smoketest.py
    └── translate
        └── example.mlir
/workspaces/chemical-compiler/external/llvm-project/mlir # source code for llvm-project, NEVER make any changes here
```


With CMAKE build directories in

```
/workspaces/chemical-compiler/builds/
├── chemcomp
│   └── build
└── llvm-project
    ├── build
    └── install
```

DONT MAKE any changes to the llvm-project.

To build the MLIR compiler that parses dialects we can use:

```bash
cmake --build /workspaces/chemical-compiler/builds/chemcomp/build/ --target standalone-opt
```

To run all tests:

```bash
cmake --build /workspaces/chemical-compiler/builds/chemcomp/build/ --target check-standalone
```


# Actions 

## Understanding the features

- READ /workspaces/chemical-compiler/chemcomp/docs/crn_explanation.md to understand how chemical reactions work
- READ /workspaces/chemical-compiler/chemcomp/docs/project_goals.md to understand the overall goals of the project and the structure of the codebase.
- READ /workspaces/chemical-compiler/chemcomp/docs/cr_explanations.md to understand the chemical reaction dialect and its operations.

## MLIR Dialect compilation

### MLIR Concepts

- READ /workspaces/chemical-compiler/chemcomp/specs/MLIR-dialect-design.md to understand the structure of a MLIR dialect implementation and how the different files work together. This will be important for implementing the chemical reaction dialect.

### Examples

For inspiration on how a complex tablegen MLIR dialect are implemented look into:
- READ /workspaces/chemical-compiler/external/llvm-project/mlir/include/mlir/Dialect/Tosa/IR/TosaTypesBase.td
- READ /workspaces/chemical-compiler/external/llvm-project/mlir/include/mlir/Dialect/Tensor/IR/TensorOps.td
- READ /workspaces/chemical-compiler/external/llvm-project/mlir/test/Dialect/Tensor/ops.mlir

### Ongoing tests

These tests are currently being implemented and will fail on check-standalone cmake target. Make sure they match the design and that they pass by the end of the implementation.

- READ /workspaces/chemical-compiler/chemcomp/test/Dialects/ACR/roundtrip.mlir
- READ /workspaces/chemical-compiler/chemcomp/test/Dialects/ACR/types.mlir
- READ /workspaces/chemical-compiler/chemcomp/test/Dialects/ACRN/roundtrip.mlir

## Implementing the chemical reaction dialect

- Provide a plan
- Make sure tests in /workspaces/chemical-compiler/chemcomp/test/Dialects reflect the discussed design
- Implement the chemical reaction dialect in the include/Standalone/Dialect/ACR directory. This will involve updating the necessary TableGen files to define the dialect, its operations, and its types. A fast way to check for compilation errors is to run `cmake --build /workspaces/chemical-compiler/builds/chemcomp/build/ --target MLIRACR` after making changes to the dialect.
- Implement the necessary definitions in the lib/Dialect/ACR directory. This will involve implementing the logic for the dialect, its operations, and its types when need. A fast way to check for compilation errors is to run `cmake --build /workspaces/chemical-compiler/builds/chemcomp/build/ --target MLIRACR` after making changes to the implementation.
- Make sure all tests pass by running `cmake --build /workspaces/chemical-compiler/builds/chemcomp/build/ --target check-standalone`
