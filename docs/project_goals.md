# ChemComp Project Goals


## ChemComp - Chemical Compiler Tool

ChemComp is set to implement a compilation flow that converts mathematical operations to chemical/biological operations of well-known reaction behavior, providing translation of relevant classes of ODEs. 


## Folder structure

```
.
├── chemcomp
│   ├── CMakeLists.txt
│   ├── LICENSE.txt
│   ├── README.md
│   ├── include
│   │   ├── CMakeLists.txt
│   │   ├── Standalone
│   │   └── Standalone-c
│   ├── lib
│   │   ├── CAPI
│   │   ├── CMakeLists.txt
│   │   └── Standalone
│   ├── python
│   │   ├── CMakeLists.txt
│   │   ├── StandaloneExtension.cpp
│   │   └── mlir_standalone
│   ├── standalone-opt
│   │   ├── CMakeLists.txt
│   │   └── standalone-opt.cpp
│   ├── standalone-plugin
│   │   ├── CMakeLists.txt
│   │   └── standalone-plugin.cpp
│   ├── standalone-translate
│   │   ├── CMakeLists.txt
│   │   └── standalone-translate.cpp
│   └── test
│       ├── CAPI
│       ├── CMakeLists.txt
│       ├── Dialects
│       ├── Standalone
│       ├── lit.cfg.py
│       ├── lit.site.cfg.py.in
│       ├── python
│       └── translate
├── copilot
│   └── project_goals.md # This file
├── external
│   └── llvm-project
├── builds
│   ├── chemcomp
│   │   └── build
│   └── llvm-project
│       ├── build
│       └── install
└── todo.md
```

## Rationale

ChemComp is an MLIR based compiler that will implement abstractions to compile ODEs into CRNs.

We are currently implementing the following dialects that enable progressive lowering of ODEs to CRNs:

- ode: A dialect that represents ODEs.
- acr: A dialect that represents abstract chemical reactions.
- acrn: A dialect that represents CRNs using vectors and tensors

## Ongoing

- Implementing the crn dialect