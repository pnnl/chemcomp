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
│   │   ├── ChemComp
│   │   │   ├── Conversion    # CRNToDB, shared Passes.td/Passes.h
│   │   │   ├── Dialect       # ACR
│   │   │   └── Target        # JSON
│   │   └── ChemComp-c
│   ├── lib
│   │   ├── CAPI
│   │   ├── CMakeLists.txt
│   │   ├── Conversion        # mirrors include/ChemComp/Conversion
│   │   ├── Dialect
│   │   └── Target
│   ├── python
│   │   ├── CMakeLists.txt
│   │   ├── ChemCompExtension.cpp
│   │   └── mlir_chemcomp
│   ├── tools
│   │   ├── CMakeLists.txt
│   │   ├── chemcomp-opt
│   │   ├── chemcomp-plugin
│   │   ├── chemcomp-tblgen
│   │   └── chemcomp-translate
│   ├── utils
│   └── test
│       ├── CAPI
│       ├── CMakeLists.txt
│       ├── Conversion
│       ├── Dialects
│       ├── Tools
│       ├── chemcomp-tblgen
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