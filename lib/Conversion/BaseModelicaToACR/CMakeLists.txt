# Built only with CHEMCOMP_ENABLE_MARCO=ON: links the MARCO Base Modelica
# dialect (LGPLv3).
add_mlir_library(MLIRBaseModelicaToACR
    BaseModelicaToACR.cpp
    KnownCRNRules.cpp
    ODEPolynomial.cpp

    ADDITIONAL_HEADER_DIRS
    ${PROJECT_SOURCE_DIR}/include/ChemComp/Conversion/BaseModelicaToACR

    DEPENDS
    ChemCompBaseModelicaToACRPassIncGen

    LINK_LIBS PUBLIC
    MLIRACR
    MLIRIR
    MLIRPass
    MLIRSupport
    marco::MLIRBaseModelica
)
