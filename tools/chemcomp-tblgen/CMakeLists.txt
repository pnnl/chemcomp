set(LLVM_LINK_COMPONENTS
  Support
  TableGen
)

add_tablegen(chemcomp-tblgen CHEMCOMP
  DESTINATION "${CMAKE_INSTALL_BINDIR}"
  EXPORT CHEMCOMP
  chemcomp-tblgen.cpp
  MotifDBGen.cpp
)

set_target_properties(chemcomp-tblgen PROPERTIES FOLDER "Tablegenning")
