add_subdirectory(CRNToDB)
add_subdirectory(ODEToACR)

if(CHEMCOMP_ENABLE_MARCO)
  add_subdirectory(BaseModelicaToACR)
endif()
