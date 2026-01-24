if(PIC_USE_LEGACY_MPI)
  message(STATUS "USING LEGACY")
  add_subdirectory(mpi)
else()
  add_subdirectory(mpi_f08)
endif()
target_sources(${main_lib} PRIVATE pic_mpi_lib.F90)

# Groups and scope management (must come before darrays)
add_subdirectory(groups)

# Distributed arrays (DDI-compatible)
add_subdirectory(darrays)

# DDI compatibility wrapper (must come after darrays and groups)
add_subdirectory(ddi_compat)
