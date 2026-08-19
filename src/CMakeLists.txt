#
# One backend directory, chosen here. The modules inside present the same public
# interface, so nothing downstream changes -- `pic_mpi_lib` selects between them
# with the preprocessor and a consumer keeps writing `use pic_mpi_lib`.
#
if(NOT PIC_ENABLE_MPI)
  add_subdirectory(serial)
elseif(PIC_USE_LEGACY_MPI)
  add_subdirectory(mpi)
else()
  add_subdirectory(mpi_f08)
endif()

target_sources(${main_lib} PRIVATE pic_mpi_lib.F90)

# Distributed arrays (DDI-compatible)
if(PIC_ENABLE_DARRAYS)
  add_subdirectory(darrays)
endif()
