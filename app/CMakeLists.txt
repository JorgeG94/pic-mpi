# Only reached when PMPI_BUILD_STANDALONE=ON; see the top-level CMakeLists. Both
# drivers only use `pic_mpi_lib`, so they build against whichever backend the
# library selected, serial included.
foreach(app main legacy_main)
  add_executable(${app} ${app}.f90)
  target_link_libraries(${app} PRIVATE ${main_lib})
endforeach()
