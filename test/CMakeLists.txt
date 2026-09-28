#
# The MPI tests need MPI -- to link, and to launch with `mpirun`. A build with
# PIC_ENABLE_MPI=OFF has neither, so they are skipped there and a test of the
# single-rank backend takes their place. That one runs directly: the point of
# the configuration is that there is nothing to launch it with.
#
if(NOT PIC_ENABLE_MPI)
  add_executable(test_serial test_serial.f90)
  target_link_libraries(test_serial PRIVATE ${main_lib})
  add_test(NAME test_serial COMMAND test_serial)
  return()
endif()

# Test configuration for pic-mpi Tests are MPI programs and need to be run with
# mpirun

# Find MPI for test execution
find_package(MPI REQUIRED)

# Macro to add an MPI test
macro(add_mpi_test test_name num_procs)
  # A test that has to reach past the pic-mpi API to the MPI underneath needs to
  # know which backend is active, so it is a .F90 and gets the same USE_LEGACY
  # the library does. Plain .f90 tests are unaffected.
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.F90)
    add_executable(${test_name} ${test_name}.F90)
  else()
    add_executable(${test_name} ${test_name}.f90)
  endif()
  target_link_libraries(${test_name} PRIVATE ${main_lib})
  if(PIC_USE_LEGACY_MPI)
    target_compile_definitions(${test_name} PRIVATE USE_LEGACY)
  endif()

  add_test(
    NAME ${test_name}
    COMMAND
      ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${num_procs}
      ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${test_name}> ${MPIEXEC_POSTFLAGS})

  # Allow oversubscription for CI environments
  set_tests_properties(
    ${test_name} PROPERTIES ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
endmacro()

# Add tests with their required number of MPI processes
add_mpi_test(test_rma_win_allocate 2)
add_mpi_test(test_rma_get_put 2)
add_mpi_test(test_rma_nonblocking 2)
add_mpi_test(test_rma_accumulate 2)
add_mpi_test(test_collectives 2)
# add_mpi_test(test_communicators 2)
add_mpi_test(test_point_to_point 2)

# Four ranks: the split-into-halves check needs two per half to be meaningful.
add_mpi_test(test_comm_from_handle 4)
