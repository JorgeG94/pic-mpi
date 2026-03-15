# Test configuration for pic-mpi Tests are MPI programs and need to be run with
# mpirun

# Find MPI for test execution
find_package(MPI REQUIRED)

# Macro to add an MPI test
macro(add_mpi_test test_name num_procs)
  add_executable(${test_name} ${test_name}.f90)
  target_link_libraries(${test_name} PRIVATE ${main_lib})
  target_include_directories(${test_name} PRIVATE ${CMAKE_BINARY_DIR}/modules)

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
#add_mpi_test(test_communicators 2)
add_mpi_test(test_point_to_point 2)
add_mpi_test(test_darrays 2)
