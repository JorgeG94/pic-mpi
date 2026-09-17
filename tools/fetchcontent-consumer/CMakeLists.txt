#
# A minimal consumer, used by .github/workflows/consumer.yml.
#
# The point of this project is what it does NOT contain: no include path, no
# knowledge of where pic-mpi keeps its .mod files. It links `pic-mpi` and
# nothing else. If `use pic_mpi_lib` compiles here, the library exported its
# module directory; if it only compiled because the .mod happened to land in a
# directory this project already searched, that is the bug this guards.
#
cmake_minimum_required(VERSION 3.28)
project(pic-mpi-consumer LANGUAGES Fortran)

include(FetchContent)
FetchContent_Declare(pic-mpi SOURCE_DIR "${PIC_MPI_SOURCE_DIR}")
FetchContent_MakeAvailable(pic-mpi)

add_executable(consumer main.f90)
target_link_libraries(consumer PRIVATE pic-mpi)
