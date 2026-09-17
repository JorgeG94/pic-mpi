!> Minimal pic-mpi consumer; see the CMakeLists.txt beside it.
program consumer
   use pic_mpi_lib, only: pic_mpi_init, pic_mpi_finalize
   implicit none
   call pic_mpi_init()
   call pic_mpi_finalize()
end program consumer
