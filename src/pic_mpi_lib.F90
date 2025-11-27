module pic_mpi_lib
#ifndef USE_LEGACY
use pic_mpi
#else
use pic_mpi_f08
#endif
implicit none


!integer, parameter, public :: MAX_PROCESSOR_NAME=256

end module pic_mpi_lib
