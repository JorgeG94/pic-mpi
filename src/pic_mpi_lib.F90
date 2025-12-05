!> Top-level MPI wrapper module with preprocessor selection
!!
!! This module provides a unified interface to MPI by selecting between
!! the modern mpi_f08 interface (pic_mpi_f08) and the legacy MPI interface
!! (pic_mpi) at compile time using the USE_LEGACY preprocessor flag.
!!
!! Usage:
!!   - Without USE_LEGACY: Uses modern mpi_f08 bindings (recommended)
!!   - With USE_LEGACY: Uses traditional integer-based MPI bindings
!!
!! Both implementations provide the same object-oriented API for consistency.
!!
!! @author Jorge Luis Galvez Vallejo
!! @date 2025
module pic_mpi_lib
#ifdef USE_LEGACY
   use pic_mpi !! Legacy MPI interface (integer handles)
#else
   use pic_mpi_f08 !! Modern MPI interface (mpi_f08)
#endif
   implicit none

end module pic_mpi_lib
