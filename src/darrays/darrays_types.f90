!> Type definitions for distributed arrays
!!
!! Provides the darray_t type for DDI-style distributed 2D arrays.
!! Columns are distributed across MPI ranks with each rank owning
!! a contiguous block of columns. Supports dp, sp, i32, and i64 data types.
module darrays_types
   use pic_types, only: int32, int64, sp, dp
   use pic_mpi_lib, only: win_t
   implicit none
   private

   public :: darray_t
   public :: DTYPE_DP, DTYPE_SP, DTYPE_I32, DTYPE_I64

   ! Data type identifiers
   integer(int32), parameter :: DTYPE_DP = 1   !! Double precision real
   integer(int32), parameter :: DTYPE_SP = 2   !! Single precision real
   integer(int32), parameter :: DTYPE_I32 = 3  !! 32-bit integer
   integer(int32), parameter :: DTYPE_I64 = 4  !! 64-bit integer

   !> Distributed array descriptor
   !!
   !! Represents a 2D array distributed by columns across MPI ranks.
   !! Each rank owns a contiguous block of columns. Only one data pointer
   !! is active at a time, determined by the dtype field.
   type :: darray_t
      integer(int32) :: handle = -1         !! Unique array handle
      integer(int32) :: dtype = 0           !! Data type (DTYPE_DP, DTYPE_SP, etc.)
      integer(int32) :: nrows = 0           !! Total number of rows
      integer(int32) :: ncols = 0           !! Total number of columns
      integer(int32) :: my_first_col = 0    !! First column owned (0-indexed)
      integer(int32) :: my_ncols = 0        !! Number of columns owned
      integer(int64) :: local_size = 0      !! Size of local data (nrows * my_ncols)
      ! Data pointers - only one is active based on dtype
      real(dp), pointer :: data_dp(:) => null()
      real(sp), pointer :: data_sp(:) => null()
      integer(int32), pointer :: data_i32(:) => null()
      integer(int64), pointer :: data_i64(:) => null()
      type(win_t) :: win                    !! MPI window for RMA access
      logical :: active = .false.           !! Is this array slot in use?
   end type darray_t

end module darrays_types
