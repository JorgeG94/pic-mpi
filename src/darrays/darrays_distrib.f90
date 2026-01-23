!> Column distribution helpers for distributed arrays
!!
!! Provides functions to calculate how columns are distributed across
!! MPI ranks and to find owners and local offsets for remote access.
module darrays_distrib
   use pic_types, only: int32, int64, dp
   use mpi_f08, only: MPI_ADDRESS_KIND
   implicit none
   private

   public :: calculate_distribution
   public :: get_owner
   public :: get_local_offset

contains

   !> Calculate column distribution for a rank
   !!
   !! Distributes ncols columns across nranks as evenly as possible.
   !! Extra columns go to lower-numbered ranks.
   !!
   !! @param ncols Total number of columns
   !! @param nranks Number of MPI ranks
   !! @param rank This rank's number (0-indexed)
   !! @param first_col Output: first column owned (0-indexed)
   !! @param my_ncols Output: number of columns owned
   pure subroutine calculate_distribution(ncols, nranks, rank, first_col, my_ncols)
      integer(int32), intent(in) :: ncols, nranks, rank
      integer(int32), intent(out) :: first_col, my_ncols
      integer(int32) :: base_cols, extra

      base_cols = ncols/nranks
      extra = mod(ncols, nranks)

      ! Ranks 0..extra-1 get one extra column
      if (rank < extra) then
         my_ncols = base_cols + 1
         first_col = rank*(base_cols + 1)
      else
         my_ncols = base_cols
         first_col = extra*(base_cols + 1) + (rank - extra)*base_cols
      end if
   end subroutine calculate_distribution

   !> Find which rank owns a given column
   !!
   !! @param ncols Total number of columns
   !! @param nranks Number of MPI ranks
   !! @param col Column index (0-indexed)
   !! @return Rank that owns the column
   pure function get_owner(ncols, nranks, col) result(owner)
      integer(int32), intent(in) :: ncols, nranks, col
      integer(int32) :: owner
      integer(int32) :: base_cols, extra, boundary

      base_cols = ncols/nranks
      extra = mod(ncols, nranks)

      ! Boundary where ranks switch from (base+1) to base columns
      boundary = extra*(base_cols + 1)

      if (col < boundary) then
         owner = col/(base_cols + 1)
      else
         owner = extra + (col - boundary)/base_cols
      end if
   end function get_owner

   !> Get local offset in owner's window for a given (row, col)
   !!
   !! Returns the displacement in the owner's local_data array
   !! for the element at (row, col) in the global array.
   !! Data is stored column-major: offset = row + (local_col * nrows)
   !!
   !! @param nrows Number of rows in array
   !! @param ncols Number of columns in array
   !! @param nranks Number of MPI ranks
   !! @param row Row index (0-indexed)
   !! @param col Column index (0-indexed)
   !! @return Displacement in owner's window (in elements)
   pure function get_local_offset(nrows, ncols, nranks, row, col) result(offset)
      integer(int32), intent(in) :: nrows, ncols, nranks, row, col
      integer(MPI_ADDRESS_KIND) :: offset
      integer(int32) :: owner_first_col, owner_ncols, local_col

      call calculate_distribution(ncols, nranks, get_owner(ncols, nranks, col), &
                                  owner_first_col, owner_ncols)
      local_col = col - owner_first_col
      offset = int(row, MPI_ADDRESS_KIND) + int(local_col, MPI_ADDRESS_KIND)*int(nrows, MPI_ADDRESS_KIND)
   end function get_local_offset

end module darrays_distrib
