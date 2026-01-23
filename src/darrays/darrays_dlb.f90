!> Dynamic load balancing for distributed arrays
!!
!! Provides atomic counter-based load balancing using MPI-3 RMA.
!! Rank 0 hosts the counter, all ranks can atomically increment it.
module darrays_dlb
   use pic_types, only: int32, int64, dp
   use pic_mpi_lib, only: comm_t, win_t, win_allocate
   use mpi_f08, only: MPI_ADDRESS_KIND
   implicit none
   private

   public :: dlb_init, dlb_finalize
   public :: dlb_reset, dlb_next

   type(comm_t), save :: dlb_comm
   type(win_t), save :: counter_win
   integer(int64), pointer, save :: counter_ptr(:) => null()
   logical, save :: dlb_initialized = .false.

contains

   !> Initialize dynamic load balancing
   !!
   !! Creates a shared counter window on rank 0.
   !! @param comm MPI communicator to use
   subroutine dlb_init(comm)
      type(comm_t), intent(in) :: comm

      dlb_comm = comm%duplicate()

      ! Allocate counter on all ranks (rank 0's is the master)
      call win_allocate(dlb_comm, 1_int32, counter_ptr, counter_win)
      counter_ptr(1) = 0_int64

      call dlb_comm%barrier()
      dlb_initialized = .true.
   end subroutine dlb_init

   !> Finalize dynamic load balancing
   subroutine dlb_finalize()
      if (.not. dlb_initialized) return

      call dlb_comm%barrier()
      call counter_win%finalize()
      nullify (counter_ptr)
      call dlb_comm%finalize()
      dlb_initialized = .false.
   end subroutine dlb_finalize

   !> Reset the load balancing counter to zero
   !!
   !! Must be called collectively by all ranks.
   subroutine dlb_reset()
      call counter_win%fence()
      if (dlb_comm%leader()) then
         counter_ptr(1) = 0_int64
      end if
      call counter_win%fence()
   end subroutine dlb_reset

   !> Get next work unit atomically
   !!
   !! Atomically increments counter on rank 0 and returns old value.
   !! Each call returns a unique value across all ranks.
   !!
   !! @param counter Output: unique counter value (0-indexed)
   subroutine dlb_next(counter)
      integer(int64), intent(out) :: counter
      integer(int64) :: increment

      increment = 1_int64

      call counter_win%lock(0)
      call counter_win%fetch_and_add_i64(0, 0_MPI_ADDRESS_KIND, increment, counter)
      call counter_win%flush(0)
      call counter_win%unlock(0)
   end subroutine dlb_next

end module darrays_dlb
