!> Distributed arrays module (DDI-compatible)
!!
!! Provides DDI-style distributed 2D arrays built on MPI-3 RMA.
!! Re-exports all functionality from submodules plus convenience wrappers.
!!
!! Usage:
!!   use darrays
!!   call darrays_init(comm)
!!   call darray_create(nrows, ncols, handle)
!!   call darray_put(handle, 1, nrows, jlo, jhi, data)
!!   call darray_get(handle, 1, nrows, jlo, jhi, data)
!!   call darray_destroy(handle)
!!   call darrays_finalize()
!!
module darrays
   use pic_types, only: int32, int64, sp, dp
   use pic_mpi_lib, only: comm_t, allreduce

   ! Re-export types and constants
   use darrays_types, only: darray_t, DTYPE_DP, DTYPE_SP, DTYPE_I32, DTYPE_I64

   ! Re-export core operations
   use darrays_core, only: darrays_init, darrays_finalize, &
                           darray_create, darray_destroy, darray_distrib, &
                           darray_get, darray_put, darray_acc

   ! Re-export load balancing
   use darrays_dlb, only: dlb_init, dlb_finalize, dlb_reset, dlb_next

   implicit none
   private

   ! Types and constants
   public :: darray_t
   public :: DTYPE_DP, DTYPE_SP, DTYPE_I32, DTYPE_I64

   ! Initialization
   public :: darrays_init, darrays_finalize

   ! Array operations
   public :: darray_create, darray_destroy, darray_distrib
   public :: darray_get, darray_put, darray_acc

   ! Load balancing
   public :: dlb_init, dlb_finalize, dlb_reset, dlb_next

   ! Convenience wrappers (DDI-compatible names)
   public :: darray_gsumf, darray_gsumi, darray_sync, darray_nproc

   ! Module-level communicator for convenience wrappers
   type(comm_t), save :: wrapper_comm
   logical, save :: wrapper_initialized = .false.

   public :: darrays_set_comm

contains

   !> Set communicator for convenience wrappers
   !!
   !! Call after darrays_init to enable gsumf/gsumi/sync/nproc.
   subroutine darrays_set_comm(comm)
      type(comm_t), intent(in) :: comm
      wrapper_comm = comm%duplicate()
      wrapper_initialized = .true.
   end subroutine darrays_set_comm

   !> Global sum for double precision array (DDI_GSUMF equivalent)
   !!
   !! In-place reduction across all ranks.
   subroutine darray_gsumf(buffer, n)
      real(dp), intent(inout) :: buffer(*)
      integer(int32), intent(in) :: n
      real(dp), allocatable :: temp(:)

      allocate (temp(n))
      temp = buffer(1:n)
      call allreduce(wrapper_comm, temp, n)
      buffer(1:n) = temp
      deallocate (temp)
   end subroutine darray_gsumf

   !> Global sum for integer array (DDI_GSUMI equivalent)
   !!
   !! In-place reduction across all ranks.
   subroutine darray_gsumi(buffer, n)
      integer(int32), intent(inout) :: buffer(*)
      integer(int32), intent(in) :: n
      integer(int32), allocatable :: temp(:)

      allocate (temp(n))
      temp = buffer(1:n)
      call allreduce(wrapper_comm, temp, n)
      buffer(1:n) = temp
      deallocate (temp)
   end subroutine darray_gsumi

   !> Synchronization barrier (DDI_SYNC equivalent)
   !!
   !! @param tag Ignored (for DDI compatibility)
   subroutine darray_sync(tag)
      integer(int32), intent(in) :: tag
      call wrapper_comm%barrier()
   end subroutine darray_sync

   !> Get number of processes and current rank (DDI_NPROC equivalent)
   !!
   !! @param np Output: total number of processes
   !! @param me Output: current rank (0-indexed)
   subroutine darray_nproc(np, me)
      integer(int32), intent(out) :: np, me
      np = wrapper_comm%size()
      me = wrapper_comm%rank()
   end subroutine darray_nproc

end module darrays
