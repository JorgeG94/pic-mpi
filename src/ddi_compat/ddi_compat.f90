!> DDI-compatible API wrapper
!!
!! Provides exact DDI API for GAMESS compatibility.
!! Routes all calls through pic-mpi's darrays and groups modules.
!!
!! Usage:
!!   use ddi_compat
!!   call ddi_init()
!!   call ddi_create(idim, jdim, handle)
!!   call ddi_put(handle, ilo, ihi, jlo, jhi, buffer)
!!   call ddi_get(handle, ilo, ihi, jlo, jhi, buffer)
!!   call ddi_destroy(handle)
!!   call ddi_finalize()
!!
module ddi_compat
   use iso_c_binding, only: c_loc, c_f_pointer
   use pic_types, only: int32, int64, sp, dp
   use pic_mpi_lib, only: pic_mpi_init, pic_mpi_finalize, comm_world, comm_t, &
                          allreduce
   use mpi_f08, only: MPI_Bcast, MPI_Send, MPI_Recv, MPI_DOUBLE_PRECISION, &
                      MPI_INTEGER, MPI_STATUS_IGNORE

   ! Import darrays functionality
   use darrays, only: darrays_init, darrays_finalize, &
                      darray_create, darray_destroy, darray_distrib, &
                      darray_get, darray_put, darray_acc, &
                      dlb_init, dlb_finalize, dlb_reset, dlb_next

   ! Import groups functionality
   use groups, only: groups_init, groups_finalize, &
                     ddi_group_create, ddi_scope, ddi_ascope, &
                     get_working_comm, groups_get_comm, &
                     DDI_COMM_WORLD, DDI_COMM_GROUP, DDI_COMM_MASTER

   implicit none
   private

   ! Re-export DDI communicator constants
   public :: DDI_COMM_WORLD, DDI_COMM_GROUP, DDI_COMM_MASTER

   ! Initialization and finalization
   public :: ddi_init, ddi_finalize, ddi_memory

   ! Array operations (generic interfaces for dp, i32, i64)
   public :: ddi_create, ddi_destroy, ddi_distrib
   public :: ddi_get, ddi_put, ddi_acc

   ! Communication
   public :: ddi_gsumf, ddi_gsumi, ddi_bcast
   public :: ddi_send, ddi_recv

   ! Synchronization
   public :: ddi_sync, ddi_nproc

   ! Load balancing
   public :: ddi_dlbreset, ddi_dlbnext

   ! Groups and scope
   public :: ddi_group_create, ddi_scope, ddi_ascope

   ! Generic interfaces for array operations
   interface ddi_create
      module procedure ddi_create_dp
      module procedure ddi_create_i32
      module procedure ddi_create_i64
   end interface ddi_create

   interface ddi_get
      module procedure ddi_get_dp
      module procedure ddi_get_i32
      module procedure ddi_get_i64
   end interface ddi_get

   interface ddi_put
      module procedure ddi_put_dp
      module procedure ddi_put_i32
      module procedure ddi_put_i64
   end interface ddi_put

   interface ddi_acc
      module procedure ddi_acc_dp
      module procedure ddi_acc_i32
      module procedure ddi_acc_i64
   end interface ddi_acc

   ! Module state
   type(comm_t), save :: world_comm
   logical, save :: ddi_initialized = .false.

contains

   !> Initialize DDI
   !!
   !! Initializes MPI, groups, and darrays modules.
   subroutine ddi_init()
      call pic_mpi_init()
      world_comm = comm_world()
      call groups_init(world_comm)
      call darrays_init(world_comm)
      call dlb_init(world_comm)
      ddi_initialized = .true.
   end subroutine ddi_init

   !> Finalize DDI
   !!
   !! Cleans up all DDI resources and finalizes MPI.
   subroutine ddi_finalize()
      if (.not. ddi_initialized) return

      call dlb_finalize()
      call darrays_finalize()
      call groups_finalize()
      call world_comm%finalize()
      call pic_mpi_finalize()
      ddi_initialized = .false.
   end subroutine ddi_finalize

   !> Set DDI memory (no-op for MPI-based implementation)
   !!
   !! @param memddi Memory size in megawords (ignored)
   subroutine ddi_memory(memddi)
      integer(int64), intent(in) :: memddi
      ! No-op: MPI handles memory allocation dynamically
   end subroutine ddi_memory

   !> Create a distributed array (double precision)
   !!
   !! @param idim Number of rows
   !! @param jdim Number of columns
   !! @param handle Output: array handle
   !! @param init_val Initial value (determines array type)
   subroutine ddi_create_dp(idim, jdim, handle, init_val)
      integer(int32), intent(in) :: idim, jdim
      integer(int32), intent(out) :: handle
      real(dp), intent(in) :: init_val

      call darray_create(idim, jdim, handle, init_val)
   end subroutine ddi_create_dp

   !> Create a distributed array (32-bit integer)
   subroutine ddi_create_i32(idim, jdim, handle, init_val)
      integer(int32), intent(in) :: idim, jdim
      integer(int32), intent(out) :: handle
      integer(int32), intent(in) :: init_val

      call darray_create(idim, jdim, handle, init_val)
   end subroutine ddi_create_i32

   !> Create a distributed array (64-bit integer)
   subroutine ddi_create_i64(idim, jdim, handle, init_val)
      integer(int32), intent(in) :: idim, jdim
      integer(int32), intent(out) :: handle
      integer(int64), intent(in) :: init_val

      call darray_create(idim, jdim, handle, init_val)
   end subroutine ddi_create_i64

   !> Destroy a distributed array
   !!
   !! @param handle Array handle
   subroutine ddi_destroy(handle)
      integer(int32), intent(in) :: handle
      call darray_destroy(handle)
   end subroutine ddi_destroy

   !> Query distribution for a rank
   !!
   !! Returns the range of rows and columns owned by a rank.
   !! All indices are 1-based (Fortran convention).
   !!
   !! @param handle Array handle
   !! @param rank Rank to query (0-indexed)
   !! @param ilo Output: first row (always 1)
   !! @param ihi Output: last row (always nrows)
   !! @param jlo Output: first column owned by rank
   !! @param jhi Output: last column owned by rank
   subroutine ddi_distrib(handle, rank, ilo, ihi, jlo, jhi)
      integer(int32), intent(in) :: handle, rank
      integer(int32), intent(out) :: ilo, ihi, jlo, jhi
      call darray_distrib(handle, rank, ilo, ihi, jlo, jhi)
   end subroutine ddi_distrib

   !> Get data from a distributed array (double precision)
   !!
   !! Retrieves a patch of data. Indices are 1-based.
   !!
   !! @param handle Array handle
   !! @param ilo First row (1-indexed)
   !! @param ihi Last row (1-indexed)
   !! @param jlo First column (1-indexed)
   !! @param jhi Last column (1-indexed)
   !! @param buffer Output buffer (column-major order)
   subroutine ddi_get_dp(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      real(dp), intent(out) :: buffer(*)
      call darray_get(handle, ilo, ihi, jlo, jhi, buffer)
   end subroutine ddi_get_dp

   !> Get data from a distributed array (32-bit integer)
   subroutine ddi_get_i32(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      integer(int32), intent(out) :: buffer(*)
      call darray_get(handle, ilo, ihi, jlo, jhi, buffer)
   end subroutine ddi_get_i32

   !> Get data from a distributed array (64-bit integer)
   subroutine ddi_get_i64(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      integer(int64), intent(out) :: buffer(*)
      call darray_get(handle, ilo, ihi, jlo, jhi, buffer)
   end subroutine ddi_get_i64

   !> Put data into a distributed array (double precision)
   !!
   !! Stores a patch of data. Indices are 1-based.
   !!
   !! @param handle Array handle
   !! @param ilo First row (1-indexed)
   !! @param ihi Last row (1-indexed)
   !! @param jlo First column (1-indexed)
   !! @param jhi Last column (1-indexed)
   !! @param buffer Input buffer (column-major order)
   subroutine ddi_put_dp(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      real(dp), intent(in) :: buffer(*)
      call darray_put(handle, ilo, ihi, jlo, jhi, buffer)
   end subroutine ddi_put_dp

   !> Put data into a distributed array (32-bit integer)
   subroutine ddi_put_i32(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      integer(int32), intent(in) :: buffer(*)
      call darray_put(handle, ilo, ihi, jlo, jhi, buffer)
   end subroutine ddi_put_i32

   !> Put data into a distributed array (64-bit integer)
   subroutine ddi_put_i64(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      integer(int64), intent(in) :: buffer(*)
      call darray_put(handle, ilo, ihi, jlo, jhi, buffer)
   end subroutine ddi_put_i64

   !> Accumulate data into a distributed array (double precision)
   !!
   !! Adds a patch of data atomically. Indices are 1-based.
   !!
   !! @param handle Array handle
   !! @param ilo First row (1-indexed)
   !! @param ihi Last row (1-indexed)
   !! @param jlo First column (1-indexed)
   !! @param jhi Last column (1-indexed)
   !! @param buffer Input buffer (column-major order)
   subroutine ddi_acc_dp(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      real(dp), intent(in) :: buffer(*)
      call darray_acc(handle, ilo, ihi, jlo, jhi, buffer)
   end subroutine ddi_acc_dp

   !> Accumulate data into a distributed array (32-bit integer)
   subroutine ddi_acc_i32(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      integer(int32), intent(in) :: buffer(*)
      call darray_acc(handle, ilo, ihi, jlo, jhi, buffer)
   end subroutine ddi_acc_i32

   !> Accumulate data into a distributed array (64-bit integer)
   subroutine ddi_acc_i64(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      integer(int64), intent(in) :: buffer(*)
      call darray_acc(handle, ilo, ihi, jlo, jhi, buffer)
   end subroutine ddi_acc_i64

   !> Global sum for double precision array
   !!
   !! In-place sum reduction across all ranks in current scope.
   !!
   !! @param tag Message tag (ignored, for DDI compatibility)
   !! @param buffer In/out: data to sum
   !! @param n Number of elements
   subroutine ddi_gsumf(tag, buffer, n)
      integer(int32), intent(in) :: tag
      real(dp), intent(inout) :: buffer(*)
      integer(int32), intent(in) :: n
      type(comm_t) :: work_comm
      real(dp), allocatable :: temp(:)

      work_comm = get_working_comm()
      allocate (temp(n))
      temp = buffer(1:n)
      call allreduce(work_comm, temp, n)
      buffer(1:n) = temp
      deallocate (temp)
   end subroutine ddi_gsumf

   !> Global sum for integer array
   !!
   !! In-place sum reduction across all ranks in current scope.
   !!
   !! @param tag Message tag (ignored, for DDI compatibility)
   !! @param buffer In/out: data to sum
   !! @param n Number of elements
   subroutine ddi_gsumi(tag, buffer, n)
      integer(int32), intent(in) :: tag
      integer(int32), intent(inout) :: buffer(*)
      integer(int32), intent(in) :: n
      type(comm_t) :: work_comm
      integer(int32), allocatable :: temp(:)

      work_comm = get_working_comm()
      allocate (temp(n))
      temp = buffer(1:n)
      call allreduce(work_comm, temp, n)
      buffer(1:n) = temp
      deallocate (temp)
   end subroutine ddi_gsumi

   !> Broadcast double precision data from root to all ranks
   !!
   !! @param tag Message tag (ignored)
   !! @param dtype Data type: 'F' for double, 'I' for integer
   !! @param buffer Data buffer
   !! @param n Number of elements
   !! @param root Root rank (0-indexed)
   subroutine ddi_bcast(tag, dtype, buffer, n, root)
      integer(int32), intent(in) :: tag
      character(len=1), intent(in) :: dtype
      real(dp), intent(inout), target :: buffer(*)
      integer(int32), intent(in) :: n
      integer(int32), intent(in) :: root
      type(comm_t) :: work_comm
      integer(int32), pointer :: int_buf(:)
      integer(int32) :: ierr

      work_comm = get_working_comm()

      select case (dtype)
      case ('F', 'f')
         call MPI_Bcast(buffer, n, MPI_DOUBLE_PRECISION, root, work_comm%get(), ierr)
      case ('I', 'i')
         ! For integer, reinterpret the buffer
         call c_f_pointer(c_loc(buffer(1)), int_buf, [n])
         call MPI_Bcast(int_buf, n, MPI_INTEGER, root, work_comm%get(), ierr)
      case default
         error stop "ddi_bcast: unsupported dtype"
      end select
   end subroutine ddi_bcast

   !> Point-to-point send
   !!
   !! @param buffer Data to send
   !! @param size Number of bytes
   !! @param to Destination rank
   subroutine ddi_send(buffer, size, to)
      real(dp), intent(in) :: buffer(*)
      integer(int32), intent(in) :: size
      integer(int32), intent(in) :: to
      type(comm_t) :: work_comm
      integer(int32) :: nelems, ierr

      work_comm = get_working_comm()
      ! DDI sends bytes, convert to number of dp elements
      nelems = size/8
      if (nelems > 0) then
         call MPI_Send(buffer, nelems, MPI_DOUBLE_PRECISION, to, 0, work_comm%get(), ierr)
      end if
   end subroutine ddi_send

   !> Point-to-point receive
   !!
   !! @param buffer Output buffer
   !! @param size Number of bytes
   !! @param from Source rank
   subroutine ddi_recv(buffer, size, from)
      real(dp), intent(out) :: buffer(*)
      integer(int32), intent(in) :: size
      integer(int32), intent(in) :: from
      type(comm_t) :: work_comm
      integer(int32) :: nelems, ierr

      work_comm = get_working_comm()
      ! DDI receives bytes, convert to number of dp elements
      nelems = size/8
      if (nelems > 0) then
         call MPI_Recv(buffer, nelems, MPI_DOUBLE_PRECISION, from, 0, work_comm%get(), &
                       MPI_STATUS_IGNORE, ierr)
      end if
   end subroutine ddi_recv

   !> Synchronization barrier
   !!
   !! @param tag Tag (ignored, for DDI compatibility)
   subroutine ddi_sync(tag)
      integer(int32), intent(in) :: tag
      type(comm_t) :: work_comm
      work_comm = get_working_comm()
      call work_comm%barrier()
   end subroutine ddi_sync

   !> Get number of processes and current rank
   !!
   !! @param np Output: number of processes in current scope
   !! @param me Output: current rank (0-indexed) in current scope
   subroutine ddi_nproc(np, me)
      integer(int32), intent(out) :: np, me
      type(comm_t) :: work_comm
      work_comm = get_working_comm()
      np = work_comm%size()
      me = work_comm%rank()
   end subroutine ddi_nproc

   !> Reset dynamic load balancing counter
   subroutine ddi_dlbreset()
      call dlb_reset()
   end subroutine ddi_dlbreset

   !> Get next work unit from dynamic load balancer
   !!
   !! @param counter Output: next unique counter value (0-indexed)
   subroutine ddi_dlbnext(counter)
      integer(int64), intent(out) :: counter
      call dlb_next(counter)
   end subroutine ddi_dlbnext

end module ddi_compat
