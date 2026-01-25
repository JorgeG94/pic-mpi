!> DDI implementation module
!!
!! Internal implementation of DDI-compatible API for GAMESS.
!! Routes all calls through pic-mpi's darrays and groups modules.
!!
!! Note: For GAMESS linkage, use the external wrappers in ddi_api.f90.
!! This module provides the actual implementation.
!!
module ddi_impl
   use iso_c_binding, only: c_loc, c_f_pointer
   use pic_types, only: int32, int64, sp, dp
   use pic_mpi_lib, only: pic_mpi_init, pic_mpi_finalize, comm_world, comm_t, &
                          allreduce, request_t, waitall
   use mpi_f08, only: MPI_Bcast, MPI_Send, MPI_Recv, MPI_Isend, MPI_Irecv, &
                      MPI_Wait, MPI_DOUBLE_PRECISION, MPI_INTEGER, &
                      MPI_INTEGER8, &
                      MPI_STATUS_IGNORE, MPI_ANY_SOURCE, MPI_Request

   ! Import darrays functionality
   use darrays, only: darrays_init, darrays_finalize, darrays_sync_all, &
                      darray_create, darray_destroy, darray_distrib, &
                      darray_get, darray_put, darray_acc, &
                      dlb_init, dlb_finalize, dlb_reset, dlb_next

   ! Import groups functionality
   use groups, only: groups_init, groups_finalize, &
                     groups_ddi_group_create => ddi_group_create, &
                     groups_ddi_scope => ddi_scope, &
                     groups_ddi_ascope => ddi_ascope, &
                     get_working_comm, groups_get_comm, &
                     groups_register_comm, groups_unregister_comm, &
                     groups_get_ngroups, groups_get_group_id, &
                     DDI_COMM_WORLD, DDI_COMM_GROUP, DDI_COMM_MASTER

   implicit none
   private

   ! Re-export DDI communicator constants
   public :: DDI_COMM_WORLD, DDI_COMM_GROUP, DDI_COMM_MASTER

   ! Initialization and finalization
   public :: ddi_init_impl, ddi_finalize_impl, ddi_memory_impl

   ! Array operations (generic interfaces for dp, i32, i64)
   public :: ddi_create_dp_impl, ddi_create_i32_impl, ddi_create_i64_impl
   public :: ddi_destroy_impl, ddi_distrib_impl
   public :: ddi_get_dp_impl, ddi_get_i32_impl, ddi_get_i64_impl
   public :: ddi_put_dp_impl, ddi_put_i32_impl, ddi_put_i64_impl
   public :: ddi_acc_dp_impl, ddi_acc_i32_impl, ddi_acc_i64_impl
   public :: ddi_zero_impl

   ! Communication
   public :: ddi_gsumf_impl, ddi_gsumi_impl, ddi_bcast_impl
   public :: ddi_send_impl, ddi_recv_impl, ddi_recvany_impl
   public :: ddi_isend_impl, ddi_irecv_impl, ddi_wait_impl

   ! Synchronization
   public :: ddi_sync_impl, ddi_nproc_impl

   ! Node/Group queries
   public :: ddi_nnode_impl, ddi_ngroup_impl

   ! Load balancing
   public :: ddi_dlbreset_impl, ddi_dlbnext_impl
   public :: ddi_gdlbreset_impl, ddi_gdlbnext_impl

   ! Groups and scope
   public :: ddi_group_create_impl, ddi_scope_impl, ddi_ascope_impl

   ! Communicator-specific operations
   public :: ddi_comm_create_impl, ddi_comm_destroy_impl
   public :: ddi_nproc_comm_impl, ddi_sync_comm_impl
   public :: ddi_gsumf_comm_impl, ddi_gsumi_comm_impl, ddi_bcast_comm_impl

   ! Array-wide operations
   public :: ddi_arr_fill_impl, ddi_arr_scale_impl
   public :: ddi_arr_min_impl, ddi_arr_max_impl
   public :: ddi_arr_dot_impl, ddi_arr_add_impl, ddi_arr_acc_impl

   ! Combined operations
   public :: ddi_getacc_impl

   ! Custom group creation
   public :: ddi_group_create_custom_impl

   ! Output control
   public :: ddi_output_impl

   ! Timer functions
   public :: ddi_timer_reset_impl, ddi_timer_output_impl

   ! SMP (shared memory) operations
   public :: ddi_smp_create_impl, ddi_smp_destroy_impl
   public :: ddi_smp_nproc_impl, ddi_smp_offset_impl
   public :: ddi_smp_sync_impl, ddi_smp_bcast_impl, ddi_smp_gsumf_impl

   ! Masters operations
   public :: ddi_masters_gsumf_impl, ddi_masters_bcast_impl

   ! Communicator-specific get/put
   public :: ddi_get_comm_impl, ddi_put_comm_impl

   ! Data server queries (no-op for MPI-3 RMA - no data servers)
   public :: ddi_get_dsid_impl

   ! Scatter accumulate
   public :: ddi_scatter_acc_impl

   ! Process-based dynamic load balancing
   public :: ddi_procdlb_create_impl, ddi_procdlb_destroy_impl
   public :: ddi_procdlb_reset_impl, ddi_procdlb_next_impl

   ! Misc
   public :: ddi_ndistrib_impl, ddi_level_impl

   ! No-op stubs for rarely used functions
   public :: ddi_pend_impl, ddi_pbeg_impl, ddi_addr_test_impl
   public :: ddi_gsum_impl, ddi_gdlbreset_inner_impl, ddi_gdlbnext_inner_impl

   ! Module state
   type(comm_t), save :: world_comm
   logical, save :: ddi_initialized = .false.
   logical, save :: we_initialized_mpi = .false.  ! Track if we called MPI_Init
   integer(int32), save :: output_level = 1

   ! Debug flag - set via environment variable DDI_DEBUG=1
   logical, save :: ddi_debug = .false.
   logical, save :: ddi_debug_checked = .false.

   ! Non-blocking request registry
   integer(int32), parameter :: MAX_REQUESTS = 1000
   type(MPI_Request), save :: request_registry(MAX_REQUESTS)
   logical, save :: request_active(MAX_REQUESTS) = .false.

   ! Group-level DLB state
   type :: group_dlb_t
      integer(int64) :: counter
      logical :: initialized
   end type group_dlb_t
   type(group_dlb_t), save :: group_dlb_counters(32)

   ! Timer state
   real(dp), save :: timer_start = 0.0_dp
   real(dp), save :: timer_total = 0.0_dp

   ! SMP (shared memory) array storage
   integer(int32), parameter :: MAX_SMP_ARRAYS = 100
   type :: smp_array_t
      logical :: active = .false.
      integer(int64) :: size            ! Size in bytes
      integer(int32) :: local_size      ! Local size per process
      real(dp), pointer, contiguous :: data(:) => null()
   end type smp_array_t
   type(smp_array_t), save, target :: smp_arrays(MAX_SMP_ARRAYS)
   type(comm_t), save :: smp_comm           ! Intra-node communicator
   type(comm_t), save :: smp_masters_comm   ! SMP masters (one per node)
   logical, save :: smp_initialized = .false.

   ! Process-based DLB state
   integer(int32), parameter :: MAX_PROCDLB = 32
   type :: procdlb_t
      logical :: active = .false.
      integer(int64) :: counter
      integer(int64) :: total_tasks
      integer(int32) :: last_iproc = -1  ! Track which iproc we're currently working on
   end type procdlb_t
   type(procdlb_t), save :: procdlb_state(MAX_PROCDLB)

contains

   !> Safely get working communicator
   !!
   !! Returns the working comm if initialized, otherwise does lazy init
   !! and returns the working comm. (GAMESS compatibility)
   function safe_get_working_comm() result(work_comm)
      type(comm_t) :: work_comm

      ! Lazy init if called before ddi_init (GAMESS compatibility)
      if (.not. ddi_initialized) then
         call ddi_init_impl()
      end if

      work_comm = get_working_comm()
   end function safe_get_working_comm

   !> Check if DDI debug is enabled (without printing)
   function ddi_debug_enabled() result(enabled)
      logical :: enabled
      character(len=32) :: env_val

      if (.not. ddi_debug_checked) then
         call get_environment_variable("DDI_DEBUG", env_val)
         ddi_debug = (len_trim(env_val) > 0)
         ddi_debug_checked = .true.
      end if
      enabled = ddi_debug
   end function ddi_debug_enabled

   !> Print debug message if DDI_DEBUG is set
   subroutine ddi_debug_print(msg)
      use mpi_f08, only: MPI_Initialized, MPI_COMM_WORLD, MPI_Comm_rank
      character(len=*), intent(in) :: msg
      integer(int32) :: me, ierr
      logical :: mpi_is_init

      if (ddi_debug_enabled()) then
         ! Use MPI_COMM_WORLD directly to avoid recursion through safe_get_working_comm
         call MPI_Initialized(mpi_is_init, ierr)
         if (mpi_is_init) then
            call MPI_Comm_rank(MPI_COMM_WORLD, me, ierr)
         else
            me = 0
         end if
         write (*, '(A,I4,A,A)') "[DDI DEBUG rank ", me, "] ", trim(msg)
         flush (6)
      end if
   end subroutine ddi_debug_print

   !> Initialize DDI
   !!
   !! Initializes MPI, groups, and darrays modules.
   subroutine ddi_init_impl()
      use mpi_f08, only: MPI_Initialized, MPI_Init, MPI_COMM_WORLD, MPI_Comm_rank
      integer(int32) :: i, ierr, me
      logical :: mpi_is_init

      ! Guard against double initialization
      if (ddi_initialized) return

      ! Only init MPI if not already initialized (GAMESS may have done it)
      call MPI_Initialized(mpi_is_init, ierr)
      if (.not. mpi_is_init) then
         call MPI_Init(ierr)
         we_initialized_mpi = .true.
      end if

      call MPI_Comm_rank(MPI_COMM_WORLD, me, ierr)
      world_comm = comm_world()

      call ddi_debug_print("ddi_init ENTER")

      call groups_init(world_comm)
      call darrays_init(world_comm)
      call dlb_init(world_comm)

      ! Initialize request registry
      do i = 1, MAX_REQUESTS
         request_active(i) = .false.
      end do

      ! Initialize group DLB state
      do i = 1, 32
         group_dlb_counters(i)%counter = 0
         group_dlb_counters(i)%initialized = .false.
      end do

      ddi_initialized = .true.
      call ddi_debug_print("ddi_init EXIT - initialized=true")
   end subroutine ddi_init_impl

   !> Finalize DDI
   !!
   !! Cleans up all DDI resources. Only finalizes MPI if we initialized it.
   subroutine ddi_finalize_impl()
      call ddi_debug_print("ddi_finalize ENTER")

      if (.not. ddi_initialized) then
         call ddi_debug_print("ddi_finalize: not initialized, skipping")
         return
      end if

      call dlb_finalize()
      call darrays_finalize()
      call groups_finalize()
      call world_comm%finalize()

      ! Only finalize MPI if we initialized it (not if GAMESS did)
      if (we_initialized_mpi) then
         call pic_mpi_finalize()
         we_initialized_mpi = .false.
      end if

      ddi_initialized = .false.
      call ddi_debug_print("ddi_finalize EXIT")
   end subroutine ddi_finalize_impl

   !> Set DDI memory (no-op for MPI-based implementation)
   !!
   !! @param memddi Memory size in megawords (ignored)
   subroutine ddi_memory_impl(memddi)
      integer(int64), intent(in) :: memddi
      ! No-op: MPI handles memory allocation dynamically
   end subroutine ddi_memory_impl

   !> Set output verbosity level
   !!
   !! @param level Output level (0=silent, 1=normal, 2=verbose)
   subroutine ddi_output_impl(level)
      integer(int32), intent(in) :: level
      output_level = level
   end subroutine ddi_output_impl

   !> Create a distributed array (double precision)
   subroutine ddi_create_dp_impl(idim, jdim, handle, init_val)
      integer(int32), intent(in) :: idim, jdim
      integer(int32), intent(out) :: handle
      real(dp), intent(in) :: init_val
      character(len=128) :: msg

      ! Lazy init if called before ddi_init (GAMESS compatibility)
      if (.not. ddi_initialized) call ddi_init_impl()

      if (ddi_debug_enabled()) then
         write (msg, '(A,I0,A,I0)') "ddi_create ENTER idim=", idim, " jdim=", jdim
         call ddi_debug_print(msg)
      end if
      call darray_create(idim, jdim, handle, init_val)
      if (ddi_debug_enabled()) then
         write (msg, '(A,I0)') "ddi_create EXIT handle=", handle
         call ddi_debug_print(msg)
      end if
   end subroutine ddi_create_dp_impl

   !> Create a distributed array (32-bit integer)
   subroutine ddi_create_i32_impl(idim, jdim, handle, init_val)
      integer(int32), intent(in) :: idim, jdim
      integer(int32), intent(out) :: handle
      integer(int32), intent(in) :: init_val

      ! Lazy init if called before ddi_init (GAMESS compatibility)
      if (.not. ddi_initialized) call ddi_init_impl()

      call darray_create(idim, jdim, handle, init_val)
   end subroutine ddi_create_i32_impl

   !> Create a distributed array (64-bit integer)
   subroutine ddi_create_i64_impl(idim, jdim, handle, init_val)
      integer(int32), intent(in) :: idim, jdim
      integer(int32), intent(out) :: handle
      integer(int64), intent(in) :: init_val

      ! Lazy init if called before ddi_init (GAMESS compatibility)
      if (.not. ddi_initialized) call ddi_init_impl()

      call darray_create(idim, jdim, handle, init_val)
   end subroutine ddi_create_i64_impl

   !> Destroy a distributed array
   subroutine ddi_destroy_impl(handle)
      integer(int32), intent(in) :: handle
      character(len=64) :: msg

      if (ddi_debug_enabled()) then
         write (msg, '(A,I0)') "ddi_destroy ENTER handle=", handle
         call ddi_debug_print(msg)
      end if
      call darray_destroy(handle)
      call ddi_debug_print("ddi_destroy EXIT")
   end subroutine ddi_destroy_impl

   !> Query distribution for a rank
   subroutine ddi_distrib_impl(handle, rank, ilo, ihi, jlo, jhi)
      integer(int32), intent(in) :: handle, rank
      integer(int32), intent(out) :: ilo, ihi, jlo, jhi
      call darray_distrib(handle, rank, ilo, ihi, jlo, jhi)
   end subroutine ddi_distrib_impl

   !> Get data from a distributed array (double precision)
   subroutine ddi_get_dp_impl(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      real(dp), intent(out) :: buffer(*)

      call darray_get(handle, ilo, ihi, jlo, jhi, buffer)
   end subroutine ddi_get_dp_impl

   !> Get data from a distributed array (32-bit integer)
   subroutine ddi_get_i32_impl(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      integer(int32), intent(out) :: buffer(*)
      call darray_get(handle, ilo, ihi, jlo, jhi, buffer)
   end subroutine ddi_get_i32_impl

   !> Get data from a distributed array (64-bit integer)
   subroutine ddi_get_i64_impl(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      integer(int64), intent(out) :: buffer(*)
      call darray_get(handle, ilo, ihi, jlo, jhi, buffer)
   end subroutine ddi_get_i64_impl

   !> Put data into a distributed array (double precision)
   subroutine ddi_put_dp_impl(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      real(dp), intent(in) :: buffer(*)

      call darray_put(handle, ilo, ihi, jlo, jhi, buffer)
   end subroutine ddi_put_dp_impl

   !> Put data into a distributed array (32-bit integer)
   subroutine ddi_put_i32_impl(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      integer(int32), intent(in) :: buffer(*)
      call darray_put(handle, ilo, ihi, jlo, jhi, buffer)
   end subroutine ddi_put_i32_impl

   !> Put data into a distributed array (64-bit integer)
   subroutine ddi_put_i64_impl(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      integer(int64), intent(in) :: buffer(*)
      call darray_put(handle, ilo, ihi, jlo, jhi, buffer)
   end subroutine ddi_put_i64_impl

   !> Accumulate data into a distributed array (double precision)
   subroutine ddi_acc_dp_impl(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      real(dp), intent(in) :: buffer(*)

      call darray_acc(handle, ilo, ihi, jlo, jhi, buffer)
   end subroutine ddi_acc_dp_impl

   !> Accumulate data into a distributed array (32-bit integer)
   subroutine ddi_acc_i32_impl(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      integer(int32), intent(in) :: buffer(*)
      call darray_acc(handle, ilo, ihi, jlo, jhi, buffer)
   end subroutine ddi_acc_i32_impl

   !> Accumulate data into a distributed array (64-bit integer)
   subroutine ddi_acc_i64_impl(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      integer(int64), intent(in) :: buffer(*)
      call darray_acc(handle, ilo, ihi, jlo, jhi, buffer)
   end subroutine ddi_acc_i64_impl

   !> Zero an entire distributed array
   !!
   !! @param handle Array handle
   subroutine ddi_zero_impl(handle)
      integer(int32), intent(in) :: handle
      integer(int32) :: ilo, ihi, jlo, jhi, nrows_patch
      real(dp), allocatable :: zeros(:)
      type(comm_t) :: work_comm
      integer(int32) :: me

      work_comm = safe_get_working_comm()
      me = work_comm%rank()

      ! Get my local distribution
      call darray_distrib(handle, me, ilo, ihi, jlo, jhi)

      ! Only proceed if we own columns
      if (jhi >= jlo) then
         nrows_patch = ihi - ilo + 1
         allocate (zeros(nrows_patch*(jhi - jlo + 1)))
         zeros = 0.0_dp
         call darray_put(handle, ilo, ihi, jlo, jhi, zeros)
         deallocate (zeros)
      end if

      call work_comm%barrier()
   end subroutine ddi_zero_impl

   !> Global sum for double precision array
   subroutine ddi_gsumf_impl(tag, buffer, n)
      integer(int32), intent(in) :: tag
      real(dp), intent(inout) :: buffer(*)
      integer(int32), intent(in) :: n
      type(comm_t) :: work_comm
      real(dp), allocatable :: temp(:)

      work_comm = safe_get_working_comm()
      allocate (temp(n))
      temp = buffer(1:n)
      call allreduce(work_comm, temp, n)
      buffer(1:n) = temp
      deallocate (temp)
   end subroutine ddi_gsumf_impl

   !> Global sum for integer array
   subroutine ddi_gsumi_impl(tag, buffer, n)
      use mpi_f08, only: MPI_Allreduce, MPI_SUM, MPI_IN_PLACE
      integer(int32), intent(in) :: tag
      integer(int64), intent(inout) :: buffer(*)  ! int64 for GAMESS -i8 mode
      integer(int32), intent(in) :: n
      type(comm_t) :: work_comm
      integer(int32) :: ierr

      work_comm = safe_get_working_comm()
      ! Use MPI_INTEGER8 for GAMESS -i8 mode (64-bit integers)
      call MPI_Allreduce(MPI_IN_PLACE, buffer, n, MPI_INTEGER8, MPI_SUM, work_comm%get(), ierr)
   end subroutine ddi_gsumi_impl

   !> Broadcast data from root to all ranks
   subroutine ddi_bcast_impl(tag, dtype, buffer, n, root)
      integer(int32), intent(in) :: tag
      character(len=1), intent(in) :: dtype
      real(dp), intent(inout), target :: buffer(*)
      integer(int32), intent(in) :: n
      integer(int32), intent(in) :: root
      type(comm_t) :: work_comm
      integer(int64), pointer :: int_buf(:)  ! int64 for GAMESS -i8 mode
      integer(int32) :: ierr

      work_comm = safe_get_working_comm()

      select case (dtype)
      case ('F', 'f')
         call MPI_Bcast(buffer, n, MPI_DOUBLE_PRECISION, root, work_comm%get(), ierr)
      case ('I', 'i')
         ! For integer, reinterpret the buffer
         ! GAMESS with -i8 uses 64-bit integers, so use MPI_INTEGER8
         call c_f_pointer(c_loc(buffer(1)), int_buf, [n])
         call MPI_Bcast(int_buf, n, MPI_INTEGER8, root, work_comm%get(), ierr)
      case default
         error stop "ddi_bcast: unsupported dtype"
      end select

      call ddi_debug_print("ddi_bcast EXIT")
   end subroutine ddi_bcast_impl

   !> Point-to-point send
   subroutine ddi_send_impl(buffer, size, to)
      real(dp), intent(in) :: buffer(*)
      integer(int32), intent(in) :: size
      integer(int32), intent(in) :: to
      type(comm_t) :: work_comm
      integer(int32) :: nelems, ierr

      work_comm = safe_get_working_comm()
      ! DDI sends bytes, convert to number of dp elements
      nelems = size/8
      if (nelems > 0) then
         call MPI_Send(buffer, nelems, MPI_DOUBLE_PRECISION, to, 0, work_comm%get(), ierr)
      end if
   end subroutine ddi_send_impl

   !> Point-to-point receive
   subroutine ddi_recv_impl(buffer, size, from)
      real(dp), intent(out) :: buffer(*)
      integer(int32), intent(in) :: size
      integer(int32), intent(in) :: from
      type(comm_t) :: work_comm
      integer(int32) :: nelems, ierr

      work_comm = safe_get_working_comm()
      ! DDI receives bytes, convert to number of dp elements
      nelems = size/8
      if (nelems > 0) then
         call MPI_Recv(buffer, nelems, MPI_DOUBLE_PRECISION, from, 0, work_comm%get(), &
                       MPI_STATUS_IGNORE, ierr)
      end if
   end subroutine ddi_recv_impl

   !> Receive from any source
   subroutine ddi_recvany_impl(buffer, size, from)
      use mpi_f08, only: MPI_Status
      real(dp), intent(out) :: buffer(*)
      integer(int32), intent(in) :: size
      integer(int32), intent(out) :: from
      type(comm_t) :: work_comm
      integer(int32) :: nelems, ierr
      type(MPI_Status) :: status

      work_comm = safe_get_working_comm()
      nelems = size/8
      if (nelems > 0) then
         call MPI_Recv(buffer, nelems, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, 0, &
                       work_comm%get(), status, ierr)
         from = status%MPI_SOURCE
      else
         from = -1
      end if
   end subroutine ddi_recvany_impl

   !> Non-blocking send
   subroutine ddi_isend_impl(buffer, size, to, req_handle)
      real(dp), intent(in) :: buffer(*)
      integer(int32), intent(in) :: size
      integer(int32), intent(in) :: to
      integer(int32), intent(out) :: req_handle
      type(comm_t) :: work_comm
      integer(int32) :: nelems, ierr, i

      work_comm = safe_get_working_comm()
      nelems = size/8

      ! Find free request slot
      req_handle = -1
      do i = 1, MAX_REQUESTS
         if (.not. request_active(i)) then
            req_handle = i
            exit
         end if
      end do

      if (req_handle < 0) then
         error stop "ddi_isend: no free request slots"
      end if

      if (nelems > 0) then
         call MPI_Isend(buffer, nelems, MPI_DOUBLE_PRECISION, to, 0, &
                        work_comm%get(), request_registry(req_handle), ierr)
         request_active(req_handle) = .true.
      end if
   end subroutine ddi_isend_impl

   !> Non-blocking receive
   subroutine ddi_irecv_impl(buffer, size, from, req_handle)
      real(dp), intent(out) :: buffer(*)
      integer(int32), intent(in) :: size
      integer(int32), intent(in) :: from
      integer(int32), intent(out) :: req_handle
      type(comm_t) :: work_comm
      integer(int32) :: nelems, ierr, i

      work_comm = safe_get_working_comm()
      nelems = size/8

      ! Find free request slot
      req_handle = -1
      do i = 1, MAX_REQUESTS
         if (.not. request_active(i)) then
            req_handle = i
            exit
         end if
      end do

      if (req_handle < 0) then
         error stop "ddi_irecv: no free request slots"
      end if

      if (nelems > 0) then
         call MPI_Irecv(buffer, nelems, MPI_DOUBLE_PRECISION, from, 0, &
                        work_comm%get(), request_registry(req_handle), ierr)
         request_active(req_handle) = .true.
      end if
   end subroutine ddi_irecv_impl

   !> Wait for non-blocking operation to complete
   subroutine ddi_wait_impl(req_handle)
      integer(int32), intent(in) :: req_handle
      integer(int32) :: ierr

      if (req_handle < 1 .or. req_handle > MAX_REQUESTS) return
      if (.not. request_active(req_handle)) return

      call MPI_Wait(request_registry(req_handle), MPI_STATUS_IGNORE, ierr)
      request_active(req_handle) = .false.
   end subroutine ddi_wait_impl

   !> Synchronization barrier
   !!
   !! Also syncs all array windows to ensure RMA visibility.
   subroutine ddi_sync_impl(tag)
      integer(int32), intent(in) :: tag
      type(comm_t) :: work_comm
      character(len=64) :: msg

      if (ddi_debug_enabled()) then
         write (msg, '(A,I0)') "ddi_sync ENTER tag=", tag
         call ddi_debug_print(msg)
      end if

      ! Sync all array windows for memory consistency
      call darrays_sync_all()

      work_comm = safe_get_working_comm()
      call work_comm%barrier()
      call ddi_debug_print("ddi_sync EXIT")
   end subroutine ddi_sync_impl

   !> Get number of processes and current rank
   subroutine ddi_nproc_impl(np, me)
      use mpi_f08, only: MPI_Initialized, MPI_COMM_WORLD, MPI_Comm_size, MPI_Comm_rank
      integer(int32), intent(out) :: np, me
      integer(int32) :: ierr
      logical :: mpi_is_init

      ! This may be called before collective DDI init, so avoid triggering
      ! collective operations. Just query MPI_COMM_WORLD directly.
      call MPI_Initialized(mpi_is_init, ierr)
      if (.not. mpi_is_init) then
         np = 1
         me = 0
         return
      end if

      call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD, me, ierr)
   end subroutine ddi_nproc_impl

   !> Get number of nodes and node rank
   !!
   !! Uses MPI shared memory splitting to determine node topology.
   !! @param nnodes Output: total number of nodes
   !! @param my_node Output: this rank's node ID (0-indexed)
   subroutine ddi_nnode_impl(nnodes, my_node)
      use mpi_f08, only: MPI_Initialized, MPI_COMM_WORLD, MPI_Comm_size, MPI_Comm_rank
      integer(int32), intent(out) :: nnodes, my_node
      integer(int32) :: ierr, world_size, world_rank
      logical :: mpi_is_init

      ! This function may be called by only one rank (MASWRK in GAMESS)
      ! so we must NOT do any collective operations here.
      ! Just return simple defaults based on MPI_COMM_WORLD.

      call MPI_Initialized(mpi_is_init, ierr)
      if (.not. mpi_is_init) then
         ! MPI not even initialized - serial run
         nnodes = 1
         my_node = 0
         return
      end if

      ! Get basic info from MPI_COMM_WORLD (no collectives)
      call MPI_Comm_size(MPI_COMM_WORLD, world_size, ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD, world_rank, ierr)

      ! Assume all ranks are on same node (conservative for non-collective call)
      ! A proper implementation would need collective ops to determine node topology
      nnodes = 1
      my_node = 0
   end subroutine ddi_nnode_impl

   !> Get number of groups and current group ID
   subroutine ddi_ngroup_impl(ngroups, my_group)
      integer(int32), intent(out) :: ngroups, my_group

      ngroups = groups_get_ngroups()
      my_group = groups_get_group_id()
   end subroutine ddi_ngroup_impl

   !> Reset dynamic load balancing counter
   subroutine ddi_dlbreset_impl()
      ! Ensure DDI (and MPI) is initialized before DLB operations
      if (.not. ddi_initialized) call ddi_init_impl()
      call dlb_reset()
   end subroutine ddi_dlbreset_impl

   !> Get next work unit from dynamic load balancer
   subroutine ddi_dlbnext_impl(counter)
      integer(int64), intent(out) :: counter
      ! Ensure DDI (and MPI) is initialized before DLB operations
      if (.not. ddi_initialized) call ddi_init_impl()
      call dlb_next(counter)
   end subroutine ddi_dlbnext_impl

   !> Reset group-level dynamic load balancing counter
   !!
   !! Resets the DLB counter for the current group scope.
   subroutine ddi_gdlbreset_impl()
      type(comm_t) :: work_comm
      integer(int32) :: comm_id
      integer(int32), allocatable :: reset_val(:)

      work_comm = safe_get_working_comm()
      comm_id = 1  ! Use slot index based on current scope

      ! Reset counter on rank 0 of the current scope
      allocate (reset_val(1))
      reset_val(1) = 0
      call work_comm%barrier()
      if (work_comm%leader()) then
         group_dlb_counters(comm_id)%counter = 0
      end if
      group_dlb_counters(comm_id)%initialized = .true.
      call work_comm%barrier()
      deallocate (reset_val)
   end subroutine ddi_gdlbreset_impl

   !> Get next work unit from group-level dynamic load balancer
   !!
   !! Returns unique counter values within the current group scope.
   subroutine ddi_gdlbnext_impl(counter)
      use mpi_f08, only: MPI_Bcast, MPI_INTEGER8
      integer(int64), intent(out) :: counter
      type(comm_t) :: work_comm
      integer(int32) :: comm_id, ierr

      work_comm = safe_get_working_comm()
      comm_id = 1

      ! Leader assigns counters and broadcasts the base
      if (work_comm%leader()) then
         counter = group_dlb_counters(comm_id)%counter
         group_dlb_counters(comm_id)%counter = group_dlb_counters(comm_id)%counter + int(work_comm%size(), int64)
      end if

      ! Broadcast base counter from leader
      call MPI_Bcast(counter, 1, MPI_INTEGER8, 0, work_comm%get(), ierr)

      ! Each rank gets its unique counter value
      counter = counter + int(work_comm%rank(), int64)
   end subroutine ddi_gdlbnext_impl

   !> Create DDI-style groups (wrapper)
   subroutine ddi_group_create_impl(ngroups, world_id, group_id, master_id)
      integer(int32), intent(in) :: ngroups
      integer(int32), intent(out) :: world_id, group_id, master_id
      character(len=64) :: msg

      if (ddi_debug_enabled()) then
         write (msg, '(A,I0)') "ddi_group_create ENTER ngroups=", ngroups
         call ddi_debug_print(msg)
      end if
      call groups_ddi_group_create(ngroups, world_id, group_id, master_id)
      if (ddi_debug_enabled()) then
         write (msg, '(A,I0,A,I0,A,I0)') "ddi_group_create EXIT world=", world_id, " group=", group_id, " master=", master_id
         call ddi_debug_print(msg)
      end if
   end subroutine ddi_group_create_impl

   !> Switch scope (synchronous)
   subroutine ddi_scope_impl(comm_id)
      integer(int32), intent(in) :: comm_id
      character(len=64) :: msg

      if (ddi_debug_enabled()) then
         write (msg, '(A,I0)') "ddi_scope ENTER comm_id=", comm_id
         call ddi_debug_print(msg)
      end if
      call groups_ddi_scope(comm_id)
      call ddi_debug_print("ddi_scope EXIT")
   end subroutine ddi_scope_impl

   !> Switch scope (asynchronous)
   subroutine ddi_ascope_impl(comm_id)
      integer(int32), intent(in) :: comm_id
      character(len=64) :: msg

      if (ddi_debug_enabled()) then
         write (msg, '(A,I0)') "ddi_ascope ENTER comm_id=", comm_id
         call ddi_debug_print(msg)
      end if
      call groups_ddi_ascope(comm_id)
      call ddi_debug_print("ddi_ascope EXIT")
   end subroutine ddi_ascope_impl

   !> Reset timing counters
   subroutine ddi_timer_reset_impl()
      use mpi_f08, only: MPI_Wtime
      timer_start = MPI_Wtime()
      timer_total = 0.0_dp
   end subroutine ddi_timer_reset_impl

   !> Output timing statistics
   subroutine ddi_timer_output_impl()
      use mpi_f08, only: MPI_Wtime
      type(comm_t) :: work_comm
      real(dp) :: elapsed

      work_comm = safe_get_working_comm()
      elapsed = MPI_Wtime() - timer_start

      if (work_comm%leader() .and. output_level > 0) then
         print '(A,F12.3,A)', " DDI Timer: ", elapsed, " seconds"
      end if
   end subroutine ddi_timer_output_impl

   ! ========================================================================
   ! Communicator-specific operations
   ! ========================================================================

   !> Create a sub-communicator by splitting
   !!
   !! @param parent_id Parent communicator ID
   !! @param color Ranks with same color go into same communicator
   !! @param new_comm_id Output: new communicator ID
   subroutine ddi_comm_create_impl(parent_id, color, new_comm_id)
      integer(int32), intent(in) :: parent_id, color
      integer(int32), intent(out) :: new_comm_id
      type(comm_t) :: parent_comm, new_comm

      parent_comm = groups_get_comm(parent_id)
      new_comm = parent_comm%split_by(color)
      call groups_register_comm(new_comm, new_comm_id)
   end subroutine ddi_comm_create_impl

   !> Destroy a user-created communicator
   !!
   !! @param comm_id Communicator ID to destroy (must be >= 4)
   subroutine ddi_comm_destroy_impl(comm_id)
      integer(int32), intent(in) :: comm_id
      call groups_unregister_comm(comm_id)
   end subroutine ddi_comm_destroy_impl

   !> Get number of processes in a specific communicator
   !!
   !! @param comm_id Communicator ID
   !! @param np Output: number of processes
   !! @param me Output: current rank (0-indexed)
   subroutine ddi_nproc_comm_impl(comm_id, np, me)
      integer(int32), intent(in) :: comm_id
      integer(int32), intent(out) :: np, me
      type(comm_t) :: comm

      comm = groups_get_comm(comm_id)
      np = comm%size()
      me = comm%rank()
   end subroutine ddi_nproc_comm_impl

   !> Barrier in a specific communicator
   !!
   !! @param comm_id Communicator ID
   subroutine ddi_sync_comm_impl(comm_id)
      integer(int32), intent(in) :: comm_id
      type(comm_t) :: comm

      comm = groups_get_comm(comm_id)
      call comm%barrier()
   end subroutine ddi_sync_comm_impl

   !> Global sum for double precision in a specific communicator
   !!
   !! @param comm_id Communicator ID
   !! @param buffer In/out: data to sum
   !! @param n Number of elements
   subroutine ddi_gsumf_comm_impl(comm_id, buffer, n)
      integer(int32), intent(in) :: comm_id
      real(dp), intent(inout) :: buffer(*)
      integer(int32), intent(in) :: n
      type(comm_t) :: comm
      real(dp), allocatable :: temp(:)

      comm = groups_get_comm(comm_id)
      allocate (temp(n))
      temp = buffer(1:n)
      call allreduce(comm, temp, n)
      buffer(1:n) = temp
      deallocate (temp)
   end subroutine ddi_gsumf_comm_impl

   !> Global sum for integer in a specific communicator
   !!
   !! @param comm_id Communicator ID
   !! @param buffer In/out: data to sum
   !! @param n Number of elements
   subroutine ddi_gsumi_comm_impl(comm_id, buffer, n)
      use mpi_f08, only: MPI_Allreduce, MPI_SUM, MPI_IN_PLACE
      integer(int32), intent(in) :: comm_id
      integer(int64), intent(inout) :: buffer(*)  ! int64 for GAMESS -i8 compatibility
      integer(int32), intent(in) :: n
      type(comm_t) :: comm
      integer(int32) :: ierr

      comm = groups_get_comm(comm_id)
      call MPI_Allreduce(MPI_IN_PLACE, buffer, n, MPI_INTEGER8, MPI_SUM, comm%get(), ierr)
   end subroutine ddi_gsumi_comm_impl

   !> Broadcast in a specific communicator
   !!
   !! @param comm_id Communicator ID
   !! @param dtype Data type: 'F' for double, 'I' for integer
   !! @param buffer Data buffer
   !! @param n Number of elements
   !! @param root Root rank (0-indexed)
   subroutine ddi_bcast_comm_impl(comm_id, dtype, buffer, n, root)
      integer(int32), intent(in) :: comm_id
      character(len=1), intent(in) :: dtype
      real(dp), intent(inout), target :: buffer(*)
      integer(int32), intent(in) :: n
      integer(int32), intent(in) :: root
      type(comm_t) :: comm
      integer(int64), pointer :: int_buf(:)  ! int64 for GAMESS -i8 mode
      integer(int32) :: ierr

      comm = groups_get_comm(comm_id)

      select case (dtype)
      case ('F', 'f')
         call MPI_Bcast(buffer, n, MPI_DOUBLE_PRECISION, root, comm%get(), ierr)
      case ('I', 'i')
         ! GAMESS with -i8 uses 64-bit integers, so use MPI_INTEGER8
         call c_f_pointer(c_loc(buffer(1)), int_buf, [n])
         call MPI_Bcast(int_buf, n, MPI_INTEGER8, root, comm%get(), ierr)
      case default
         error stop "ddi_bcast_comm: unsupported dtype"
      end select
   end subroutine ddi_bcast_comm_impl

   ! ========================================================================
   ! Array-wide operations
   ! ========================================================================

   !> Fill an entire distributed array with a constant value
   !!
   !! @param handle Array handle
   !! @param value Value to fill with
   subroutine ddi_arr_fill_impl(handle, value)
      integer(int32), intent(in) :: handle
      real(dp), intent(in) :: value
      integer(int32) :: ilo, ihi, jlo, jhi, nrows_patch
      real(dp), allocatable :: fill_buf(:)
      type(comm_t) :: work_comm
      integer(int32) :: me

      work_comm = safe_get_working_comm()
      me = work_comm%rank()

      ! Get my local distribution
      call darray_distrib(handle, me, ilo, ihi, jlo, jhi)

      ! Only proceed if we own columns
      if (jhi >= jlo) then
         nrows_patch = ihi - ilo + 1
         allocate (fill_buf(nrows_patch*(jhi - jlo + 1)))
         fill_buf = value
         call darray_put(handle, ilo, ihi, jlo, jhi, fill_buf)
         deallocate (fill_buf)
      end if

      call work_comm%barrier()
   end subroutine ddi_arr_fill_impl

   !> Scale an entire distributed array by a scalar
   !!
   !! @param handle Array handle
   !! @param scale Scale factor
   subroutine ddi_arr_scale_impl(handle, scale)
      integer(int32), intent(in) :: handle
      real(dp), intent(in) :: scale
      integer(int32) :: ilo, ihi, jlo, jhi, nrows_patch, i
      real(dp), allocatable :: buf(:)
      type(comm_t) :: work_comm
      integer(int32) :: me

      work_comm = safe_get_working_comm()
      me = work_comm%rank()

      call darray_distrib(handle, me, ilo, ihi, jlo, jhi)

      if (jhi >= jlo) then
         nrows_patch = ihi - ilo + 1
         allocate (buf(nrows_patch*(jhi - jlo + 1)))
         call darray_get(handle, ilo, ihi, jlo, jhi, buf)
         do i = 1, size(buf)
            buf(i) = buf(i)*scale
         end do
         call darray_put(handle, ilo, ihi, jlo, jhi, buf)
         deallocate (buf)
      end if

      call work_comm%barrier()
   end subroutine ddi_arr_scale_impl

   !> Find minimum element in a distributed array
   !!
   !! @param handle Array handle
   !! @param minval Output: minimum value
   subroutine ddi_arr_min_impl(handle, minval)
      use mpi_f08, only: MPI_MIN, MPI_Allreduce
      integer(int32), intent(in) :: handle
      real(dp), intent(out) :: minval
      integer(int32) :: ilo, ihi, jlo, jhi, nrows_patch, i
      real(dp), allocatable :: buf(:)
      real(dp) :: local_min, global_min
      type(comm_t) :: work_comm
      integer(int32) :: me, ierr

      work_comm = safe_get_working_comm()
      me = work_comm%rank()

      call darray_distrib(handle, me, ilo, ihi, jlo, jhi)

      local_min = huge(local_min)
      if (jhi >= jlo) then
         nrows_patch = ihi - ilo + 1
         allocate (buf(nrows_patch*(jhi - jlo + 1)))
         call darray_get(handle, ilo, ihi, jlo, jhi, buf)
         do i = 1, size(buf)
            if (buf(i) < local_min) local_min = buf(i)
         end do
         deallocate (buf)
      end if

      call MPI_Allreduce(local_min, global_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
                         work_comm%get(), ierr)
      minval = global_min
   end subroutine ddi_arr_min_impl

   !> Find maximum element in a distributed array
   !!
   !! @param handle Array handle
   !! @param maxval Output: maximum value
   subroutine ddi_arr_max_impl(handle, maxval)
      use mpi_f08, only: MPI_MAX, MPI_Allreduce
      integer(int32), intent(in) :: handle
      real(dp), intent(out) :: maxval
      integer(int32) :: ilo, ihi, jlo, jhi, nrows_patch, i
      real(dp), allocatable :: buf(:)
      real(dp) :: local_max, global_max
      type(comm_t) :: work_comm
      integer(int32) :: me, ierr

      work_comm = safe_get_working_comm()
      me = work_comm%rank()

      call darray_distrib(handle, me, ilo, ihi, jlo, jhi)

      local_max = -huge(local_max)
      if (jhi >= jlo) then
         nrows_patch = ihi - ilo + 1
         allocate (buf(nrows_patch*(jhi - jlo + 1)))
         call darray_get(handle, ilo, ihi, jlo, jhi, buf)
         do i = 1, size(buf)
            if (buf(i) > local_max) local_max = buf(i)
         end do
         deallocate (buf)
      end if

      call MPI_Allreduce(local_max, global_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                         work_comm%get(), ierr)
      maxval = global_max
   end subroutine ddi_arr_max_impl

   !> Compute dot product of two distributed arrays
   !!
   !! @param handle_a First array handle
   !! @param handle_b Second array handle
   !! @param dotprod Output: dot product
   subroutine ddi_arr_dot_impl(handle_a, handle_b, dotprod)
      integer(int32), intent(in) :: handle_a, handle_b
      real(dp), intent(out) :: dotprod
      integer(int32) :: ilo, ihi, jlo, jhi, nrows_patch, i
      real(dp), allocatable :: buf_a(:), buf_b(:)
      type(comm_t) :: work_comm
      integer(int32) :: me

      work_comm = safe_get_working_comm()
      me = work_comm%rank()

      call darray_distrib(handle_a, me, ilo, ihi, jlo, jhi)

      dotprod = 0.0_dp
      if (jhi >= jlo) then
         nrows_patch = ihi - ilo + 1
         allocate (buf_a(nrows_patch*(jhi - jlo + 1)))
         allocate (buf_b(nrows_patch*(jhi - jlo + 1)))
         call darray_get(handle_a, ilo, ihi, jlo, jhi, buf_a)
         call darray_get(handle_b, ilo, ihi, jlo, jhi, buf_b)
         do i = 1, size(buf_a)
            dotprod = dotprod + buf_a(i)*buf_b(i)
         end do
         deallocate (buf_a, buf_b)
      end if

      ! In-place allreduce to sum across all ranks
      call allreduce(work_comm, dotprod)
   end subroutine ddi_arr_dot_impl

   !> Add two distributed arrays: C = A + B
   !!
   !! @param handle_a First array handle (A)
   !! @param handle_b Second array handle (B)
   !! @param handle_c Result array handle (C)
   subroutine ddi_arr_add_impl(handle_a, handle_b, handle_c)
      integer(int32), intent(in) :: handle_a, handle_b, handle_c
      integer(int32) :: ilo, ihi, jlo, jhi, nrows_patch, i
      real(dp), allocatable :: buf_a(:), buf_b(:)
      type(comm_t) :: work_comm
      integer(int32) :: me

      work_comm = safe_get_working_comm()
      me = work_comm%rank()

      call darray_distrib(handle_a, me, ilo, ihi, jlo, jhi)

      if (jhi >= jlo) then
         nrows_patch = ihi - ilo + 1
         allocate (buf_a(nrows_patch*(jhi - jlo + 1)))
         allocate (buf_b(nrows_patch*(jhi - jlo + 1)))
         call darray_get(handle_a, ilo, ihi, jlo, jhi, buf_a)
         call darray_get(handle_b, ilo, ihi, jlo, jhi, buf_b)
         do i = 1, size(buf_a)
            buf_a(i) = buf_a(i) + buf_b(i)
         end do
         call darray_put(handle_c, ilo, ihi, jlo, jhi, buf_a)
         deallocate (buf_a, buf_b)
      end if

      call work_comm%barrier()
   end subroutine ddi_arr_add_impl

   !> Accumulate one distributed array into another: B = B + alpha*A
   !!
   !! @param handle_a Source array handle (A)
   !! @param handle_b Target array handle (B)
   !! @param alpha Scale factor
   subroutine ddi_arr_acc_impl(handle_a, handle_b, alpha)
      integer(int32), intent(in) :: handle_a, handle_b
      real(dp), intent(in) :: alpha
      integer(int32) :: ilo, ihi, jlo, jhi, nrows_patch, i
      real(dp), allocatable :: buf_a(:), buf_b(:)
      type(comm_t) :: work_comm
      integer(int32) :: me

      work_comm = safe_get_working_comm()
      me = work_comm%rank()

      call darray_distrib(handle_a, me, ilo, ihi, jlo, jhi)

      if (jhi >= jlo) then
         nrows_patch = ihi - ilo + 1
         allocate (buf_a(nrows_patch*(jhi - jlo + 1)))
         allocate (buf_b(nrows_patch*(jhi - jlo + 1)))
         call darray_get(handle_a, ilo, ihi, jlo, jhi, buf_a)
         call darray_get(handle_b, ilo, ihi, jlo, jhi, buf_b)
         do i = 1, size(buf_a)
            buf_b(i) = buf_b(i) + alpha*buf_a(i)
         end do
         call darray_put(handle_b, ilo, ihi, jlo, jhi, buf_b)
         deallocate (buf_a, buf_b)
      end if

      call work_comm%barrier()
   end subroutine ddi_arr_acc_impl

   !> Get and accumulate combined operation
   !!
   !! Gets a patch from the array and atomically accumulates buffer into it.
   !! Returns the original values before accumulation.
   !!
   !! @param handle Array handle
   !! @param ilo First row (1-indexed)
   !! @param ihi Last row (1-indexed)
   !! @param jlo First column (1-indexed)
   !! @param jhi Last column (1-indexed)
   !! @param buffer In: values to accumulate, Out: original values
   subroutine ddi_getacc_impl(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      real(dp), intent(inout) :: buffer(*)
      integer(int32) :: nrows_patch, ncols_patch, n
      real(dp), allocatable :: original(:)

      nrows_patch = ihi - ilo + 1
      ncols_patch = jhi - jlo + 1
      n = nrows_patch*ncols_patch

      allocate (original(n))

      ! Get original values
      call darray_get(handle, ilo, ihi, jlo, jhi, original)

      ! Accumulate buffer into array
      call darray_acc(handle, ilo, ihi, jlo, jhi, buffer)

      ! Return original values
      buffer(1:n) = original

      deallocate (original)
   end subroutine ddi_getacc_impl

   !> Create groups with custom sizes
   !!
   !! @param ngroups Number of groups to create
   !! @param group_sizes Array of size ngroups with ranks per group
   !! @param world_id Output: world communicator ID
   !! @param group_id Output: group communicator ID
   !! @param master_id Output: masters communicator ID
   subroutine ddi_group_create_custom_impl(ngroups, group_sizes, world_id, group_id, master_id)
      use mpi_f08, only: MPI_UNDEFINED
      integer(int32), intent(in) :: ngroups
      integer(int32), intent(in) :: group_sizes(*)
      integer(int32), intent(out) :: world_id, group_id, master_id
      type(comm_t) :: world_comm, new_group_comm, new_master_comm
      integer(int32) :: world_size, world_rank
      integer(int32) :: my_group, cumulative, i, color

      world_comm = safe_get_working_comm()
      world_size = world_comm%size()
      world_rank = world_comm%rank()

      ! Determine which group this rank belongs to
      cumulative = 0
      my_group = -1
      do i = 1, ngroups
         if (world_rank >= cumulative .and. world_rank < cumulative + group_sizes(i)) then
            my_group = i - 1  ! 0-indexed
            exit
         end if
         cumulative = cumulative + group_sizes(i)
      end do

      if (my_group < 0) then
         error stop "ddi_group_create_custom: rank not assigned to any group"
      end if

      ! Create group communicator
      new_group_comm = world_comm%split_by(my_group)
      call groups_register_comm(new_group_comm, group_id)

      ! Create masters communicator (rank 0 in each group)
      if (new_group_comm%leader()) then
         color = 0
      else
         color = MPI_UNDEFINED
      end if
      new_master_comm = world_comm%split_by(color)
      call groups_register_comm(new_master_comm, master_id)

      world_id = DDI_COMM_WORLD
   end subroutine ddi_group_create_custom_impl

   ! ========================================================================
   ! SMP (Shared Memory) Operations
   ! ========================================================================

   !> Initialize SMP communicator if not already done
   subroutine init_smp_comm()
      type(comm_t) :: work_comm
      integer(int32) :: smp_rank, color

      if (smp_initialized) return

      work_comm = safe_get_working_comm()
      smp_comm = work_comm%split()  ! Split by shared memory

      ! Create masters communicator (contains only SMP rank 0 from each node)
      smp_rank = smp_comm%rank()
      if (smp_rank == 0) then
         color = 0  ! Masters join comm
      else
         color = 1  ! Non-masters get separate comm (which we discard)
      end if
      smp_masters_comm = work_comm%split_by(color)

      smp_initialized = .true.
   end subroutine init_smp_comm

   !> Create an SMP (shared memory) array
   !!
   !! @param nelems Number of double precision elements
   !! @param handle Output: array handle
   subroutine ddi_smp_create_impl(nelems, handle)
      integer(int64), intent(in) :: nelems
      integer(int32), intent(out) :: handle
      integer(int32) :: i

      call init_smp_comm()

      ! Find free slot
      handle = -1
      do i = 1, MAX_SMP_ARRAYS
         if (.not. smp_arrays(i)%active) then
            handle = i
            exit
         end if
      end do

      if (handle < 0) then
         error stop "ddi_smp_create: no free SMP array slots"
      end if

      ! Each process gets a full copy (replicated, not truly shared)
      ! This wastes memory but provides correct semantics
      allocate (smp_arrays(handle)%data(nelems))
      smp_arrays(handle)%data = 0.0_dp
      smp_arrays(handle)%size = nelems
      smp_arrays(handle)%local_size = int(nelems, int32)
      smp_arrays(handle)%active = .true.
   end subroutine ddi_smp_create_impl

   !> Destroy an SMP array
   !!
   !! @param handle Array handle
   subroutine ddi_smp_destroy_impl(handle)
      integer(int32), intent(in) :: handle

      if (handle < 1 .or. handle > MAX_SMP_ARRAYS) return
      if (.not. smp_arrays(handle)%active) return

      if (associated(smp_arrays(handle)%data)) then
         deallocate (smp_arrays(handle)%data)
      end if
      smp_arrays(handle)%active = .false.
   end subroutine ddi_smp_destroy_impl

   !> Get SMP process count and rank
   !!
   !! @param np Output: number of SMP processes
   !! @param me Output: SMP rank (0-indexed)
   subroutine ddi_smp_nproc_impl(np, me)
      use mpi_f08, only: MPI_Initialized
      integer(int32), intent(out) :: np, me
      integer(int32) :: ierr
      logical :: mpi_is_init

      ! Check if MPI is initialized first
      call MPI_Initialized(mpi_is_init, ierr)
      if (.not. mpi_is_init) then
         np = 1
         me = 0
         return
      end if

      call init_smp_comm()
      np = smp_comm%size()
      me = smp_comm%rank()
   end subroutine ddi_smp_nproc_impl

   !> Get offset into SMP array for this process
   !!
   !! @param handle SMP array handle
   !! @param addr_base Address of GAMESS's ADDR array
   !! @param offset Output: offset such that ADDR(offset) points to our data
   subroutine ddi_smp_offset_impl(handle, addr_base, offset)
      use iso_c_binding, only: c_loc, c_ptr
      integer(int32), intent(in) :: handle
      integer(int64), intent(in) :: addr_base
      integer(int64), intent(out) :: offset
      type(c_ptr) :: data_ptr
      integer(int64) :: data_addr

      call init_smp_comm()

      if (handle < 1 .or. handle > MAX_SMP_ARRAYS) then
         offset = 1
         return
      end if

      if (.not. smp_arrays(handle)%active) then
         offset = 1
         return
      end if

      if (.not. associated(smp_arrays(handle)%data)) then
         offset = 1
         return
      end if

      ! Get address of our allocated data
      data_ptr = c_loc(smp_arrays(handle)%data(1))
      data_addr = transfer(data_ptr, data_addr)

      ! Compute offset such that ADDR(offset) points to our data
      ! ADDR is CHARACTER*1, so offset is in bytes
      ! Fortran arrays are 1-indexed: ADDR(offset) is at address (addr_base + offset - 1)
      ! We want: addr_base + offset - 1 = data_addr
      ! So: offset = data_addr - addr_base + 1
      offset = data_addr - addr_base + 1

      print *, "DDI_SMP_OFFSET: handle=", handle, " addr_base=", addr_base, &
         " data_addr=", data_addr, " offset=", offset
   end subroutine ddi_smp_offset_impl

   !> Synchronize SMP processes
   subroutine ddi_smp_sync_impl()
      call init_smp_comm()
      call smp_comm%barrier()
   end subroutine ddi_smp_sync_impl

   !> Broadcast within SMP communicator
   !!
   !! @param buffer Data buffer
   !! @param n Number of elements
   !! @param root Root rank
   subroutine ddi_smp_bcast_impl(buffer, n, root)
      real(dp), intent(inout) :: buffer(*)
      integer(int32), intent(in) :: n, root
      integer(int32) :: ierr

      call init_smp_comm()
      call MPI_Bcast(buffer, n, MPI_DOUBLE_PRECISION, root, smp_comm%get(), ierr)
   end subroutine ddi_smp_bcast_impl

   !> Global sum within SMP communicator
   !!
   !! @param buffer Data buffer (in/out)
   !! @param n Number of elements
   subroutine ddi_smp_gsumf_impl(buffer, n)
      real(dp), intent(inout) :: buffer(*)
      integer(int32), intent(in) :: n
      real(dp), allocatable :: temp(:)

      call init_smp_comm()
      allocate (temp(n))
      temp = buffer(1:n)
      call allreduce(smp_comm, temp, n)
      buffer(1:n) = temp
      deallocate (temp)
   end subroutine ddi_smp_gsumf_impl

   ! ========================================================================
   ! Masters Operations
   ! ========================================================================

   !> Global sum among SMP masters (one per node)
   !!
   !! Only SMP masters (SMP_ME == 0) should call this.
   !!
   !! @param buffer Data buffer (in/out)
   !! @param n Number of elements
   subroutine ddi_masters_gsumf_impl(buffer, n)
      real(dp), intent(inout) :: buffer(*)
      integer(int32), intent(in) :: n
      real(dp), allocatable :: temp(:)
      integer(int32) :: smp_rank

      ! Initialize SMP comms if needed
      call init_smp_comm()

      smp_rank = smp_comm%rank()

      ! Only SMP masters participate
      if (smp_rank == 0 .and. smp_masters_comm%size() > 1) then
         allocate (temp(n))
         temp = buffer(1:n)
         call allreduce(smp_masters_comm, temp, n)
         buffer(1:n) = temp
         deallocate (temp)
      end if
   end subroutine ddi_masters_gsumf_impl

   !> Broadcast among SMP masters (one per node)
   !!
   !! Only SMP masters (SMP_ME == 0) should call this.
   !!
   !! @param buffer Data buffer
   !! @param n Number of elements
   !! @param root Root rank (within masters comm)
   subroutine ddi_masters_bcast_impl(buffer, n, root)
      real(dp), intent(inout) :: buffer(*)
      integer(int32), intent(in) :: n, root
      integer(int32) :: ierr, smp_rank

      ! Initialize SMP comms if needed
      call init_smp_comm()

      smp_rank = smp_comm%rank()

      ! Only SMP masters participate
      if (smp_rank == 0 .and. smp_masters_comm%size() > 1) then
         call MPI_Bcast(buffer, n, MPI_DOUBLE_PRECISION, root, smp_masters_comm%get(), ierr)
      end if
   end subroutine ddi_masters_bcast_impl

   ! ========================================================================
   ! Communicator-specific Get/Put
   ! ========================================================================

   !> Get data using a specific communicator
   !!
   !! @param handle Array handle
   !! @param ilo First row
   !! @param ihi Last row
   !! @param jlo First column
   !! @param jhi Last column
   !! @param buffer Output buffer
   !! @param comm_id Communicator ID
   subroutine ddi_get_comm_impl(handle, ilo, ihi, jlo, jhi, buffer, comm_id)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi, comm_id
      real(dp), intent(out) :: buffer(*)

      ! For now, just use the regular get - the comm_id affects scope
      ! which should already be set via ddi_scope
      call darray_get(handle, ilo, ihi, jlo, jhi, buffer)
   end subroutine ddi_get_comm_impl

   !> Put data using a specific communicator
   !!
   !! @param handle Array handle
   !! @param ilo First row
   !! @param ihi Last row
   !! @param jlo First column
   !! @param jhi Last column
   !! @param buffer Input buffer
   !! @param comm_id Communicator ID
   subroutine ddi_put_comm_impl(handle, ilo, ihi, jlo, jhi, buffer, comm_id)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi, comm_id
      real(dp), intent(in) :: buffer(*)

      ! For now, just use the regular put
      call darray_put(handle, ilo, ihi, jlo, jhi, buffer)
   end subroutine ddi_put_comm_impl

   !> Get data server ID (no-op for MPI-3 RMA - no data servers)
   !!
   !! @param dsid Output: data server ID (-1 = no data servers)
   subroutine ddi_get_dsid_impl(dsid)
      integer(int32), intent(out) :: dsid
      ! No data servers in MPI-3 RMA implementation
      dsid = -1
   end subroutine ddi_get_dsid_impl

   !> Scatter accumulate operation
   !!
   !! Accumulates scattered data into a distributed array.
   !! @param handle Array handle
   !! @param buff Data buffer
   !! @param nelem Number of elements
   !! @param rows Row indices (1-indexed)
   !! @param cols Column indices (1-indexed)
   subroutine ddi_scatter_acc_impl(handle, buff, nelem, rows, cols)
      integer(int32), intent(in) :: handle, nelem
      real(dp), intent(in) :: buff(*)
      integer(int32), intent(in) :: rows(*), cols(*)
      integer(int32) :: i
      real(dp) :: val(1)

      ! Scatter accumulate each element individually
      ! This is inefficient but correct - a production version would batch by owner
      do i = 1, nelem
         val(1) = buff(i)
         call darray_acc(handle, rows(i), rows(i), cols(i), cols(i), val)
      end do
   end subroutine ddi_scatter_acc_impl

   ! ========================================================================
   ! Process-based Dynamic Load Balancing
   ! ========================================================================

   !> Create a process-based DLB counter
   !!
   !! @param handle Output: DLB handle
   subroutine ddi_procdlb_create_impl(handle)
      integer(int32), intent(out) :: handle
      integer(int32) :: i

      handle = -1
      do i = 1, MAX_PROCDLB
         if (.not. procdlb_state(i)%active) then
            handle = i
            exit
         end if
      end do

      if (handle < 0) then
         error stop "ddi_procdlb_create: no free slots"
      end if

      procdlb_state(handle)%active = .true.
      procdlb_state(handle)%counter = 0
      procdlb_state(handle)%total_tasks = 0
      procdlb_state(handle)%last_iproc = -1
   end subroutine ddi_procdlb_create_impl

   !> Destroy a process-based DLB counter
   !!
   !! @param handle DLB handle
   subroutine ddi_procdlb_destroy_impl(handle)
      integer(int32), intent(in) :: handle

      if (handle < 1 .or. handle > MAX_PROCDLB) return
      procdlb_state(handle)%active = .false.
   end subroutine ddi_procdlb_destroy_impl

   !> Reset a process-based DLB counter
   !!
   !! @param handle DLB handle
   subroutine ddi_procdlb_reset_impl(handle)
      integer(int32), intent(in) :: handle

      if (handle < 1 .or. handle > MAX_PROCDLB) return
      if (.not. procdlb_state(handle)%active) return

      procdlb_state(handle)%counter = 0
      procdlb_state(handle)%last_iproc = -1
   end subroutine ddi_procdlb_reset_impl

   !> Get next task from process-based DLB
   !!
   !! Returns sequential task indices (0, 1, 2, ...) for processing columns
   !! owned by rank iproc. Counter resets when iproc changes.
   !!
   !! This is a simplified implementation without distributed work sharing.
   !! Each rank gets sequential indices within each owner's columns.
   !!
   !! @param handle DLB handle
   !! @param iproc Rank whose columns we're processing
   !! @param counter Output: next task index (0-indexed)
   subroutine ddi_procdlb_next_impl(handle, iproc, counter)
      integer(int32), intent(in) :: handle, iproc
      integer(int64), intent(out) :: counter

      if (handle < 1 .or. handle > MAX_PROCDLB) then
         counter = -1
         return
      end if
      if (.not. procdlb_state(handle)%active) then
         counter = -1
         return
      end if

      ! Reset counter when we switch to a different iproc (new rank's columns)
      if (iproc /= procdlb_state(handle)%last_iproc) then
         procdlb_state(handle)%counter = 0
         procdlb_state(handle)%last_iproc = iproc
      end if

      ! Return current counter value and increment
      counter = procdlb_state(handle)%counter
      procdlb_state(handle)%counter = procdlb_state(handle)%counter + 1_int64
   end subroutine ddi_procdlb_next_impl

   ! ========================================================================
   ! Miscellaneous
   ! ========================================================================

   !> Get node-level distribution info
   !!
   !! @param handle Array handle
   !! @param node Node index
   !! @param ilo Output: first row
   !! @param ihi Output: last row
   !! @param jlo Output: first column
   !! @param jhi Output: last column
   subroutine ddi_ndistrib_impl(handle, node, ilo, ihi, jlo, jhi)
      integer(int32), intent(in) :: handle, node
      integer(int32), intent(out) :: ilo, ihi, jlo, jhi

      ! For now, map node to rank and use regular distrib
      ! A proper implementation would aggregate across all ranks on a node
      call darray_distrib(handle, node, ilo, ihi, jlo, jhi)
   end subroutine ddi_ndistrib_impl

   !> Get DDI implementation level
   !!
   !! @param level Output: DDI level (0=none, 1=basic, 2=full)
   subroutine ddi_level_impl(level)
      integer(int32), intent(out) :: level
      ! Return 2 for full DDI support
      level = 2
   end subroutine ddi_level_impl

   ! ========================================================================
   ! No-op Stubs (rarely used internal/debug functions)
   ! ========================================================================

   !> End profiling section (no-op)
   subroutine ddi_pend_impl()
      ! No-op: profiling not implemented
   end subroutine ddi_pend_impl

   !> Begin profiling section (no-op)
   !!
   !! DDI_PBEG is a deprecated alias for DDI_Init in original DDI
   subroutine ddi_pbeg_impl(section)
      integer(int32), intent(in) :: section
      call ddi_init_impl()
   end subroutine ddi_pbeg_impl

   !> Test address validity (no-op, always returns valid)
   subroutine ddi_addr_test_impl(addr, valid)
      integer(int64), intent(in) :: addr
      integer(int32), intent(out) :: valid
      ! Always return valid
      valid = 1
   end subroutine ddi_addr_test_impl

   !> Generic global sum (no-op, use ddi_gsumf or ddi_gsumi instead)
   subroutine ddi_gsum_impl(tag, dtype, buffer, n)
      integer(int32), intent(in) :: tag, n
      character(len=1), intent(in) :: dtype
      real(dp), intent(inout) :: buffer(*)

      ! Route to appropriate typed version
      select case (dtype)
      case ('F', 'f')
         call ddi_gsumf_impl(tag, buffer, n)
      case ('I', 'i')
         ! For integers, would need separate handling
         ! Just do float version for now
         call ddi_gsumf_impl(tag, buffer, n)
      end select
   end subroutine ddi_gsum_impl

   !> Inner group DLB reset (no-op, use ddi_gdlbreset instead)
   subroutine ddi_gdlbreset_inner_impl()
      call ddi_gdlbreset_impl()
   end subroutine ddi_gdlbreset_inner_impl

   !> Inner group DLB next (no-op, use ddi_gdlbnext instead)
   subroutine ddi_gdlbnext_inner_impl(counter)
      integer(int64), intent(out) :: counter
      call ddi_gdlbnext_impl(counter)
   end subroutine ddi_gdlbnext_inner_impl

end module ddi_impl
