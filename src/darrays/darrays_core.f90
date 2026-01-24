!> Core distributed array operations
!!
!! Provides create, destroy, get, put, accumulate operations for
!! DDI-style distributed 2D arrays built on MPI-3 RMA.
!! Supports dp, sp, i32, and i64 data types.
!! Scope-aware: arrays are created in and use the working communicator.
module darrays_core
   use pic_types, only: int32, int64, sp, dp
   use pic_mpi_lib, only: comm_t, win_t, win_allocate, request_t, waitall
   use mpi_f08, only: MPI_ADDRESS_KIND
   use darrays_types, only: darray_t, DTYPE_DP, DTYPE_SP, DTYPE_I32, DTYPE_I64
   use darrays_distrib, only: calculate_distribution, get_owner, get_local_offset
   use groups, only: get_working_comm, groups_get_comm, groups_get_working_comm_id
   implicit none
   private

   public :: darrays_init, darrays_finalize
   public :: darray_create, darray_destroy, darray_distrib
   public :: darray_get, darray_put, darray_acc
   public :: darrays_get_comm

   integer(int32), parameter :: MAX_ARRAYS = 100
   type(darray_t), target, save :: registry(MAX_ARRAYS)
   logical, save :: initialized = .false.

   ! Generic interface for create - disambiguated by init_val type (required)
   interface darray_create
      module procedure darray_create_dp
      module procedure darray_create_sp
      module procedure darray_create_i32
      module procedure darray_create_i64
   end interface darray_create

   ! Generic interfaces for get/put/acc (disambiguated by buffer type)
   interface darray_get
      module procedure darray_get_dp
      module procedure darray_get_sp
      module procedure darray_get_i32
      module procedure darray_get_i64
   end interface darray_get

   interface darray_put
      module procedure darray_put_dp
      module procedure darray_put_sp
      module procedure darray_put_i32
      module procedure darray_put_i64
   end interface darray_put

   interface darray_acc
      module procedure darray_acc_dp
      module procedure darray_acc_sp
      module procedure darray_acc_i32
      module procedure darray_acc_i64
   end interface darray_acc

contains

   !> Initialize the distributed arrays module
   !!
   !! Note: The comm parameter is now ignored. Arrays use the working
   !! communicator from the groups module. Call groups_init() first.
   subroutine darrays_init(comm)
      type(comm_t), intent(in) :: comm
      integer(int32) :: i

      ! Note: comm parameter kept for API compatibility but not used
      ! Arrays now use working communicator from groups module

      do i = 1, MAX_ARRAYS
         registry(i)%active = .false.
         registry(i)%handle = -1
         registry(i)%comm_id = 0
      end do
      initialized = .true.
   end subroutine darrays_init

   !> Finalize the distributed arrays module
   subroutine darrays_finalize()
      integer(int32) :: i

      if (.not. initialized) return

      do i = 1, MAX_ARRAYS
         if (registry(i)%active) then
            call darray_destroy(i)
         end if
      end do

      initialized = .false.
   end subroutine darrays_finalize

   !> Get communicator for an array by its comm_id
   !!
   !! @param comm_id Communicator ID stored in array
   !! @return The communicator
   function darrays_get_comm(comm_id) result(comm)
      integer(int32), intent(in) :: comm_id
      type(comm_t) :: comm
      comm = groups_get_comm(comm_id)
   end function darrays_get_comm

   !> Find a free registry slot
   function find_free_slot() result(slot)
      integer(int32) :: slot
      integer(int32) :: i

      slot = -1
      do i = 1, MAX_ARRAYS
         if (.not. registry(i)%active) then
            slot = i
            exit
         end if
      end do

      if (slot < 0) then
         error stop "darrays: no free registry slots"
      end if
   end function find_free_slot

   ! ========================================================================
   ! Double precision (dp) operations
   ! ========================================================================

   subroutine darray_create_dp(nrows, ncols, handle, init_val)
      integer(int32), intent(in) :: nrows, ncols
      integer(int32), intent(out) :: handle
      real(dp), intent(in) :: init_val
      type(darray_t), pointer :: arr
      type(comm_t) :: work_comm
      integer(int32) :: slot, length

      slot = find_free_slot()
      arr => registry(slot)
      arr%handle = slot
      arr%dtype = DTYPE_DP
      arr%nrows = nrows
      arr%ncols = ncols
      arr%comm_id = groups_get_working_comm_id()
      work_comm = get_working_comm()

      call calculate_distribution(ncols, work_comm%size(), work_comm%rank(), &
                                  arr%my_first_col, arr%my_ncols)

      arr%local_size = int(nrows, int64)*int(arr%my_ncols, int64)
      length = int(arr%local_size, int32)
      if (length > 0) then
         call win_allocate(work_comm, length, arr%data_dp, arr%win)
      else
         call win_allocate(work_comm, 1_int32, arr%data_dp, arr%win)
      end if

      arr%data_dp = init_val

      arr%active = .true.
      handle = slot
      call work_comm%barrier()
   end subroutine darray_create_dp

   subroutine darray_get_dp(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      real(dp), intent(out) :: buffer(*)
      type(darray_t), pointer :: arr
      type(comm_t) :: arr_comm
      type(request_t), allocatable :: requests(:)
      integer(MPI_ADDRESS_KIND) :: disp
      integer(int32) :: j, col, owner, nrows_patch, buf_offset, req_count, nranks, ncols_req

      arr => registry(handle)
      arr_comm = groups_get_comm(arr%comm_id)
      nrows_patch = ihi - ilo + 1
      nranks = arr_comm%size()
      ncols_req = jhi - jlo + 1

      allocate (requests(ncols_req))
      call arr%win%lock_all()

      req_count = 0
      do j = jlo, jhi
         col = j - 1
         owner = get_owner(arr%ncols, nranks, col)
         disp = get_local_offset(arr%nrows, arr%ncols, nranks, ilo - 1, col)
         buf_offset = (j - jlo)*nrows_patch + 1
         req_count = req_count + 1
         call arr%win%rget_dp(owner, disp, nrows_patch, buffer(buf_offset), requests(req_count))
      end do

      call arr%win%flush_all()
      call waitall(requests(1:req_count))
      call arr%win%unlock_all()
      deallocate (requests)
   end subroutine darray_get_dp

   subroutine darray_put_dp(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      real(dp), intent(in) :: buffer(*)
      type(darray_t), pointer :: arr
      type(comm_t) :: arr_comm
      type(request_t), allocatable :: requests(:)
      integer(MPI_ADDRESS_KIND) :: disp
      integer(int32) :: j, col, owner, nrows_patch, buf_offset, req_count, nranks

      arr => registry(handle)
      arr_comm = groups_get_comm(arr%comm_id)
      nrows_patch = ihi - ilo + 1
      nranks = arr_comm%size()

      allocate (requests(jhi - jlo + 1))
      call arr%win%lock_all()

      req_count = 0
      do j = jlo, jhi
         col = j - 1
         owner = get_owner(arr%ncols, nranks, col)
         disp = get_local_offset(arr%nrows, arr%ncols, nranks, ilo - 1, col)
         buf_offset = (j - jlo)*nrows_patch + 1
         req_count = req_count + 1
         call arr%win%rput_dp(owner, disp, nrows_patch, buffer(buf_offset), requests(req_count))
      end do

      call arr%win%flush_all()
      call waitall(requests(1:req_count))
      call arr%win%unlock_all()
      deallocate (requests)
   end subroutine darray_put_dp

   subroutine darray_acc_dp(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      real(dp), intent(in) :: buffer(*)
      type(darray_t), pointer :: arr
      type(comm_t) :: arr_comm
      integer(MPI_ADDRESS_KIND) :: disp
      integer(int32) :: j, col, owner, nrows_patch, buf_offset, nranks

      arr => registry(handle)
      arr_comm = groups_get_comm(arr%comm_id)
      nrows_patch = ihi - ilo + 1
      nranks = arr_comm%size()

      do j = jlo, jhi
         col = j - 1
         owner = get_owner(arr%ncols, nranks, col)
         disp = get_local_offset(arr%nrows, arr%ncols, nranks, ilo - 1, col)
         buf_offset = (j - jlo)*nrows_patch + 1
         call arr%win%lock(owner)
         call arr%win%accumulate_dp(owner, disp, nrows_patch, buffer(buf_offset))
         call arr%win%unlock(owner)
      end do
   end subroutine darray_acc_dp

   ! ========================================================================
   ! Single precision (sp) operations
   ! ========================================================================

   subroutine darray_create_sp(nrows, ncols, handle, init_val)
      integer(int32), intent(in) :: nrows, ncols
      integer(int32), intent(out) :: handle
      real(sp), intent(in) :: init_val
      type(darray_t), pointer :: arr
      type(comm_t) :: work_comm
      integer(int32) :: slot, length

      slot = find_free_slot()
      arr => registry(slot)
      arr%handle = slot
      arr%dtype = DTYPE_SP
      arr%nrows = nrows
      arr%ncols = ncols
      arr%comm_id = groups_get_working_comm_id()
      work_comm = get_working_comm()

      call calculate_distribution(ncols, work_comm%size(), work_comm%rank(), &
                                  arr%my_first_col, arr%my_ncols)

      arr%local_size = int(nrows, int64)*int(arr%my_ncols, int64)
      length = int(arr%local_size, int32)
      if (length > 0) then
         call win_allocate(work_comm, length, arr%data_sp, arr%win)
      else
         call win_allocate(work_comm, 1_int32, arr%data_sp, arr%win)
      end if

      arr%data_sp = init_val

      arr%active = .true.
      handle = slot
      call work_comm%barrier()
   end subroutine darray_create_sp

   subroutine darray_get_sp(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      real(sp), intent(out) :: buffer(*)
      type(darray_t), pointer :: arr
      type(comm_t) :: arr_comm
      type(request_t), allocatable :: requests(:)
      integer(MPI_ADDRESS_KIND) :: disp
      integer(int32) :: j, col, owner, nrows_patch, buf_offset, req_count, nranks, ncols_req

      arr => registry(handle)
      arr_comm = groups_get_comm(arr%comm_id)
      nrows_patch = ihi - ilo + 1
      nranks = arr_comm%size()
      ncols_req = jhi - jlo + 1

      allocate (requests(ncols_req))
      call arr%win%lock_all()

      req_count = 0
      do j = jlo, jhi
         col = j - 1
         owner = get_owner(arr%ncols, nranks, col)
         disp = get_local_offset(arr%nrows, arr%ncols, nranks, ilo - 1, col)
         buf_offset = (j - jlo)*nrows_patch + 1
         req_count = req_count + 1
         call arr%win%rget_sp(owner, disp, nrows_patch, buffer(buf_offset), requests(req_count))
      end do

      call arr%win%flush_all()
      call waitall(requests(1:req_count))
      call arr%win%unlock_all()
      deallocate (requests)
   end subroutine darray_get_sp

   subroutine darray_put_sp(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      real(sp), intent(in) :: buffer(*)
      type(darray_t), pointer :: arr
      type(comm_t) :: arr_comm
      type(request_t), allocatable :: requests(:)
      integer(MPI_ADDRESS_KIND) :: disp
      integer(int32) :: j, col, owner, nrows_patch, buf_offset, req_count, nranks

      arr => registry(handle)
      arr_comm = groups_get_comm(arr%comm_id)
      nrows_patch = ihi - ilo + 1
      nranks = arr_comm%size()

      allocate (requests(jhi - jlo + 1))
      call arr%win%lock_all()

      req_count = 0
      do j = jlo, jhi
         col = j - 1
         owner = get_owner(arr%ncols, nranks, col)
         disp = get_local_offset(arr%nrows, arr%ncols, nranks, ilo - 1, col)
         buf_offset = (j - jlo)*nrows_patch + 1
         req_count = req_count + 1
         call arr%win%rput_sp(owner, disp, nrows_patch, buffer(buf_offset), requests(req_count))
      end do

      call arr%win%flush_all()
      call waitall(requests(1:req_count))
      call arr%win%unlock_all()
      deallocate (requests)
   end subroutine darray_put_sp

   subroutine darray_acc_sp(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      real(sp), intent(in) :: buffer(*)
      type(darray_t), pointer :: arr
      type(comm_t) :: arr_comm
      integer(MPI_ADDRESS_KIND) :: disp
      integer(int32) :: j, col, owner, nrows_patch, buf_offset, nranks

      arr => registry(handle)
      arr_comm = groups_get_comm(arr%comm_id)
      nrows_patch = ihi - ilo + 1
      nranks = arr_comm%size()

      do j = jlo, jhi
         col = j - 1
         owner = get_owner(arr%ncols, nranks, col)
         disp = get_local_offset(arr%nrows, arr%ncols, nranks, ilo - 1, col)
         buf_offset = (j - jlo)*nrows_patch + 1
         call arr%win%lock(owner)
         call arr%win%accumulate_sp(owner, disp, nrows_patch, buffer(buf_offset))
         call arr%win%unlock(owner)
      end do
   end subroutine darray_acc_sp

   ! ========================================================================
   ! Integer32 (i32) operations
   ! ========================================================================

   subroutine darray_create_i32(nrows, ncols, handle, init_val)
      integer(int32), intent(in) :: nrows, ncols
      integer(int32), intent(out) :: handle
      integer(int32), intent(in) :: init_val
      type(darray_t), pointer :: arr
      type(comm_t) :: work_comm
      integer(int32) :: slot, length

      slot = find_free_slot()
      arr => registry(slot)
      arr%handle = slot
      arr%dtype = DTYPE_I32
      arr%nrows = nrows
      arr%ncols = ncols
      arr%comm_id = groups_get_working_comm_id()
      work_comm = get_working_comm()

      call calculate_distribution(ncols, work_comm%size(), work_comm%rank(), &
                                  arr%my_first_col, arr%my_ncols)

      arr%local_size = int(nrows, int64)*int(arr%my_ncols, int64)
      length = int(arr%local_size, int32)
      if (length > 0) then
         call win_allocate(work_comm, length, arr%data_i32, arr%win)
      else
         call win_allocate(work_comm, 1_int32, arr%data_i32, arr%win)
      end if

      arr%data_i32 = init_val

      arr%active = .true.
      handle = slot
      call work_comm%barrier()
   end subroutine darray_create_i32

   subroutine darray_get_i32(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      integer(int32), intent(out) :: buffer(*)
      type(darray_t), pointer :: arr
      type(comm_t) :: arr_comm
      type(request_t), allocatable :: requests(:)
      integer(MPI_ADDRESS_KIND) :: disp
      integer(int32) :: j, col, owner, nrows_patch, buf_offset, req_count, nranks, ncols_req

      arr => registry(handle)
      arr_comm = groups_get_comm(arr%comm_id)
      nrows_patch = ihi - ilo + 1
      nranks = arr_comm%size()
      ncols_req = jhi - jlo + 1

      allocate (requests(ncols_req))
      call arr%win%lock_all()

      req_count = 0
      do j = jlo, jhi
         col = j - 1
         owner = get_owner(arr%ncols, nranks, col)
         disp = get_local_offset(arr%nrows, arr%ncols, nranks, ilo - 1, col)
         buf_offset = (j - jlo)*nrows_patch + 1
         req_count = req_count + 1
         call arr%win%rget_i32(owner, disp, nrows_patch, buffer(buf_offset), requests(req_count))
      end do

      call arr%win%flush_all()
      call waitall(requests(1:req_count))
      call arr%win%unlock_all()
      deallocate (requests)
   end subroutine darray_get_i32

   subroutine darray_put_i32(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      integer(int32), intent(in) :: buffer(*)
      type(darray_t), pointer :: arr
      type(comm_t) :: arr_comm
      type(request_t), allocatable :: requests(:)
      integer(MPI_ADDRESS_KIND) :: disp
      integer(int32) :: j, col, owner, nrows_patch, buf_offset, req_count, nranks

      arr => registry(handle)
      arr_comm = groups_get_comm(arr%comm_id)
      nrows_patch = ihi - ilo + 1
      nranks = arr_comm%size()

      allocate (requests(jhi - jlo + 1))
      call arr%win%lock_all()

      req_count = 0
      do j = jlo, jhi
         col = j - 1
         owner = get_owner(arr%ncols, nranks, col)
         disp = get_local_offset(arr%nrows, arr%ncols, nranks, ilo - 1, col)
         buf_offset = (j - jlo)*nrows_patch + 1
         req_count = req_count + 1
         call arr%win%rput_i32(owner, disp, nrows_patch, buffer(buf_offset), requests(req_count))
      end do

      call arr%win%flush_all()
      call waitall(requests(1:req_count))
      call arr%win%unlock_all()
      deallocate (requests)
   end subroutine darray_put_i32

   subroutine darray_acc_i32(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      integer(int32), intent(in) :: buffer(*)
      type(darray_t), pointer :: arr
      type(comm_t) :: arr_comm
      integer(MPI_ADDRESS_KIND) :: disp
      integer(int32) :: j, col, owner, nrows_patch, buf_offset, nranks

      arr => registry(handle)
      arr_comm = groups_get_comm(arr%comm_id)
      nrows_patch = ihi - ilo + 1
      nranks = arr_comm%size()

      do j = jlo, jhi
         col = j - 1
         owner = get_owner(arr%ncols, nranks, col)
         disp = get_local_offset(arr%nrows, arr%ncols, nranks, ilo - 1, col)
         buf_offset = (j - jlo)*nrows_patch + 1
         call arr%win%lock(owner)
         call arr%win%accumulate_i32(owner, disp, nrows_patch, buffer(buf_offset))
         call arr%win%unlock(owner)
      end do
   end subroutine darray_acc_i32

   ! ========================================================================
   ! Integer64 (i64) operations
   ! ========================================================================

   subroutine darray_create_i64(nrows, ncols, handle, init_val)
      integer(int32), intent(in) :: nrows, ncols
      integer(int32), intent(out) :: handle
      integer(int64), intent(in) :: init_val
      type(darray_t), pointer :: arr
      type(comm_t) :: work_comm
      integer(int32) :: slot, length

      slot = find_free_slot()
      arr => registry(slot)
      arr%handle = slot
      arr%dtype = DTYPE_I64
      arr%nrows = nrows
      arr%ncols = ncols
      arr%comm_id = groups_get_working_comm_id()
      work_comm = get_working_comm()

      call calculate_distribution(ncols, work_comm%size(), work_comm%rank(), &
                                  arr%my_first_col, arr%my_ncols)

      arr%local_size = int(nrows, int64)*int(arr%my_ncols, int64)
      length = int(arr%local_size, int32)
      if (length > 0) then
         call win_allocate(work_comm, length, arr%data_i64, arr%win)
      else
         call win_allocate(work_comm, 1_int32, arr%data_i64, arr%win)
      end if

      arr%data_i64 = init_val

      arr%active = .true.
      handle = slot
      call work_comm%barrier()
   end subroutine darray_create_i64

   subroutine darray_get_i64(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      integer(int64), intent(out) :: buffer(*)
      type(darray_t), pointer :: arr
      type(comm_t) :: arr_comm
      type(request_t), allocatable :: requests(:)
      integer(MPI_ADDRESS_KIND) :: disp
      integer(int32) :: j, col, owner, nrows_patch, buf_offset, req_count, nranks, ncols_req

      arr => registry(handle)
      arr_comm = groups_get_comm(arr%comm_id)
      nrows_patch = ihi - ilo + 1
      nranks = arr_comm%size()
      ncols_req = jhi - jlo + 1

      allocate (requests(ncols_req))
      call arr%win%lock_all()

      req_count = 0
      do j = jlo, jhi
         col = j - 1
         owner = get_owner(arr%ncols, nranks, col)
         disp = get_local_offset(arr%nrows, arr%ncols, nranks, ilo - 1, col)
         buf_offset = (j - jlo)*nrows_patch + 1
         req_count = req_count + 1
         call arr%win%rget_i64(owner, disp, nrows_patch, buffer(buf_offset), requests(req_count))
      end do

      call arr%win%flush_all()
      call waitall(requests(1:req_count))
      call arr%win%unlock_all()
      deallocate (requests)
   end subroutine darray_get_i64

   subroutine darray_put_i64(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      integer(int64), intent(in) :: buffer(*)
      type(darray_t), pointer :: arr
      type(comm_t) :: arr_comm
      type(request_t), allocatable :: requests(:)
      integer(MPI_ADDRESS_KIND) :: disp
      integer(int32) :: j, col, owner, nrows_patch, buf_offset, req_count, nranks

      arr => registry(handle)
      arr_comm = groups_get_comm(arr%comm_id)
      nrows_patch = ihi - ilo + 1
      nranks = arr_comm%size()

      allocate (requests(jhi - jlo + 1))
      call arr%win%lock_all()

      req_count = 0
      do j = jlo, jhi
         col = j - 1
         owner = get_owner(arr%ncols, nranks, col)
         disp = get_local_offset(arr%nrows, arr%ncols, nranks, ilo - 1, col)
         buf_offset = (j - jlo)*nrows_patch + 1
         req_count = req_count + 1
         call arr%win%rput_i64(owner, disp, nrows_patch, buffer(buf_offset), requests(req_count))
      end do

      call arr%win%flush_all()
      call waitall(requests(1:req_count))
      call arr%win%unlock_all()
      deallocate (requests)
   end subroutine darray_put_i64

   subroutine darray_acc_i64(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      integer(int64), intent(in) :: buffer(*)
      type(darray_t), pointer :: arr
      type(comm_t) :: arr_comm
      integer(MPI_ADDRESS_KIND) :: disp
      integer(int32) :: j, col, owner, nrows_patch, buf_offset, nranks

      arr => registry(handle)
      arr_comm = groups_get_comm(arr%comm_id)
      nrows_patch = ihi - ilo + 1
      nranks = arr_comm%size()

      do j = jlo, jhi
         col = j - 1
         owner = get_owner(arr%ncols, nranks, col)
         disp = get_local_offset(arr%nrows, arr%ncols, nranks, ilo - 1, col)
         buf_offset = (j - jlo)*nrows_patch + 1
         call arr%win%lock(owner)
         call arr%win%accumulate_i64(owner, disp, nrows_patch, buffer(buf_offset))
         call arr%win%unlock(owner)
      end do
   end subroutine darray_acc_i64

   ! ========================================================================
   ! Common operations (type-independent)
   ! ========================================================================

   !> Destroy a distributed array
   subroutine darray_destroy(handle)
      integer(int32), intent(in) :: handle
      type(darray_t), pointer :: arr
      type(comm_t) :: arr_comm

      if (handle < 1 .or. handle > MAX_ARRAYS) return
      arr => registry(handle)
      if (.not. arr%active) return

      arr_comm = groups_get_comm(arr%comm_id)
      call arr_comm%barrier()
      call arr%win%finalize()

      ! Nullify the appropriate data pointer based on dtype
      select case (arr%dtype)
      case (DTYPE_DP)
         nullify (arr%data_dp)
      case (DTYPE_SP)
         nullify (arr%data_sp)
      case (DTYPE_I32)
         nullify (arr%data_i32)
      case (DTYPE_I64)
         nullify (arr%data_i64)
      end select

      arr%active = .false.
      arr%handle = -1
      arr%dtype = 0
      arr%comm_id = 0
   end subroutine darray_destroy

   !> Query distribution for a rank
   subroutine darray_distrib(handle, rank, ilo, ihi, jlo, jhi)
      integer(int32), intent(in) :: handle, rank
      integer(int32), intent(out) :: ilo, ihi, jlo, jhi
      type(darray_t), pointer :: arr
      type(comm_t) :: arr_comm
      integer(int32) :: first_col, ncols_owned

      arr => registry(handle)
      arr_comm = groups_get_comm(arr%comm_id)

      call calculate_distribution(arr%ncols, arr_comm%size(), rank, &
                                  first_col, ncols_owned)

      ilo = 1
      ihi = arr%nrows
      jlo = first_col + 1
      jhi = first_col + ncols_owned
   end subroutine darray_distrib

end module darrays_core
