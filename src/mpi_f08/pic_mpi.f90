!> Modern MPI wrapper module using mpi_f08 interface
!!
!! This module provides a high-level object-oriented interface to MPI
!! using the modern mpi_f08 bindings. It wraps MPI communicators and
!! requests into derived types with type-bound procedures.
!!
module pic_mpi_f08
   use pic_types, only: int32, int64, sp, dp
   use mpi_f08, only: MPI_Comm, MPI_Status, MPI_Request, MPI_COMM_NULL, MPI_COMM_WORLD, &
                      MPI_COMM_TYPE_SHARED, MPI_INFO_NULL, MPI_UNDEFINED, &
                      MPI_INTEGER, MPI_INTEGER8, MPI_STATUS_IGNORE, MPI_REQUEST_NULL, MPI_Comm_rank, &
                      MPI_Comm_size, MPI_Comm_dup, MPI_Barrier, &
                      MPI_Comm_split_type, MPI_Comm_split, MPI_Send, MPI_Recv, &
                      MPI_Isend, MPI_Irecv, MPI_Wait, MPI_Waitall, MPI_Test, &
                      MPI_Probe, MPI_Get_count, MPI_Iprobe, MPI_Comm_free, &
                      MPI_Abort, MPI_Allgather, MPI_Get_processor_name, &
                      MPI_Bcast, MPI_Init, MPI_Init_thread, MPI_Query_thread, MPI_Finalize, &
                      MPI_THREAD_SINGLE, MPI_THREAD_FUNNELED, MPI_THREAD_SERIALIZED, MPI_THREAD_MULTIPLE, &
                      MPI_ANY_SOURCE, MPI_ANY_TAG, &
                      MPI_MAX_PROCESSOR_NAME, MPI_LOGICAL, &
                      operator(==), operator(/=), MPI_DOUBLE_PRECISION, MPI_REAL, &
                      MPI_Win, MPI_Op, MPI_WIN_NULL, MPI_Win_create, MPI_Win_create_dynamic, &
                      MPI_Win_free, MPI_Win_fence, MPI_Win_lock, MPI_Win_unlock, &
                      MPI_Win_lock_all, MPI_Win_unlock_all, MPI_Win_flush, MPI_Win_flush_all, &
                      MPI_Rget, MPI_Rput, MPI_Win_allocate, &
                      MPI_Get, MPI_Put, MPI_Accumulate, MPI_Fetch_and_op, &
                      MPI_Allreduce, MPI_ADDRESS_KIND, MPI_LOCK_SHARED, MPI_SUM, MPI_MIN, MPI_MAX, MPI_IN_PLACE
   implicit none
   private

   public :: comm_t, comm_world, comm_null
   public :: send, recv, isend, irecv
   public :: request_t, wait, waitall, test
   public :: iprobe, abort_comm, allgather, get_processor_name, bcast
   public :: pic_mpi_init, pic_mpi_finalize, pic_mpi_query_thread_level
   public :: win_t, win_create, win_create_dynamic, win_allocate
   public :: allreduce

   ! Export MPI types and constants needed by applications
   public :: MPI_Status, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_MAX_PROCESSOR_NAME, MPI_Request
   public :: MPI_THREAD_SINGLE, MPI_THREAD_FUNNELED, MPI_THREAD_SERIALIZED, MPI_THREAD_MULTIPLE
   public :: MPI_SUM, MPI_MIN, MPI_MAX

   type :: request_t
   !! Wraps MPI_Request to provide object-oriented interface for
   !! non-blocking communication operations (isend, irecv)
      private
      type(MPI_Request) :: m_request = MPI_REQUEST_NULL !! Internal MPI request handle
      logical :: is_valid = .false. !! Validity flag
   contains
      procedure :: is_null => request_is_null !! Check if request is null
      procedure :: get => request_get !! Get underlying MPI_Request
      procedure :: free => request_free !! Free the request
   end type request_t

   !> MPI-3 Window type for one-sided communication (RMA)
   !!
   !! Wraps MPI_Win to provide object-oriented interface for
   !! Remote Memory Access (RMA) operations needed for DDI
   type :: win_t
      private
      type(MPI_Win) :: m_win = MPI_WIN_NULL
      logical :: is_valid = .false.
   contains
      procedure :: is_null => win_is_null
      procedure :: get_handle => win_get_handle
      procedure :: fence => win_fence
      procedure :: lock => win_lock
      procedure :: unlock => win_unlock
      procedure :: lock_all => win_lock_all
      procedure :: unlock_all => win_unlock_all
      procedure :: flush => win_flush
      procedure :: flush_all => win_flush_all
      ! Double precision (dp) RMA operations
      procedure :: get_dp => win_get_dp
      procedure :: put_dp => win_put_dp
      procedure :: rget_dp => win_rget_dp
      procedure :: rput_dp => win_rput_dp
      procedure :: accumulate_dp => win_accumulate_dp
      ! Single precision (sp) RMA operations
      procedure :: get_sp => win_get_sp
      procedure :: put_sp => win_put_sp
      procedure :: rget_sp => win_rget_sp
      procedure :: rput_sp => win_rput_sp
      procedure :: accumulate_sp => win_accumulate_sp
      ! Integer32 RMA operations
      procedure :: get_i32 => win_get_i32
      procedure :: put_i32 => win_put_i32
      procedure :: rget_i32 => win_rget_i32
      procedure :: rput_i32 => win_rput_i32
      procedure :: accumulate_i32 => win_accumulate_i32
      ! Integer64 RMA operations
      procedure :: get_i64 => win_get_i64
      procedure :: put_i64 => win_put_i64
      procedure :: rget_i64 => win_rget_i64
      procedure :: rput_i64 => win_rput_i64
      procedure :: accumulate_i64 => win_accumulate_i64
      procedure :: fetch_and_add_i64 => win_fetch_and_add_i64
      procedure :: finalize => win_finalize
   end type win_t

   !> MPI communicator wrapper type
   !!
   type :: comm_t
   !! Provides object-oriented interface to MPI communicators with
   !! type-bound procedures for common operations. Automatically caches
   !! rank and size information for efficient access.
      private
      type(MPI_Comm) :: m_comm = MPI_COMM_NULL !! Internal MPI communicator
      integer(int32) :: m_rank = -1 !! Cached rank in this communicator
      integer(int32) :: m_size = -1 !! Cached size of this communicator
      logical :: is_valid = .false. !! Validity flag
   contains
      procedure :: rank => comm_rank !! Get rank in communicator
      procedure :: size => m_size_func !! Get size of communicator
      procedure :: leader => comm_leader !! Check if this rank is leader (rank 0)
      procedure :: is_null => comm_is_null !! Check if communicator is null
      procedure :: get => comm_get !! Get underlying MPI_Comm

      procedure :: barrier => comm_barrier !! Synchronization barrier

      procedure :: split => comm_split_shared !! Split into shared memory communicators
      procedure :: split_by => comm_split_by_color !! Split communicator by color
      procedure :: discard_leader => comm_discard_leader !! Create communicator without leader
      procedure :: discard_to => comm_discard_to !! Create communicator with first N ranks
      procedure :: duplicate => comm_duplicate !! Duplicate communicator

      procedure :: finalize => comm_finalize !! Free communicator resources
   end type comm_t

   interface comm_world
      module procedure create_world_comm
   end interface

   interface comm_null
      module procedure create_null_comm
   end interface

   interface send
      module procedure :: comm_send_integer
      module procedure :: comm_send_integer_array
      module procedure :: comm_send_integer64
      module procedure :: comm_send_integer64_array
      module procedure :: comm_send_real_dp
      module procedure :: comm_send_real_dp_array
      module procedure :: comm_send_real_dp_array_2d
      module procedure :: comm_send_logical
   end interface send

   interface recv
      module procedure :: comm_recv_integer
      module procedure :: comm_recv_integer_array
      module procedure :: comm_recv_integer64
      module procedure :: comm_recv_integer64_array
      module procedure :: comm_recv_real_dp
      module procedure :: comm_recv_real_dp_array
      module procedure :: comm_recv_real_dp_array_2d
      module procedure :: comm_recv_logical
   end interface recv

   interface iprobe
      module procedure :: comm_iprobe
   end interface iprobe

   interface allgather
      module procedure :: comm_allgather_integer
   end interface allgather

   interface bcast
      module procedure :: comm_bcast_integer
      module procedure :: comm_bcast_integer64
      module procedure :: comm_bcast_real_dp
      module procedure :: comm_bcast_real_dp_array
   end interface bcast

   interface isend
      module procedure :: comm_isend_integer
      module procedure :: comm_isend_integer_array
      module procedure :: comm_isend_integer64
      module procedure :: comm_isend_integer64_array
      module procedure :: comm_isend_real_dp
      module procedure :: comm_isend_real_dp_array
      module procedure :: comm_isend_real_dp_array_2d
      module procedure :: comm_isend_logical
   end interface isend

   interface irecv
      module procedure :: comm_irecv_integer
      module procedure :: comm_irecv_integer_array
      module procedure :: comm_irecv_integer64
      module procedure :: comm_irecv_integer64_array
      module procedure :: comm_irecv_real_dp
      module procedure :: comm_irecv_real_dp_array
      module procedure :: comm_irecv_real_dp_array_2d
      module procedure :: comm_irecv_logical
   end interface irecv

   interface wait
      module procedure :: request_wait
   end interface wait

   interface waitall
      module procedure :: request_waitall
   end interface waitall

   interface test
      module procedure :: request_test
   end interface test

   interface win_create
      module procedure create_win_dp_array
   end interface win_create

   interface win_create_dynamic
      module procedure create_win_dynamic
   end interface win_create_dynamic

   interface win_allocate
      module procedure create_win_allocate_dp_1d
      module procedure create_win_allocate_dp_2d
      module procedure create_win_allocate_sp_1d
      module procedure create_win_allocate_i32_1d
      module procedure create_win_allocate_i64_1d
   end interface win_allocate

   interface allreduce
      module procedure :: allreduce_dp
      module procedure :: allreduce_dp_array
      module procedure :: allreduce_i32
      module procedure :: allreduce_i32_array
      module procedure :: allreduce_dp_to
      module procedure :: allreduce_dp_array_to
   end interface allreduce

contains

   function create_comm_from_mpi(mpi_comm_in) result(comm)
   !! Internal helper function that wraps an MPI_Comm into a comm_t object
   !! and caches rank and size information
      type(MPI_Comm), intent(in) :: mpi_comm_in
      type(comm_t) :: comm
      integer(int32) :: ierr

      comm%m_comm = mpi_comm_in
      if (mpi_comm_in /= MPI_COMM_NULL) then
         call MPI_Comm_rank(comm%m_comm, comm%m_rank, ierr)
         call MPI_Comm_size(comm%m_comm, comm%m_size, ierr)
         comm%is_valid = .true.
      else
         comm%is_valid = .false.
      end if

   end function create_comm_from_mpi

   function create_world_comm() result(comm)
   !! Creates a new communicator that duplicates MPI_COMM_WORLD.
   !! This is the standard way to obtain a communicator for application use.
      type(comm_t) :: comm
      type(MPI_Comm) :: dup_comm
      integer(int32) :: ierr

      call MPI_Comm_dup(MPI_COMM_WORLD, dup_comm, ierr)
      comm = create_comm_from_mpi(dup_comm)

   end function create_world_comm

   function create_null_comm() result(comm)
   !! Creates an invalid/null communicator object that can be used
   !! for initialization or to represent absence of a communicator.
      type(comm_t) :: comm

      ! Explicitly initialize to null/invalid state
      comm%m_comm = MPI_COMM_NULL
      comm%m_rank = -1
      comm%m_size = -1
      comm%is_valid = .false.
   end function create_null_comm

   pure function comm_rank(this) result(rank)
   !!
   !! Returns the 0-indexed rank of the calling process
      class(comm_t), intent(in) :: this
      integer :: rank
      rank = this%m_rank
   end function comm_rank

   pure function m_size_func(this) result(size)
   !! Returns the number of processes in the communicator
      class(comm_t), intent(in) :: this
      integer :: size
      size = this%m_size
   end function m_size_func

   pure function comm_leader(this) result(is_leader)
   !! Returns true if the calling process has rank 0
      class(comm_t), intent(in) :: this
      logical :: is_leader
      is_leader = (this%m_rank == 0)
   end function comm_leader

   pure function comm_is_null(this) result(is_null)
      class(comm_t), intent(in) :: this
      logical :: is_null
      is_null = .not. this%is_valid
   end function comm_is_null

   function comm_get(this) result(mpi_comm_out)
      class(comm_t), intent(in) :: this
      type(MPI_Comm) :: mpi_comm_out

      if (.not. this%is_valid) then
         error stop "Cannot get MPI_Comm from null Comm"
      end if
      mpi_comm_out = this%m_comm
   end function comm_get

   subroutine comm_barrier(this)
   !! Blocks until all processes in the communicator have called barrier
      class(comm_t), intent(in) :: this
      integer(int32) :: ierr
      call MPI_Barrier(this%m_comm, ierr)
   end subroutine comm_barrier

   function comm_split_shared(this) result(new_comm)
   !! Creates a new communicator containing only processes that share
   !! memory with each other (typically processes on the same node)
      class(comm_t), intent(in) :: this
      type(comm_t) :: new_comm
      type(MPI_Comm) :: mpi_comm_new
      integer(int32) :: ierr

      if (.not. this%is_valid) then
         new_comm = comm_null()
         return
      end if

      call MPI_Comm_split_type(this%get(), MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, mpi_comm_new, ierr)
      new_comm = create_comm_from_mpi(mpi_comm_new)
   end function comm_split_shared

   function comm_split_by_color(this, color) result(new_comm)
   !! Partitions the communicator into disjoint subgroups based on color.
   !! Processes with the same color end up in the same new communicator.
      class(comm_t), intent(in) :: this
      integer, intent(in) :: color
      type(comm_t) :: new_comm
      type(MPI_Comm) :: mpi_comm_new
      integer(int32) :: ierr

      if (.not. this%is_valid) then
         new_comm = comm_null()
         return
      end if

      call MPI_Comm_split(this%m_comm, color, 0, mpi_comm_new, ierr)
      new_comm = create_comm_from_mpi(mpi_comm_new)
   end function comm_split_by_color

   function comm_discard_leader(this) result(new_comm)
      class(comm_t), intent(in) :: this
      type(comm_t) :: new_comm
      integer :: color

      if (.not. this%is_valid) then
         new_comm = comm_null()
         return
      end if

      if (this%rank() == 0) then
         color = MPI_UNDEFINED
      else
         color = 0
      end if
      new_comm = this%split_by(color)
   end function comm_discard_leader

   function comm_discard_to(this, num_ranks) result(new_comm)
      class(comm_t), intent(in) :: this
      integer, intent(in) :: num_ranks
      type(comm_t) :: new_comm
      integer :: color

      if (.not. this%is_valid) then
         new_comm = comm_null()
         return
      end if

      if (this%rank() < num_ranks) then
         color = 0
      else
         color = MPI_UNDEFINED
      end if
      new_comm = this%split_by(color)
   end function comm_discard_to

   function comm_duplicate(this) result(new_comm)
      class(comm_t), intent(in) :: this
      type(comm_t) :: new_comm
      type(MPI_Comm) :: mpi_comm_new
      integer(int32) :: ierr

      if (.not. this%is_valid) then
         new_comm = comm_null()
         return
      end if

      call MPI_Comm_dup(this%m_comm, mpi_comm_new, ierr)
      new_comm = create_comm_from_mpi(mpi_comm_new)
   end function comm_duplicate

   subroutine comm_send_integer(comm, data, dest, tag)
   !! Blocking send of an integer to specified destination
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: data
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr

      call MPI_Send(data, 1, MPI_INTEGER, dest, tag, comm%m_comm, ierr)
   end subroutine comm_send_integer

   subroutine comm_send_integer_array(comm, data, dest, tag)
   !! Blocking send of an integer array to specified destination
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: data(:)
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr

      call MPI_Send(data, size(data), MPI_INTEGER, dest, tag, comm%m_comm, ierr)
   end subroutine comm_send_integer_array

   subroutine comm_send_integer64(comm, data, dest, tag)
   !! Blocking send of an integer64 to specified destination
      type(comm_t), intent(in) :: comm
      integer(int64), intent(in) :: data
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr

      call MPI_Send(data, 1, MPI_INTEGER8, dest, tag, comm%m_comm, ierr)
   end subroutine comm_send_integer64

   subroutine comm_send_integer64_array(comm, data, dest, tag)
   !! Blocking send of an integer64 array to specified destination
      type(comm_t), intent(in) :: comm
      integer(int64), intent(in) :: data(:)
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr

      call MPI_Send(data, size(data), MPI_INTEGER8, dest, tag, comm%m_comm, ierr)
   end subroutine comm_send_integer64_array

   subroutine comm_send_real_dp(comm, data, dest, tag)
      !! Blocking send of a single double precision real to specified destination
      type(comm_t), intent(in) :: comm
      real(dp), intent(in) :: data
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr

      call MPI_Send(data, 1, MPI_DOUBLE_PRECISION, dest, tag, comm%m_comm, ierr)
   end subroutine comm_send_real_dp

   subroutine comm_send_real_dp_array(comm, data, dest, tag)
      !! Blocking send of a double precision real array to specified destination
      type(comm_t), intent(in) :: comm
      real(dp), intent(in) :: data(:)
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr

      call MPI_Send(data, size(data), MPI_DOUBLE_PRECISION, dest, tag, comm%m_comm, ierr)
   end subroutine comm_send_real_dp_array

   subroutine comm_send_real_dp_array_2d(comm, data, dest, tag)
   !! Blocking send of a 2D double precision real array to specified destination
      type(comm_t), intent(in) :: comm
      real(dp), intent(in) :: data(:, :)
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr, dim1, dim2

      ! Send dimensions first
      dim1 = size(data, 1)
      dim2 = size(data, 2)
      call MPI_Send(dim1, 1, MPI_INTEGER, dest, tag, comm%m_comm, ierr)
      call MPI_Send(dim2, 1, MPI_INTEGER, dest, tag, comm%m_comm, ierr)

      ! Send data
      call MPI_Send(data, size(data), MPI_DOUBLE_PRECISION, dest, tag, comm%m_comm, ierr)
   end subroutine comm_send_real_dp_array_2d

   subroutine comm_send_logical(comm, data, dest, tag)
   !! Blocking send of a logical value to specified destination
      type(comm_t), intent(in) :: comm
      logical, intent(in) :: data
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr

      call MPI_Send(data, 1, MPI_LOGICAL, dest, tag, comm%m_comm, ierr)
   end subroutine comm_send_logical

   subroutine comm_recv_integer(comm, data, source, tag, status)
   !! Blocking receive of an integer from specified source.
   !! Use MPI_ANY_SOURCE or MPI_ANY_TAG for wildcards.
      type(comm_t), intent(in) :: comm
      integer(int32), intent(out) :: data
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr
      type(MPI_Status), intent(out), optional :: status

      if (present(status)) then
         call MPI_Recv(data, 1, MPI_INTEGER, source, tag, comm%m_comm, status, ierr)
      else
         call MPI_Recv(data, 1, MPI_INTEGER, source, tag, comm%m_comm, MPI_STATUS_IGNORE, ierr)
      end if
   end subroutine comm_recv_integer

   subroutine comm_recv_integer_array(comm, data, source, tag, status)
   !! Blocking receive of an integer array from specified source.
      type(comm_t), intent(in) :: comm
      integer(int32), allocatable, intent(out) :: data(:)
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(MPI_Status), intent(out) :: status
      integer(int32) :: count
      integer(int32) :: ierr

      ! First probe to get message size
      call MPI_Probe(source, tag, comm%m_comm, status, ierr)
      call MPI_Get_count(status, MPI_INTEGER, count, ierr)

      ! Allocate and receive
      allocate (data(count))
      call MPI_Recv(data, count, MPI_INTEGER, source, tag, comm%m_comm, MPI_STATUS_IGNORE, ierr)
   end subroutine comm_recv_integer_array

   subroutine comm_recv_integer64(comm, data, source, tag, status)
   !! Blocking receive of an integer64 from specified source.
   !! Use MPI_ANY_SOURCE or MPI_ANY_TAG for wildcards.
      type(comm_t), intent(in) :: comm
      integer(int64), intent(out) :: data
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr
      type(MPI_Status), intent(out), optional :: status

      if (present(status)) then
         call MPI_Recv(data, 1, MPI_INTEGER8, source, tag, comm%m_comm, status, ierr)
      else
         call MPI_Recv(data, 1, MPI_INTEGER8, source, tag, comm%m_comm, MPI_STATUS_IGNORE, ierr)
      end if
   end subroutine comm_recv_integer64

   subroutine comm_recv_integer64_array(comm, data, source, tag, status)
   !! Blocking receive of an integer64 array from specified source.
   !! Array is automatically allocated to the correct size.
      type(comm_t), intent(in) :: comm
      integer(int64), allocatable, intent(out) :: data(:)
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(MPI_Status), intent(out) :: status
      integer(int32) :: count
      integer(int32) :: ierr

      ! First probe to get message size
      call MPI_Probe(source, tag, comm%m_comm, status, ierr)
      call MPI_Get_count(status, MPI_INTEGER8, count, ierr)

      ! Allocate and receive
      allocate (data(count))
      call MPI_Recv(data, count, MPI_INTEGER8, source, tag, comm%m_comm, MPI_STATUS_IGNORE, ierr)
   end subroutine comm_recv_integer64_array

   subroutine comm_recv_real_dp(comm, data, source, tag, status)
      !! Blocking receive of a single double precision real from specified source.
      type(comm_t), intent(in) :: comm
      real(dp), intent(out) :: data
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr
      type(MPI_Status), intent(out), optional :: status

      if (present(status)) then
         call MPI_Recv(data, 1, MPI_DOUBLE_PRECISION, source, tag, comm%m_comm, status, ierr)
      else
         call MPI_Recv(data, 1, MPI_DOUBLE_PRECISION, source, tag, comm%m_comm, MPI_STATUS_IGNORE, ierr)
      end if
   end subroutine comm_recv_real_dp

   subroutine comm_recv_real_dp_array(comm, data, source, tag, status)
      !! Blocking receive of a double precision real array from specified source.
      type(comm_t), intent(in) :: comm
      real(dp), allocatable, intent(out) :: data(:)
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(MPI_Status) :: status
      integer(int32) :: count
      integer(int32) :: ierr

      ! First probe to get message size
      call MPI_Probe(source, tag, comm%m_comm, status, ierr)
      call MPI_Get_count(status, MPI_DOUBLE_PRECISION, count, ierr)

      ! Allocate and receive
      allocate (data(count))
      call MPI_Recv(data, count, MPI_DOUBLE_PRECISION, source, tag, comm%m_comm, MPI_STATUS_IGNORE, ierr)
   end subroutine comm_recv_real_dp_array

   subroutine comm_recv_real_dp_array_2d(comm, data, source, tag, status)
   !! Blocking receive of a 2D allocatable double precision real array
      type(comm_t), intent(in) :: comm
      real(dp), intent(inout), allocatable :: data(:, :)
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(MPI_Status), intent(out) :: status
      integer(int32) :: ierr, count, dim1, dim2

      ! Receive dimensions first
      call MPI_Recv(dim1, 1, MPI_INTEGER, source, tag, comm%m_comm, MPI_STATUS_IGNORE, ierr)
      call MPI_Recv(dim2, 1, MPI_INTEGER, source, tag, comm%m_comm, MPI_STATUS_IGNORE, ierr)

      ! Allocate array with received dimensions
      if (.not. allocated(data)) then
         allocate (data(dim1, dim2))
      end if

      ! Receive data
      count = dim1*dim2
      call MPI_Recv(data, count, MPI_DOUBLE_PRECISION, source, tag, comm%m_comm, status, ierr)
   end subroutine comm_recv_real_dp_array_2d

   subroutine comm_recv_logical(comm, data, source, tag, status)
   !! Blocking receive of a logical value from specified source
      type(comm_t), intent(in) :: comm
      logical, intent(out) :: data
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(MPI_Status), intent(out), optional :: status
      integer(int32) :: ierr

      if (present(status)) then
         call MPI_Recv(data, 1, MPI_LOGICAL, source, tag, comm%m_comm, status, ierr)
      else
         call MPI_Recv(data, 1, MPI_LOGICAL, source, tag, comm%m_comm, MPI_STATUS_IGNORE, ierr)
      end if
   end subroutine comm_recv_logical

   subroutine comm_iprobe(comm, source, tag, message_pending, status)
      !! Non-blocking probe for incoming messages
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      logical, intent(out) :: message_pending
      type(MPI_Status), intent(out) :: status
      integer(int32) :: ierr

      call MPI_Iprobe(source, tag, comm%m_comm, message_pending, status, ierr)
   end subroutine comm_iprobe

   subroutine comm_finalize(this)
      !! Frees the MPI communicator resources
      class(comm_t), intent(inout) :: this
      integer(int32) :: ierr

      if (this%is_valid .and. this%m_comm /= MPI_COMM_NULL) then
         call MPI_Comm_free(this%m_comm, ierr)
         this%is_valid = .false.
         this%m_comm = MPI_COMM_NULL
      end if
   end subroutine comm_finalize

   subroutine abort_comm(comm, errorcode)
      !! Aborts all processes in the communicator with the given error code
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: errorcode
      integer(int32) :: ierr

      call MPI_Abort(comm%m_comm, errorcode, ierr)
   end subroutine abort_comm

   subroutine comm_allgather_integer(comm, sendbuf, recvbuf)
      !! Gathers integer values from all processes in the communicator
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: sendbuf
      integer(int32), intent(out) :: recvbuf(:)
      integer(int32) :: ierr

      call MPI_Allgather(sendbuf, 1, MPI_INTEGER, recvbuf, 1, MPI_INTEGER, comm%m_comm, ierr)
   end subroutine comm_allgather_integer

   subroutine comm_bcast_integer(comm, buffer, count, root)
      !! Broadcasts integer data from root process to all processes in communicator
      type(comm_t), intent(in) :: comm
      integer(int32), intent(inout) :: buffer
      integer(int32), intent(in) :: count
      integer(int32), intent(in) :: root
      integer(int32) :: ierr

      call MPI_Bcast(buffer, count, MPI_INTEGER, root, comm%m_comm, ierr)
   end subroutine comm_bcast_integer

   subroutine comm_bcast_integer64(comm, buffer, count, root)
      !! Broadcasts integer64 data from root process to all processes in communicator
      type(comm_t), intent(in) :: comm
      integer(int64), intent(inout) :: buffer
      integer(int32), intent(in) :: count
      integer(int32), intent(in) :: root
      integer(int32) :: ierr

      call MPI_Bcast(buffer, count, MPI_INTEGER8, root, comm%m_comm, ierr)
   end subroutine comm_bcast_integer64

   subroutine comm_bcast_real_dp(comm, buffer, count, root)
      !! Broadcasts double precision data from root process to all processes in communicator
      type(comm_t), intent(in) :: comm
      real(dp), intent(inout) :: buffer
      integer(int32), intent(in) :: count
      integer(int32), intent(in) :: root
      integer(int32) :: ierr

      call MPI_Bcast(buffer, count, MPI_DOUBLE_PRECISION, root, comm%m_comm, ierr)
   end subroutine comm_bcast_real_dp

   subroutine comm_bcast_real_dp_array(comm, buffer, count, root)
      !! Broadcasts double precision array from root process to all processes in communicator
      type(comm_t), intent(in) :: comm
      real(dp), intent(inout) :: buffer(:)
      integer(int32), intent(in) :: count
      integer(int32), intent(in) :: root
      integer(int32) :: ierr

      call MPI_Bcast(buffer, count, MPI_DOUBLE_PRECISION, root, comm%m_comm, ierr)
   end subroutine comm_bcast_real_dp_array

   subroutine get_processor_name(name, namelen)
      !! Retrieves the name of the processor
      character(len=*), intent(inout) :: name
      integer(int32), intent(out) :: namelen
      integer(int32) :: ierr

      call MPI_Get_processor_name(name, namelen, ierr)
   end subroutine get_processor_name

   subroutine pic_mpi_init(requested_thread_level, provided_thread_level)
      !! Initialize MPI environment with optional threading support
      !!
      !! If no thread level is requested, uses MPI_THREAD_FUNNELED by default
      !! to allow OpenMP threading in compute_mbe_energy and similar functions.
      !!
      !! Thread levels:
      !!   MPI_THREAD_SINGLE: No threading support
      !!   MPI_THREAD_FUNNELED: Only main thread makes MPI calls (good for OpenMP)
      !!   MPI_THREAD_SERIALIZED: Multiple threads can make MPI calls, but not simultaneously
      !!   MPI_THREAD_MULTIPLE: Full thread safety
      integer(int32), intent(in), optional :: requested_thread_level
      integer(int32), intent(out), optional :: provided_thread_level
      integer(int32) :: ierr, requested, provided

      ! Default to FUNNELED for OpenMP compatibility
      if (present(requested_thread_level)) then
         requested = requested_thread_level
      else
         requested = MPI_THREAD_FUNNELED
      end if

      call MPI_Init_thread(requested, provided, ierr)

      ! Return the provided level if requested
      if (present(provided_thread_level)) then
         provided_thread_level = provided
      end if

      ! Warn if we didn't get what we asked for
      if (provided < requested .and. requested /= MPI_THREAD_SINGLE) then
         if (requested == MPI_THREAD_FUNNELED) then
            write (*, '(a)') "Warning: MPI_THREAD_FUNNELED requested but not provided."
            write (*, '(a)') "OpenMP threading in compute functions may not work correctly."
         else if (requested == MPI_THREAD_SERIALIZED) then
            write (*, '(a)') "Warning: MPI_THREAD_SERIALIZED requested but not provided."
         else if (requested == MPI_THREAD_MULTIPLE) then
            write (*, '(a)') "Warning: MPI_THREAD_MULTIPLE requested but not provided."
         end if
      end if
   end subroutine pic_mpi_init

   function pic_mpi_query_thread_level() result(thread_level)
      !! Query the current MPI thread support level
      use mpi_f08, only: MPI_Query_thread
      integer(int32) :: thread_level
      integer(int32) :: ierr

      call MPI_Query_thread(thread_level, ierr)
   end function pic_mpi_query_thread_level

   subroutine pic_mpi_finalize(ierr)
      !! Finalize MPI environment
      integer(int32), optional, intent(out) :: ierr
      integer(int32) :: ierr_local
      call MPI_Finalize(ierr_local)
      if (present(ierr)) ierr = ierr_local
   end subroutine pic_mpi_finalize

   ! ========================================================================
   ! Request type methods
   ! ========================================================================

   pure function request_is_null(this) result(is_null)
      !! Checks if the request is null/invalid
      class(request_t), intent(in) :: this
      logical :: is_null
      is_null = .not. this%is_valid
   end function request_is_null

   function request_get(this) result(mpi_request_out)
      !! Retrieves the underlying MPI_Request handle
      class(request_t), intent(in) :: this
      type(MPI_Request) :: mpi_request_out

      if (.not. this%is_valid) then
         error stop "Cannot get MPI_Request from null request"
      end if
      mpi_request_out = this%m_request
   end function request_get

   subroutine request_free(this)
      !! Frees the MPI request resources
      class(request_t), intent(inout) :: this
      this%m_request = MPI_REQUEST_NULL
      this%is_valid = .false.
   end subroutine request_free

   ! ========================================================================
   ! Non-blocking send operations
   ! ========================================================================

   subroutine comm_isend_integer(comm, data, dest, tag, request)
   !! Initiates a non-blocking send operation. The request must be
   !! waited on using wait() or test() before the buffer can be reused.
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: data
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      call MPI_Isend(data, 1, MPI_INTEGER, dest, tag, comm%m_comm, request%m_request, ierr)
      request%is_valid = .true.
   end subroutine comm_isend_integer

   subroutine comm_isend_integer_array(comm, data, dest, tag, request)
   !! Initiates a non-blocking send operation. The request must be
   !! waited on using wait() or test() before the buffer can be reused.
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: data(:)
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      call MPI_Isend(data, size(data), MPI_INTEGER, dest, tag, comm%m_comm, request%m_request, ierr)
      request%is_valid = .true.
   end subroutine comm_isend_integer_array

   subroutine comm_isend_integer64(comm, data, dest, tag, request)
   !! Initiates a non-blocking send operation. The request must be
   !! waited on using wait() or test() before the buffer can be reused.
      type(comm_t), intent(in) :: comm
      integer(int64), intent(in) :: data
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      call MPI_Isend(data, 1, MPI_INTEGER8, dest, tag, comm%m_comm, request%m_request, ierr)
      request%is_valid = .true.
   end subroutine comm_isend_integer64

   subroutine comm_isend_integer64_array(comm, data, dest, tag, request)
   !! Initiates a non-blocking send operation. The request must be
   !! waited on using wait() or test() before the buffer can be reused.
      type(comm_t), intent(in) :: comm
      integer(int64), intent(in) :: data(:)
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      call MPI_Isend(data, size(data), MPI_INTEGER8, dest, tag, comm%m_comm, request%m_request, ierr)
      request%is_valid = .true.
   end subroutine comm_isend_integer64_array

   subroutine comm_isend_real_dp(comm, data, dest, tag, request)
   !! Initiates a non-blocking send operation. The request must be
   !! waited on using wait() or test() before the buffer can be reused.
      type(comm_t), intent(in) :: comm
      real(dp), intent(in) :: data
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      call MPI_Isend(data, 1, MPI_DOUBLE_PRECISION, dest, tag, comm%m_comm, request%m_request, ierr)
      request%is_valid = .true.
   end subroutine comm_isend_real_dp

   subroutine comm_isend_real_dp_array(comm, data, dest, tag, request)
      !! Initiates a non-blocking send operation. The request must be
      !! waited on using wait() or test() before the buffer can be reused.
      type(comm_t), intent(in) :: comm
      real(dp), intent(in) :: data(:)
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      call MPI_Isend(data, size(data), MPI_DOUBLE_PRECISION, dest, tag, comm%m_comm, request%m_request, ierr)
      request%is_valid = .true.
   end subroutine comm_isend_real_dp_array

   subroutine comm_isend_real_dp_array_2d(comm, data, dest, tag, request)
   !! Non-blocking send of a 2D double precision real array
      type(comm_t), intent(in) :: comm
      real(dp), intent(in) :: data(:, :)
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr, dim1, dim2

      ! Send dimensions first (blocking - simple approach)
      dim1 = size(data, 1)
      dim2 = size(data, 2)
      call MPI_Send(dim1, 1, MPI_INTEGER, dest, tag, comm%m_comm, ierr)
      call MPI_Send(dim2, 1, MPI_INTEGER, dest, tag, comm%m_comm, ierr)

      ! Send data (non-blocking)
      call MPI_Isend(data, size(data), MPI_DOUBLE_PRECISION, dest, tag, comm%m_comm, request%m_request, ierr)
      request%is_valid = .true.
   end subroutine comm_isend_real_dp_array_2d

   subroutine comm_isend_logical(comm, data, dest, tag, request)
   !! Non-blocking send of a logical value
      type(comm_t), intent(in) :: comm
      logical, intent(in) :: data
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      call MPI_Isend(data, 1, MPI_LOGICAL, dest, tag, comm%m_comm, request%m_request, ierr)
      request%is_valid = .true.
   end subroutine comm_isend_logical

   ! ========================================================================
   ! Non-blocking receive operations
   ! ========================================================================

   subroutine comm_irecv_integer(comm, data, source, tag, request)
      !! Initiates a non-blocking receive operation. The request must be
      !! waited on using wait() or test() before the buffer can be used.
      type(comm_t), intent(in) :: comm
      integer(int32), intent(out) :: data
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      call MPI_Irecv(data, 1, MPI_INTEGER, source, tag, comm%m_comm, request%m_request, ierr)
      request%is_valid = .true.
   end subroutine comm_irecv_integer

   subroutine comm_irecv_integer_array(comm, data, source, tag, request)
      !! Initiates a non-blocking receive operation. The request must be
      !! waited on using wait() or test() before the buffer can be used.
      type(comm_t), intent(in) :: comm
      integer(int32), intent(out) :: data(:)
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      call MPI_Irecv(data, size(data), MPI_INTEGER, source, tag, comm%m_comm, request%m_request, ierr)
      request%is_valid = .true.
   end subroutine comm_irecv_integer_array

   subroutine comm_irecv_integer64(comm, data, source, tag, request)
   !! Initiates a non-blocking receive operation. The request must be
   !! waited on using wait() or test() before the buffer can be used.
      type(comm_t), intent(in) :: comm
      integer(int64), intent(out) :: data
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      call MPI_Irecv(data, 1, MPI_INTEGER8, source, tag, comm%m_comm, request%m_request, ierr)
      request%is_valid = .true.
   end subroutine comm_irecv_integer64

   subroutine comm_irecv_integer64_array(comm, data, source, tag, request)
   !! Initiates a non-blocking receive operation. The request must be
   !! waited on using wait() or test() before the buffer can be used.
      type(comm_t), intent(in) :: comm
      integer(int64), intent(out) :: data(:)
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      call MPI_Irecv(data, size(data), MPI_INTEGER8, source, tag, comm%m_comm, request%m_request, ierr)
      request%is_valid = .true.
   end subroutine comm_irecv_integer64_array

   subroutine comm_irecv_real_dp(comm, data, source, tag, request)
      type(comm_t), intent(in) :: comm
      real(dp), intent(out) :: data
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      call MPI_Irecv(data, 1, MPI_DOUBLE_PRECISION, source, tag, comm%m_comm, request%m_request, ierr)
      request%is_valid = .true.
   end subroutine comm_irecv_real_dp

   subroutine comm_irecv_real_dp_array(comm, data, source, tag, request)
      type(comm_t), intent(in) :: comm
      real(dp), intent(out) :: data(:)
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      call MPI_Irecv(data, size(data), MPI_DOUBLE_PRECISION, source, tag, comm%m_comm, request%m_request, ierr)
      request%is_valid = .true.
   end subroutine comm_irecv_real_dp_array

   subroutine comm_irecv_real_dp_array_2d(comm, data, source, tag, request)
   !! Non-blocking receive of a 2D allocatable double precision real array
      type(comm_t), intent(in) :: comm
      real(dp), intent(inout), allocatable :: data(:, :)
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr, dim1, dim2

      ! Receive dimensions first (blocking - needed to allocate)
      call MPI_Recv(dim1, 1, MPI_INTEGER, source, tag, comm%m_comm, MPI_STATUS_IGNORE, ierr)
      call MPI_Recv(dim2, 1, MPI_INTEGER, source, tag, comm%m_comm, MPI_STATUS_IGNORE, ierr)

      ! Allocate array with received dimensions
      if (.not. allocated(data)) then
         allocate (data(dim1, dim2))
      end if

      ! Receive data (non-blocking)
      call MPI_Irecv(data, dim1*dim2, MPI_DOUBLE_PRECISION, source, tag, comm%m_comm, request%m_request, ierr)
      request%is_valid = .true.
   end subroutine comm_irecv_real_dp_array_2d

   subroutine comm_irecv_logical(comm, data, source, tag, request)
   !! Non-blocking receive of a logical value
      type(comm_t), intent(in) :: comm
      logical, intent(out) :: data
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      call MPI_Irecv(data, 1, MPI_LOGICAL, source, tag, comm%m_comm, request%m_request, ierr)
      request%is_valid = .true.
   end subroutine comm_irecv_logical

   ! ========================================================================
   ! Request completion operations
   ! ========================================================================

   subroutine request_wait(request, status)
   !! Blocks until the operation associated with the request completes.
   !! The request is freed after completion.
      type(request_t), intent(inout) :: request
      type(MPI_Status), intent(out), optional :: status
      integer(int32) :: ierr

      if (.not. request%is_valid) then
         error stop "Cannot wait on null request"
      end if

      if (present(status)) then
         call MPI_Wait(request%m_request, status, ierr)
      else
         call MPI_Wait(request%m_request, MPI_STATUS_IGNORE, ierr)
      end if

      call request%free()
   end subroutine request_wait

   subroutine request_waitall(requests, statuses)
   !! Blocks until all operations in the request array complete.
   !! All requests are freed after completion.
      type(request_t), intent(inout) :: requests(:)
      type(MPI_Status), intent(out), optional :: statuses(:)
      type(MPI_Request), allocatable :: mpi_requests(:)
      type(MPI_Status), allocatable :: temp_statuses(:)
      integer(int32) :: i, count, ierr

      count = size(requests)
      allocate (mpi_requests(count))

      do i = 1, count
         if (requests(i)%is_valid) then
            mpi_requests(i) = requests(i)%m_request
         else
            mpi_requests(i) = MPI_REQUEST_NULL
         end if
      end do

      if (present(statuses)) then
         call MPI_Waitall(count, mpi_requests, statuses, ierr)
      else
         allocate (temp_statuses(count))
         call MPI_Waitall(count, mpi_requests, temp_statuses, ierr)
      end if

      do i = 1, count
         call requests(i)%free()
      end do
   end subroutine request_waitall

   subroutine request_test(request, flag, status)
      type(request_t), intent(inout) :: request
      logical, intent(out) :: flag
      type(MPI_Status), intent(out), optional :: status
      integer(int32) :: ierr

      if (.not. request%is_valid) then
         flag = .true.
         return
      end if

      if (present(status)) then
         call MPI_Test(request%m_request, flag, status, ierr)
      else
         call MPI_Test(request%m_request, flag, MPI_STATUS_IGNORE, ierr)
      end if

      if (flag) then
         call request%free()
      end if
   end subroutine request_test

   ! ========================================================================
   ! Window creation
   ! ========================================================================

   !> Create MPI window for RMA operations
   !!
   !! Creates a window exposing local memory to remote RMA operations.
   !! Used for DDI distributed arrays.
   function create_win_dp_array(comm, base, win_size) result(win)
      type(comm_t), intent(in) :: comm
      real(dp), target :: base(:)
      integer(MPI_ADDRESS_KIND), intent(in) :: win_size
      type(win_t) :: win
      integer(int32) :: ierr
      integer(int32) :: disp_unit

      disp_unit = int(storage_size(base(1))/8_int32, int32)
      call MPI_Win_create(base, win_size, disp_unit, &
                          MPI_INFO_NULL, comm%get(), win%m_win, ierr)
      win%is_valid = .true.
   end function create_win_dp_array

   !> Create dynamic MPI window
   !!
   !! For windows where memory will be attached later.
   !! Useful for load balancing counters.
   function create_win_dynamic(comm) result(win)
      type(comm_t), intent(in) :: comm
      type(win_t) :: win
      integer(int32) :: ierr

      call MPI_Win_create_dynamic(MPI_INFO_NULL, comm%get(), win%m_win, ierr)
      win%is_valid = .true.
   end function create_win_dynamic

   !> Allocate window memory and create window in one call (1D array)
   !!
   !! This is more efficient than separate allocation + win_create.
   !! The baseptr is associated with the allocated memory.
   !! Memory is freed when window is finalized.
   subroutine create_win_allocate_dp_1d(comm, length, baseptr, win)
      use iso_c_binding, only: c_ptr, c_f_pointer
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: length
      real(dp), pointer, intent(out) :: baseptr(:)
      type(win_t), intent(out) :: win
      type(c_ptr) :: cptr
      integer(MPI_ADDRESS_KIND) :: win_size
      integer(int32) :: disp_unit, ierr

      win_size = int(length, MPI_ADDRESS_KIND)*int(storage_size(1.0_dp)/8_int32, MPI_ADDRESS_KIND)
      disp_unit = int(storage_size(1.0_dp)/8_int32, int32)

      call MPI_Win_allocate(win_size, disp_unit, MPI_INFO_NULL, &
                            comm%get(), cptr, win%m_win, ierr)
      call c_f_pointer(cptr, baseptr, [length])
      win%is_valid = .true.
   end subroutine create_win_allocate_dp_1d

   !> Allocate window memory and create window in one call (2D array)
   !!
   !! This is more efficient than separate allocation + win_create.
   !! The baseptr is associated with the allocated memory.
   !! Memory is freed when window is finalized.
   subroutine create_win_allocate_dp_2d(comm, dim1, dim2, baseptr, win)
      use iso_c_binding, only: c_ptr, c_f_pointer
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: dim1, dim2
      real(dp), pointer, intent(out) :: baseptr(:, :)
      type(win_t), intent(out) :: win
      type(c_ptr) :: cptr
      integer(MPI_ADDRESS_KIND) :: win_size
      integer(int32) :: disp_unit, ierr, total_size

      total_size = dim1*dim2
      win_size = int(total_size, MPI_ADDRESS_KIND)*int(storage_size(1.0_dp)/8_int32, MPI_ADDRESS_KIND)
      disp_unit = int(storage_size(1.0_dp)/8_int32, int32)

      call MPI_Win_allocate(win_size, disp_unit, MPI_INFO_NULL, &
                            comm%get(), cptr, win%m_win, ierr)
      call c_f_pointer(cptr, baseptr, [dim1, dim2])
      win%is_valid = .true.
   end subroutine create_win_allocate_dp_2d

   !> Allocate window memory for integer64 array
   !!
   !! Used for atomic counters in dynamic load balancing.
   subroutine create_win_allocate_i64_1d(comm, length, baseptr, win)
      use iso_c_binding, only: c_ptr, c_f_pointer
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: length
      integer(int64), pointer, intent(out) :: baseptr(:)
      type(win_t), intent(out) :: win
      type(c_ptr) :: cptr
      integer(MPI_ADDRESS_KIND) :: win_size
      integer(int32) :: disp_unit, ierr

      win_size = int(length, MPI_ADDRESS_KIND)*int(storage_size(1_int64)/8_int32, MPI_ADDRESS_KIND)
      disp_unit = int(storage_size(1_int64)/8_int32, int32)

      call MPI_Win_allocate(win_size, disp_unit, MPI_INFO_NULL, &
                            comm%get(), cptr, win%m_win, ierr)
      call c_f_pointer(cptr, baseptr, [length])
      win%is_valid = .true.
   end subroutine create_win_allocate_i64_1d

   !> Allocate window memory for single precision array
   subroutine create_win_allocate_sp_1d(comm, length, baseptr, win)
      use iso_c_binding, only: c_ptr, c_f_pointer
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: length
      real(sp), pointer, intent(out) :: baseptr(:)
      type(win_t), intent(out) :: win
      type(c_ptr) :: cptr
      integer(MPI_ADDRESS_KIND) :: win_size
      integer(int32) :: disp_unit, ierr

      win_size = int(length, MPI_ADDRESS_KIND)*int(storage_size(1.0_sp)/8_int32, MPI_ADDRESS_KIND)
      disp_unit = int(storage_size(1.0_sp)/8_int32, int32)

      call MPI_Win_allocate(win_size, disp_unit, MPI_INFO_NULL, &
                            comm%get(), cptr, win%m_win, ierr)
      call c_f_pointer(cptr, baseptr, [length])
      win%is_valid = .true.
   end subroutine create_win_allocate_sp_1d

   !> Allocate window memory for integer32 array
   subroutine create_win_allocate_i32_1d(comm, length, baseptr, win)
      use iso_c_binding, only: c_ptr, c_f_pointer
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: length
      integer(int32), pointer, intent(out) :: baseptr(:)
      type(win_t), intent(out) :: win
      type(c_ptr) :: cptr
      integer(MPI_ADDRESS_KIND) :: win_size
      integer(int32) :: disp_unit, ierr

      win_size = int(length, MPI_ADDRESS_KIND)*int(storage_size(1_int32)/8_int32, MPI_ADDRESS_KIND)
      disp_unit = int(storage_size(1_int32)/8_int32, int32)

      call MPI_Win_allocate(win_size, disp_unit, MPI_INFO_NULL, &
                            comm%get(), cptr, win%m_win, ierr)
      call c_f_pointer(cptr, baseptr, [length])
      win%is_valid = .true.
   end subroutine create_win_allocate_i32_1d

   ! ========================================================================
   ! Window query methods
   ! ========================================================================

   pure function win_is_null(this) result(is_null)
      class(win_t), intent(in) :: this
      logical :: is_null
      is_null = .not. this%is_valid
   end function win_is_null

   function win_get_handle(this) result(mpi_win_out)
      class(win_t), intent(in) :: this
      type(MPI_Win) :: mpi_win_out

      if (.not. this%is_valid) then
         error stop "Cannot get MPI_Win from null window"
      end if
      mpi_win_out = this%m_win
   end function win_get_handle

   ! ========================================================================
   ! Synchronization
   ! ========================================================================

   !> Fence synchronization for active target RMA
   !!
   !! Completes all pending RMA operations.
   !! Use before/after Get/Put/Accumulate operations.
   subroutine win_fence(this, assert)
      class(win_t), intent(in) :: this
      integer(int32), intent(in), optional :: assert
      integer(int32) :: ierr, assert_val

      if (present(assert)) then
         assert_val = assert
      else
         assert_val = 0_int32
      end if

      call MPI_Win_fence(assert_val, this%m_win, ierr)
   end subroutine win_fence

   !> Lock window for passive target RMA
   !!
   !! Begins RMA access epoch for specified target rank.
   !! Must be paired with unlock.
   subroutine win_lock(this, rank, lock_type)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: rank
      integer(int32), intent(in), optional :: lock_type
      integer(int32) :: ierr, ltype

      if (present(lock_type)) then
         ltype = lock_type
      else
         ltype = MPI_LOCK_SHARED
      end if

      call MPI_Win_lock(ltype, rank, 0_int32, this%m_win, ierr)
   end subroutine win_lock

   !> Unlock window for passive target RMA
   subroutine win_unlock(this, rank)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: rank
      integer(int32) :: ierr

      call MPI_Win_unlock(rank, this%m_win, ierr)
   end subroutine win_unlock

   !> Lock window on all ranks for passive target RMA
   !!
   !! Begins RMA access epoch for all ranks simultaneously.
   !! More efficient than individual locks when accessing multiple ranks.
   !! Must be paired with unlock_all.
   subroutine win_lock_all(this, assert)
      class(win_t), intent(in) :: this
      integer(int32), intent(in), optional :: assert
      integer(int32) :: ierr, assert_val

      if (present(assert)) then
         assert_val = assert
      else
         assert_val = 0_int32
      end if

      call MPI_Win_lock_all(assert_val, this%m_win, ierr)
   end subroutine win_lock_all

   !> Unlock window on all ranks for passive target RMA
   subroutine win_unlock_all(this)
      class(win_t), intent(in) :: this
      integer(int32) :: ierr

      call MPI_Win_unlock_all(this%m_win, ierr)
   end subroutine win_unlock_all

   !> Flush pending RMA operations to a specific rank
   !!
   !! Ensures all RMA operations to target rank have completed.
   !! Use within lock_all/unlock_all epoch.
   subroutine win_flush(this, rank)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: rank
      integer(int32) :: ierr

      call MPI_Win_flush(rank, this%m_win, ierr)
   end subroutine win_flush

   !> Flush pending RMA operations to all ranks
   !!
   !! Ensures all RMA operations to all ranks have completed.
   !! Use within lock_all/unlock_all epoch.
   subroutine win_flush_all(this)
      class(win_t), intent(in) :: this
      integer(int32) :: ierr

      call MPI_Win_flush_all(this%m_win, ierr)
   end subroutine win_flush_all

   ! ========================================================================
   ! RMA Get/Put/Accumulate operations
   ! ========================================================================

   !> Get data from remote window
   !!
   !! Retrieves data from target rank's window into local buffer.
   !! Must be called between fence or lock/unlock pairs.
   subroutine win_get_dp(this, target_rank, target_disp, count, buffer)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      real(dp), intent(out) :: buffer(*)
      integer(int32) :: ierr

      call MPI_Get(buffer, count, MPI_DOUBLE_PRECISION, &
                   target_rank, target_disp, count, MPI_DOUBLE_PRECISION, &
                   this%m_win, ierr)
   end subroutine win_get_dp

   !> Put data to remote window
   !!
   !! Sends data from local buffer to target rank's window.
   !! Must be called between fence or lock/unlock pairs.
   subroutine win_put_dp(this, target_rank, target_disp, count, buffer)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      real(dp), intent(in) :: buffer(*)
      integer(int32) :: ierr

      call MPI_Put(buffer, count, MPI_DOUBLE_PRECISION, &
                   target_rank, target_disp, count, MPI_DOUBLE_PRECISION, &
                   this%m_win, ierr)
   end subroutine win_put_dp

   !> Non-blocking get data from remote window
   !!
   !! Retrieves data from target rank's window into local buffer.
   !! Returns a request handle for later completion via wait().
   !! Must be called between lock_all/unlock_all pairs.
   subroutine win_rget_dp(this, target_rank, target_disp, count, buffer, request)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      real(dp), intent(out) :: buffer(*)
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      call MPI_Rget(buffer, count, MPI_DOUBLE_PRECISION, &
                    target_rank, target_disp, count, MPI_DOUBLE_PRECISION, &
                    this%m_win, request%m_request, ierr)
      request%is_valid = .true.
   end subroutine win_rget_dp

   !> Non-blocking put data to remote window
   !!
   !! Sends data from local buffer to target rank's window.
   !! Returns a request handle for later completion via wait().
   !! Must be called between lock_all/unlock_all pairs.
   subroutine win_rput_dp(this, target_rank, target_disp, count, buffer, request)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      real(dp), intent(in) :: buffer(*)
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      call MPI_Rput(buffer, count, MPI_DOUBLE_PRECISION, &
                    target_rank, target_disp, count, MPI_DOUBLE_PRECISION, &
                    this%m_win, request%m_request, ierr)
      request%is_valid = .true.
   end subroutine win_rput_dp

   !> Accumulate data to remote window
   !!
   !! Atomically adds local buffer to target rank's window.
   !! Critical for DDI_ACC (Fock matrix accumulation).
   subroutine win_accumulate_dp(this, target_rank, target_disp, count, buffer, op)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      real(dp), intent(in) :: buffer(*)
      type(MPI_Op), intent(in), optional :: op
      type(MPI_Op) :: mpi_op
      integer(int32) :: ierr

      if (present(op)) then
         mpi_op = op
      else
         mpi_op = MPI_SUM
      end if

      call MPI_Accumulate(buffer, count, MPI_DOUBLE_PRECISION, &
                          target_rank, target_disp, count, MPI_DOUBLE_PRECISION, &
                          mpi_op, this%m_win, ierr)
   end subroutine win_accumulate_dp

   ! ========================================================================
   ! Single precision (sp) RMA operations
   ! ========================================================================

   subroutine win_get_sp(this, target_rank, target_disp, count, buffer)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      real(sp), intent(out) :: buffer(*)
      integer(int32) :: ierr

      call MPI_Get(buffer, count, MPI_REAL, &
                   target_rank, target_disp, count, MPI_REAL, &
                   this%m_win, ierr)
   end subroutine win_get_sp

   subroutine win_put_sp(this, target_rank, target_disp, count, buffer)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      real(sp), intent(in) :: buffer(*)
      integer(int32) :: ierr

      call MPI_Put(buffer, count, MPI_REAL, &
                   target_rank, target_disp, count, MPI_REAL, &
                   this%m_win, ierr)
   end subroutine win_put_sp

   subroutine win_rget_sp(this, target_rank, target_disp, count, buffer, request)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      real(sp), intent(out) :: buffer(*)
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      call MPI_Rget(buffer, count, MPI_REAL, &
                    target_rank, target_disp, count, MPI_REAL, &
                    this%m_win, request%m_request, ierr)
      request%is_valid = .true.
   end subroutine win_rget_sp

   subroutine win_rput_sp(this, target_rank, target_disp, count, buffer, request)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      real(sp), intent(in) :: buffer(*)
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      call MPI_Rput(buffer, count, MPI_REAL, &
                    target_rank, target_disp, count, MPI_REAL, &
                    this%m_win, request%m_request, ierr)
      request%is_valid = .true.
   end subroutine win_rput_sp

   subroutine win_accumulate_sp(this, target_rank, target_disp, count, buffer, op)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      real(sp), intent(in) :: buffer(*)
      type(MPI_Op), intent(in), optional :: op
      type(MPI_Op) :: mpi_op
      integer(int32) :: ierr

      if (present(op)) then
         mpi_op = op
      else
         mpi_op = MPI_SUM
      end if

      call MPI_Accumulate(buffer, count, MPI_REAL, &
                          target_rank, target_disp, count, MPI_REAL, &
                          mpi_op, this%m_win, ierr)
   end subroutine win_accumulate_sp

   ! ========================================================================
   ! Integer32 (i32) RMA operations
   ! ========================================================================

   subroutine win_get_i32(this, target_rank, target_disp, count, buffer)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      integer(int32), intent(out) :: buffer(*)
      integer(int32) :: ierr

      call MPI_Get(buffer, count, MPI_INTEGER, &
                   target_rank, target_disp, count, MPI_INTEGER, &
                   this%m_win, ierr)
   end subroutine win_get_i32

   subroutine win_put_i32(this, target_rank, target_disp, count, buffer)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      integer(int32), intent(in) :: buffer(*)
      integer(int32) :: ierr

      call MPI_Put(buffer, count, MPI_INTEGER, &
                   target_rank, target_disp, count, MPI_INTEGER, &
                   this%m_win, ierr)
   end subroutine win_put_i32

   subroutine win_rget_i32(this, target_rank, target_disp, count, buffer, request)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      integer(int32), intent(out) :: buffer(*)
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      call MPI_Rget(buffer, count, MPI_INTEGER, &
                    target_rank, target_disp, count, MPI_INTEGER, &
                    this%m_win, request%m_request, ierr)
      request%is_valid = .true.
   end subroutine win_rget_i32

   subroutine win_rput_i32(this, target_rank, target_disp, count, buffer, request)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      integer(int32), intent(in) :: buffer(*)
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      call MPI_Rput(buffer, count, MPI_INTEGER, &
                    target_rank, target_disp, count, MPI_INTEGER, &
                    this%m_win, request%m_request, ierr)
      request%is_valid = .true.
   end subroutine win_rput_i32

   subroutine win_accumulate_i32(this, target_rank, target_disp, count, buffer, op)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      integer(int32), intent(in) :: buffer(*)
      type(MPI_Op), intent(in), optional :: op
      type(MPI_Op) :: mpi_op
      integer(int32) :: ierr

      if (present(op)) then
         mpi_op = op
      else
         mpi_op = MPI_SUM
      end if

      call MPI_Accumulate(buffer, count, MPI_INTEGER, &
                          target_rank, target_disp, count, MPI_INTEGER, &
                          mpi_op, this%m_win, ierr)
   end subroutine win_accumulate_i32

   ! ========================================================================
   ! Integer64 (i64) RMA operations
   ! ========================================================================

   subroutine win_get_i64(this, target_rank, target_disp, count, buffer)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      integer(int64), intent(out) :: buffer(*)
      integer(int32) :: ierr

      call MPI_Get(buffer, count, MPI_INTEGER8, &
                   target_rank, target_disp, count, MPI_INTEGER8, &
                   this%m_win, ierr)
   end subroutine win_get_i64

   subroutine win_put_i64(this, target_rank, target_disp, count, buffer)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      integer(int64), intent(in) :: buffer(*)
      integer(int32) :: ierr

      call MPI_Put(buffer, count, MPI_INTEGER8, &
                   target_rank, target_disp, count, MPI_INTEGER8, &
                   this%m_win, ierr)
   end subroutine win_put_i64

   subroutine win_rget_i64(this, target_rank, target_disp, count, buffer, request)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      integer(int64), intent(out) :: buffer(*)
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      call MPI_Rget(buffer, count, MPI_INTEGER8, &
                    target_rank, target_disp, count, MPI_INTEGER8, &
                    this%m_win, request%m_request, ierr)
      request%is_valid = .true.
   end subroutine win_rget_i64

   subroutine win_rput_i64(this, target_rank, target_disp, count, buffer, request)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      integer(int64), intent(in) :: buffer(*)
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      call MPI_Rput(buffer, count, MPI_INTEGER8, &
                    target_rank, target_disp, count, MPI_INTEGER8, &
                    this%m_win, request%m_request, ierr)
      request%is_valid = .true.
   end subroutine win_rput_i64

   subroutine win_accumulate_i64(this, target_rank, target_disp, count, buffer, op)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      integer(int64), intent(in) :: buffer(*)
      type(MPI_Op), intent(in), optional :: op
      type(MPI_Op) :: mpi_op
      integer(int32) :: ierr

      if (present(op)) then
         mpi_op = op
      else
         mpi_op = MPI_SUM
      end if

      call MPI_Accumulate(buffer, count, MPI_INTEGER8, &
                          target_rank, target_disp, count, MPI_INTEGER8, &
                          mpi_op, this%m_win, ierr)
   end subroutine win_accumulate_i64

   !> Atomic fetch-and-add for load balancing
   !!
   !! Atomically increments remote counter and returns old value.
   !! Used for DDI_DLBNEXT (dynamic load balancing).
   subroutine win_fetch_and_add_i64(this, target_rank, target_disp, value, result)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int64), intent(in) :: value
      integer(int64), intent(out) :: result
      integer(int32) :: ierr

      call MPI_Fetch_and_op(value, result, MPI_INTEGER8, &
                            target_rank, target_disp, MPI_SUM, this%m_win, ierr)
   end subroutine win_fetch_and_add_i64

   ! ========================================================================
   ! Window cleanup
   ! ========================================================================

   subroutine win_finalize(this)
      class(win_t), intent(inout) :: this
      integer(int32) :: ierr

      if (this%is_valid .and. this%m_win /= MPI_WIN_NULL) then
         call MPI_Win_free(this%m_win, ierr)
         this%is_valid = .false.
         this%m_win = MPI_WIN_NULL
      end if
   end subroutine win_finalize

   ! ========================================================================
   ! Allreduce operations (for DDI_GSUMF/GSUMI)
   ! ========================================================================

   !> Allreduce for scalar double precision
   !!
   !! In-place global reduction. Replaces DDI_GSUMF for scalars.
   subroutine allreduce_dp(comm, buffer, op)
      type(comm_t), intent(in) :: comm
      real(dp), intent(inout) :: buffer
      type(MPI_Op), intent(in), optional :: op
      type(MPI_Op) :: mpi_op
      integer :: ierr

      if (present(op)) then
         mpi_op = op
      else
         mpi_op = MPI_SUM
      end if

      call MPI_Allreduce(MPI_IN_PLACE, buffer, 1, MPI_DOUBLE_PRECISION, &
                         mpi_op, comm%get(), ierr)
   end subroutine allreduce_dp

   !> Allreduce for double precision array
   !!
   !! In-place global reduction. Replaces DDI_GSUMF for arrays.
   !! This is THE most-called DDI function (1,301 calls in GAMESS).
   subroutine allreduce_dp_array(comm, buffer, count, op)
      type(comm_t), intent(in) :: comm
      real(dp), intent(inout) :: buffer(:)
      integer, intent(in), optional :: count
      type(MPI_Op), intent(in), optional :: op
      type(MPI_Op) :: mpi_op
      integer :: ierr, n

      if (present(count)) then
         n = count
      else
         n = size(buffer)
      end if

      if (present(op)) then
         mpi_op = op
      else
         mpi_op = MPI_SUM
      end if

      call MPI_Allreduce(MPI_IN_PLACE, buffer, n, MPI_DOUBLE_PRECISION, &
                         mpi_op, comm%get(), ierr)
   end subroutine allreduce_dp_array

   !> Allreduce for scalar integer
   !!
   !! In-place global reduction. Replaces DDI_GSUMI for scalars.
   subroutine allreduce_i32(comm, buffer, op)
      type(comm_t), intent(in) :: comm
      integer(int32), intent(inout) :: buffer
      type(MPI_Op), intent(in), optional :: op
      type(MPI_Op) :: mpi_op
      integer :: ierr

      if (present(op)) then
         mpi_op = op
      else
         mpi_op = MPI_SUM
      end if

      call MPI_Allreduce(MPI_IN_PLACE, buffer, 1, MPI_INTEGER, &
                         mpi_op, comm%get(), ierr)
   end subroutine allreduce_i32

   !> Allreduce for integer array
   !!
   !! In-place global reduction. Replaces DDI_GSUMI for arrays.
   subroutine allreduce_i32_array(comm, buffer, count, op)
      type(comm_t), intent(in) :: comm
      integer(int32), intent(inout) :: buffer(:)
      integer, intent(in), optional :: count
      type(MPI_Op), intent(in), optional :: op
      type(MPI_Op) :: mpi_op
      integer :: ierr, n

      if (present(count)) then
         n = count
      else
         n = size(buffer)
      end if

      if (present(op)) then
         mpi_op = op
      else
         mpi_op = MPI_SUM
      end if

      call MPI_Allreduce(MPI_IN_PLACE, buffer, n, MPI_INTEGER, &
                         mpi_op, comm%get(), ierr)
   end subroutine allreduce_i32_array

   !> Non-in-place allreduce for scalar double precision
   !!
   !! Reduces sendbuf and stores result in recvbuf.
   !! Useful for timestep reduction where local value must be preserved.
   subroutine allreduce_dp_to(comm, sendbuf, recvbuf, op)
      type(comm_t), intent(in) :: comm
      real(dp), intent(in) :: sendbuf
      real(dp), intent(out) :: recvbuf
      type(MPI_Op), intent(in), optional :: op
      type(MPI_Op) :: mpi_op
      integer :: ierr

      if (present(op)) then
         mpi_op = op
      else
         mpi_op = MPI_SUM
      end if

      call MPI_Allreduce(sendbuf, recvbuf, 1, MPI_DOUBLE_PRECISION, &
                         mpi_op, comm%get(), ierr)
   end subroutine allreduce_dp_to

   !> Non-in-place allreduce for double precision array
   !!
   !! Reduces sendbuf and stores result in recvbuf.
   subroutine allreduce_dp_array_to(comm, sendbuf, recvbuf, count, op)
      type(comm_t), intent(in) :: comm
      real(dp), intent(in) :: sendbuf(:)
      real(dp), intent(out) :: recvbuf(:)
      integer, intent(in), optional :: count
      type(MPI_Op), intent(in), optional :: op
      type(MPI_Op) :: mpi_op
      integer :: ierr, n

      if (present(count)) then
         n = count
      else
         n = size(sendbuf)
      end if

      if (present(op)) then
         mpi_op = op
      else
         mpi_op = MPI_SUM
      end if

      call MPI_Allreduce(sendbuf, recvbuf, n, MPI_DOUBLE_PRECISION, &
                         mpi_op, comm%get(), ierr)
   end subroutine allreduce_dp_array_to

end module pic_mpi_f08
