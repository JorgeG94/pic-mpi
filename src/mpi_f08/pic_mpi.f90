!> Modern MPI wrapper module using mpi_f08 interface
!!
!! This module provides a high-level object-oriented interface to MPI
!! using the modern mpi_f08 bindings. It wraps MPI communicators and
!! requests into derived types with type-bound procedures.
!!
module pic_mpi_f08
   use pic_types, only: int32, dp
   use mpi_f08, only: MPI_Comm, MPI_Status, MPI_Request, MPI_COMM_NULL, MPI_COMM_WORLD, &
                      MPI_COMM_TYPE_SHARED, MPI_INFO_NULL, MPI_UNDEFINED, &
                      MPI_INTEGER, MPI_STATUS_IGNORE, MPI_REQUEST_NULL, MPI_Comm_rank, &
                      MPI_Comm_size, MPI_Comm_dup, MPI_Barrier, &
                      MPI_Comm_split_type, MPI_Comm_split, MPI_Send, MPI_Recv, &
                      MPI_Isend, MPI_Irecv, MPI_Wait, MPI_Waitall, MPI_Test, &
                      MPI_Probe, MPI_Get_count, MPI_Iprobe, MPI_Comm_free, &
                      MPI_Abort, MPI_Allgather, MPI_Get_processor_name, &
                      MPI_Bcast, MPI_Init, MPI_Init_thread, MPI_Query_thread, MPI_Finalize, &
                      MPI_THREAD_SINGLE, MPI_THREAD_FUNNELED, MPI_THREAD_SERIALIZED, MPI_THREAD_MULTIPLE, &
                      MPI_ANY_SOURCE, MPI_ANY_TAG, &
                      MPI_MAX_PROCESSOR_NAME, &
                      operator(==), operator(/=), MPI_DOUBLE_PRECISION
   implicit none
   private

   public :: comm_t, comm_world, comm_null
   public :: send, recv, isend, irecv
   public :: request_t, wait, waitall, test
   public :: iprobe, abort_comm, allgather, get_processor_name, bcast
   public :: pic_mpi_init, pic_mpi_finalize, pic_mpi_query_thread_level

   ! Export MPI types and constants needed by applications
   public :: MPI_Status, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_MAX_PROCESSOR_NAME
   public :: MPI_THREAD_SINGLE, MPI_THREAD_FUNNELED, MPI_THREAD_SERIALIZED, MPI_THREAD_MULTIPLE

   !> Request type for non-blocking MPI operations
   !!
   !! Wraps MPI_Request to provide object-oriented interface for
   !! non-blocking communication operations (isend, irecv)
   type :: request_t
      private
      type(MPI_Request) :: m_request = MPI_REQUEST_NULL !! Internal MPI request handle
      logical :: is_valid = .false. !! Validity flag
   contains
      procedure :: is_null => request_is_null !! Check if request is null
      procedure :: get => request_get !! Get underlying MPI_Request
      procedure :: free => request_free !! Free the request
   end type request_t

   !> MPI communicator wrapper type
   !!
   !! Provides object-oriented interface to MPI communicators with
   !! type-bound procedures for common operations. Automatically caches
   !! rank and size information for efficient access.
   type :: comm_t
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
      module procedure :: comm_send_real_dp
      module procedure :: comm_send_real_dp_array
   end interface send

   interface recv
      module procedure :: comm_recv_integer
      module procedure :: comm_recv_integer_array
      module procedure :: comm_recv_real_dp
      module procedure :: comm_recv_real_dp_array
   end interface recv

   interface iprobe
      module procedure :: comm_iprobe
   end interface iprobe

   interface allgather
      module procedure :: comm_allgather_integer
   end interface allgather

   interface bcast
      module procedure :: comm_bcast_integer
   end interface bcast

   interface isend
      module procedure :: comm_isend_integer
      module procedure :: comm_isend_integer_array
      module procedure :: comm_isend_real_dp
      module procedure :: comm_isend_real_dp_array
   end interface isend

   interface irecv
      module procedure :: comm_irecv_integer
      module procedure :: comm_irecv_integer_array
      module procedure :: comm_irecv_real_dp
      module procedure :: comm_irecv_real_dp_array
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

contains

   !> Create comm_t from MPI_Comm handle
   !!
   !! Internal helper function that wraps an MPI_Comm into a comm_t object
   !! and caches rank and size information
   function create_comm_from_mpi(mpi_comm_in) result(comm)
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

   !> Create a duplicate of MPI_COMM_WORLD
   !!
   !! Creates a new communicator that duplicates MPI_COMM_WORLD.
   !! This is the standard way to obtain a communicator for application use.
   function create_world_comm() result(comm)
      type(comm_t) :: comm
      type(MPI_Comm) :: dup_comm
      integer(int32) :: ierr

      call MPI_Comm_dup(MPI_COMM_WORLD, dup_comm, ierr)
      comm = create_comm_from_mpi(dup_comm)

   end function create_world_comm

   !> Create a null communicator
   !!
   !! Creates an invalid/null communicator object that can be used
   !! for initialization or to represent absence of a communicator.
   function create_null_comm() result(comm)
      type(comm_t) :: comm

      ! Explicitly initialize to null/invalid state
      comm%m_comm = MPI_COMM_NULL
      comm%m_rank = -1
      comm%m_size = -1
      comm%is_valid = .false.
   end function create_null_comm

   !> Get rank of calling process in communicator
   !!
   !! Returns the 0-indexed rank of the calling process
   pure function comm_rank(this) result(rank)
      class(comm_t), intent(in) :: this
      integer :: rank
      rank = this%m_rank
   end function comm_rank

   !> Get size of communicator
   !!
   !! Returns the number of processes in the communicator
   pure function m_size_func(this) result(size)
      class(comm_t), intent(in) :: this
      integer :: size
      size = this%m_size
   end function m_size_func

   !> Check if calling process is the leader
   !!
   !! Returns true if the calling process has rank 0
   pure function comm_leader(this) result(is_leader)
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

   !> Synchronization barrier
   !!
   !! Blocks until all processes in the communicator have called barrier
   subroutine comm_barrier(this)
      class(comm_t), intent(in) :: this
      integer(int32) :: ierr
      call MPI_Barrier(this%m_comm, ierr)
   end subroutine comm_barrier

   !> Split communicator into shared memory communicators
   !!
   !! Creates a new communicator containing only processes that share
   !! memory with each other (typically processes on the same node)
   function comm_split_shared(this) result(new_comm)
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

   !> Split communicator by color
   !!
   !! Partitions the communicator into disjoint subgroups based on color.
   !! Processes with the same color end up in the same new communicator.
   function comm_split_by_color(this, color) result(new_comm)
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

   !> Send a single integer value
   !!
   !! Blocking send of an integer to specified destination
   subroutine comm_send_integer(comm, data, dest, tag)
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: data
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr

      call MPI_Send(data, 1, MPI_INTEGER, dest, tag, comm%m_comm, ierr)
   end subroutine comm_send_integer

   !> Send an integer array
   !!
   !! Blocking send of an integer array to specified destination
   subroutine comm_send_integer_array(comm, data, dest, tag)
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: data(:)
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr

      call MPI_Send(data, size(data), MPI_INTEGER, dest, tag, comm%m_comm, ierr)
   end subroutine comm_send_integer_array

   subroutine comm_send_real_dp(comm, data, dest, tag)
      type(comm_t), intent(in) :: comm
      real(dp), intent(in) :: data
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr

      call MPI_Send(data, 1, MPI_DOUBLE_PRECISION, dest, tag, comm%m_comm, ierr)
   end subroutine comm_send_real_dp

   subroutine comm_send_real_dp_array(comm, data, dest, tag)
      type(comm_t), intent(in) :: comm
      real(dp), intent(in) :: data(:)
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr

      call MPI_Send(data, size(data), MPI_DOUBLE_PRECISION, dest, tag, comm%m_comm, ierr)
   end subroutine comm_send_real_dp_array

   !> Receive a single integer value
   !!
   !! Blocking receive of an integer from specified source.
   !! Use MPI_ANY_SOURCE or MPI_ANY_TAG for wildcards.
   subroutine comm_recv_integer(comm, data, source, tag, status)
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

   subroutine comm_recv_real_dp(comm, data, source, tag, status)
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

   subroutine comm_iprobe(comm, source, tag, message_pending, status)
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      logical, intent(out) :: message_pending
      type(MPI_Status), intent(out) :: status
      integer(int32) :: ierr

      call MPI_Iprobe(source, tag, comm%m_comm, message_pending, status, ierr)
   end subroutine comm_iprobe

   subroutine comm_finalize(this)
      class(comm_t), intent(inout) :: this
      integer(int32) :: ierr

      if (this%is_valid .and. this%m_comm /= MPI_COMM_NULL) then
         call MPI_Comm_free(this%m_comm, ierr)
         this%is_valid = .false.
         this%m_comm = MPI_COMM_NULL
      end if
   end subroutine comm_finalize

   subroutine abort_comm(comm, errorcode)
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: errorcode
      integer(int32) :: ierr

      call MPI_Abort(comm%m_comm, errorcode, ierr)
   end subroutine abort_comm

   subroutine comm_allgather_integer(comm, sendbuf, recvbuf)
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: sendbuf
      integer(int32), intent(out) :: recvbuf(:)
      integer(int32) :: ierr

      call MPI_Allgather(sendbuf, 1, MPI_INTEGER, recvbuf, 1, MPI_INTEGER, comm%m_comm, ierr)
   end subroutine comm_allgather_integer

   subroutine comm_bcast_integer(comm, buffer, count, root)
      type(comm_t), intent(in) :: comm
      integer(int32), intent(inout) :: buffer
      integer(int32), intent(in) :: count
      integer(int32), intent(in) :: root
      integer(int32) :: ierr

      call MPI_Bcast(buffer, count, MPI_INTEGER, root, comm%m_comm, ierr)
   end subroutine comm_bcast_integer

   subroutine get_processor_name(name, namelen)
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
      class(request_t), intent(in) :: this
      logical :: is_null
      is_null = .not. this%is_valid
   end function request_is_null

   function request_get(this) result(mpi_request_out)
      class(request_t), intent(in) :: this
      type(MPI_Request) :: mpi_request_out

      if (.not. this%is_valid) then
         error stop "Cannot get MPI_Request from null request"
      end if
      mpi_request_out = this%m_request
   end function request_get

   subroutine request_free(this)
      class(request_t), intent(inout) :: this
      this%m_request = MPI_REQUEST_NULL
      this%is_valid = .false.
   end subroutine request_free

   ! ========================================================================
   ! Non-blocking send operations
   ! ========================================================================

   !> Non-blocking send of a single integer
   !!
   !! Initiates a non-blocking send operation. The request must be
   !! waited on using wait() or test() before the buffer can be reused.
   subroutine comm_isend_integer(comm, data, dest, tag, request)
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
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: data(:)
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      call MPI_Isend(data, size(data), MPI_INTEGER, dest, tag, comm%m_comm, request%m_request, ierr)
      request%is_valid = .true.
   end subroutine comm_isend_integer_array

   subroutine comm_isend_real_dp(comm, data, dest, tag, request)
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
      type(comm_t), intent(in) :: comm
      real(dp), intent(in) :: data(:)
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      call MPI_Isend(data, size(data), MPI_DOUBLE_PRECISION, dest, tag, comm%m_comm, request%m_request, ierr)
      request%is_valid = .true.
   end subroutine comm_isend_real_dp_array

   ! ========================================================================
   ! Non-blocking receive operations
   ! ========================================================================

   subroutine comm_irecv_integer(comm, data, source, tag, request)
      type(comm_t), intent(in) :: comm
      integer(int32), intent(out) :: data
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      call MPI_Irecv(data, 1, MPI_INTEGER, source, tag, comm%m_comm, request%m_request, ierr)
      request%is_valid = .true.
   end subroutine comm_irecv_integer

   subroutine comm_irecv_integer_array(comm, data, count, source, tag, request)
      type(comm_t), intent(in) :: comm
      integer(int32), intent(out) :: data(:)
      integer(int32), intent(in) :: count
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      call MPI_Irecv(data, count, MPI_INTEGER, source, tag, comm%m_comm, request%m_request, ierr)
      request%is_valid = .true.
   end subroutine comm_irecv_integer_array

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

   subroutine comm_irecv_real_dp_array(comm, data, count, source, tag, request)
      type(comm_t), intent(in) :: comm
      real(dp), intent(out) :: data(:)
      integer(int32), intent(in) :: count
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      call MPI_Irecv(data, count, MPI_DOUBLE_PRECISION, source, tag, comm%m_comm, request%m_request, ierr)
      request%is_valid = .true.
   end subroutine comm_irecv_real_dp_array

   ! ========================================================================
   ! Request completion operations
   ! ========================================================================

   !> Wait for a non-blocking operation to complete
   !!
   !! Blocks until the operation associated with the request completes.
   !! The request is freed after completion.
   subroutine request_wait(request, status)
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

   !> Wait for all non-blocking operations to complete
   !!
   !! Blocks until all operations in the request array complete.
   !! All requests are freed after completion.
   subroutine request_waitall(requests, statuses)
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

end module pic_mpi_f08
