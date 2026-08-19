!> Single-rank stand-in for `pic_mpi_f08`, for builds without MPI.
!!
!! WHY THIS EXISTS
!! ---------------
!! Not every consumer of PIC wants an MPI dependency. A CI runner without one,
!! a laptop build, or a compiler whose vendor MPI is awkward to get hold of --
!! in each case the alternative is no build at all, because `pic_mpi_lib` had
!! only the two MPI backends to choose from.
!!
!! This is the third: the same public interface, with one rank and no library
!! underneath. `pic_mpi_lib` selects it with `USE_SERIAL`, exactly as it
!! selects the legacy bindings with `USE_LEGACY`, so a consumer's source does
!! not change.
!!
!! WHAT THE BODIES DO
!! ------------------
!!   * `rank` is 0, `size` is 1, and every rank is the leader.
!!   * A broadcast is a no-op: the data is already on the only rank there is.
!!   * A gather or a reduction copies the contribution into the result.
!!   * `iprobe` reports nothing pending, because no peer exists to have sent
!!     anything. That keeps polling loops terminating rather than spinning.
!!   * Waiting completes immediately: nothing was ever started.
!!   * POINT-TO-POINT AND ONE-SIDED CALLS STOP THE PROGRAM. A single rank has
!!     no peer, so reaching a send or a get means the caller failed to take
!!     its `size() == 1` path. Returning quietly would turn that into wrong
!!     results much later; stopping names it where it happened.
!!
!! GENERATED, and deliberately. The signatures are taken from
!! `mpi_f08/pic_mpi.f90` so that the two cannot disagree -- 140 procedures
!! behind a dozen generics is 140 chances to typo an interface that would only
!! surface as an unresolved generic in somebody else's build. Only the bodies
!! are written by hand, by family. Regenerate with `tools/gen_serial.py`.
!!
module pic_mpi_serial
   use pic_types, only: int32, int64, sp, dp
   implicit none
   private

   public :: comm_t, comm_world, comm_null
   public :: send, recv, isend, irecv
   public :: comm_isend_real_dp_array_n, comm_irecv_real_dp_array_n  ! Direct export for host_data blocks (nvhpc bug workaround)
   public :: comm_isend_real_sp_array_n, comm_irecv_real_sp_array_n  ! Single precision equivalents
   public :: comm_recv_real_dp_array_n, comm_recv_real_sp_array_n    ! Fixed-size blocking recv variants
   public :: comm_recv_integer_array_n, comm_recv_integer64_array_n
   public :: comm_send_real_dp_array_2d_n, comm_recv_real_dp_array_2d_n
   public :: comm_isend_real_dp_array_2d_n
   public :: comm_send_integer_array_2d_n, comm_recv_integer_array_2d_n
   public :: request_t, wait, waitall, test
   public :: iprobe, probe, abort_comm, allgather, get_processor_name, bcast
   public :: pic_mpi_init, pic_mpi_finalize, pic_mpi_query_thread_level
   public :: win_t, win_create, win_create_dynamic, win_allocate
   public :: allreduce
   public :: MPI_Status, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_MAX_PROCESSOR_NAME, MPI_Request
   public :: MPI_THREAD_SINGLE, MPI_THREAD_FUNNELED, MPI_THREAD_SERIALIZED, MPI_THREAD_MULTIPLE
   public :: MPI_SUM, MPI_MIN, MPI_MAX

   !
   ! Stand-ins for the handles and constants `mpi_f08` would supply. They carry
   ! no meaning here beyond letting the same declarations compile: a consumer
   ! that stores an `MPI_Status` or compares against `MPI_ANY_SOURCE` keeps
   ! working, it just never receives anything.
   !
   type :: MPI_Comm
      integer(int32) :: v = 0
   end type MPI_Comm

   type :: MPI_Request
      integer(int32) :: v = 0
   end type MPI_Request

   type :: MPI_Status
      integer(int32) :: MPI_SOURCE = -1
      integer(int32) :: MPI_TAG = -1
      integer(int32) :: MPI_ERROR = 0
   end type MPI_Status

   type :: MPI_Op
      integer(int32) :: v = 0
   end type MPI_Op

   type :: MPI_Win
      integer(int32) :: v = 0
   end type MPI_Win

   type(MPI_Comm), parameter :: MPI_COMM_NULL = MPI_Comm(0)
   type(MPI_Op), parameter :: MPI_SUM = MPI_Op(1)
   type(MPI_Op), parameter :: MPI_MIN = MPI_Op(2)
   type(MPI_Op), parameter :: MPI_MAX = MPI_Op(3)

   integer(int32), parameter :: MPI_ANY_SOURCE = -1
   integer(int32), parameter :: MPI_ANY_TAG = -1
   integer(int32), parameter :: MPI_MAX_PROCESSOR_NAME = 256
   integer(int32), parameter :: MPI_THREAD_SINGLE = 0
   integer(int32), parameter :: MPI_THREAD_FUNNELED = 1
   integer(int32), parameter :: MPI_THREAD_SERIALIZED = 2
   integer(int32), parameter :: MPI_THREAD_MULTIPLE = 3
   integer, parameter :: MPI_ADDRESS_KIND = int64

   type(MPI_Request), parameter :: MPI_REQUEST_NULL = MPI_Request(0)
   type(MPI_Win), parameter :: MPI_WIN_NULL = MPI_Win(0)

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

   type :: win_t
      !! MPI-3 Window type for one-sided communication (RMA)
      !!
      !! Wraps MPI_Win to provide object-oriented interface for
      !! Remote Memory Access (RMA) operations needed for DDI
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

   type :: comm_t
      !! MPI communicator wrapper type
      !!
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
      module procedure :: comm_send_real_dp_array_2d_n
      module procedure :: comm_send_integer_array_2d_n
      module procedure :: comm_send_real_sp
      module procedure :: comm_send_real_sp_array
      module procedure :: comm_send_real_sp_array_2d
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
      module procedure :: comm_recv_real_sp
      module procedure :: comm_recv_real_sp_array
      module procedure :: comm_recv_real_sp_array_2d
      module procedure :: comm_recv_logical
      ! Fixed-size receive variants (count is positional, no probe-then-allocate).
      ! Use these when receiving into a pre-allocated buffer.
      module procedure :: comm_recv_integer_array_n
      module procedure :: comm_recv_integer64_array_n
      module procedure :: comm_recv_real_dp_array_n
      module procedure :: comm_recv_real_sp_array_n
      module procedure :: comm_recv_real_dp_array_2d_n
      module procedure :: comm_recv_integer_array_2d_n
   end interface recv

   interface iprobe
      module procedure :: comm_iprobe
   end interface iprobe

   interface probe
      module procedure :: comm_probe
   end interface probe

   interface allgather
      module procedure :: comm_allgather_integer
   end interface allgather

   interface bcast
      module procedure :: comm_bcast_integer
      module procedure :: comm_bcast_integer64
      module procedure :: comm_bcast_real_dp
      module procedure :: comm_bcast_real_dp_array
      module procedure :: comm_bcast_real_sp
      module procedure :: comm_bcast_real_sp_array
   end interface bcast

   interface isend
      module procedure :: comm_isend_integer
      module procedure :: comm_isend_integer_array
      module procedure :: comm_isend_integer64
      module procedure :: comm_isend_integer64_array
      module procedure :: comm_isend_real_dp
      module procedure :: comm_isend_real_dp_array
      module procedure :: comm_isend_real_dp_array_2d
      module procedure :: comm_isend_real_dp_array_2d_n
      module procedure :: comm_isend_real_sp
      module procedure :: comm_isend_real_sp_array
      module procedure :: comm_isend_real_sp_array_2d
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
      module procedure :: comm_irecv_real_sp
      module procedure :: comm_irecv_real_sp_array
      module procedure :: comm_irecv_real_sp_array_2d
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
      module procedure :: allreduce_sp
      module procedure :: allreduce_sp_array
      module procedure :: allreduce_i32
      module procedure :: allreduce_i32_array
      module procedure :: allreduce_dp_to
      module procedure :: allreduce_dp_array_to
      module procedure :: allreduce_sp_to
      module procedure :: allreduce_sp_array_to
   end interface allreduce


contains

   function create_comm_from_mpi(mpi_comm_in) result(comm)
   !! Internal helper function that wraps an MPI_Comm into a comm_t object
   !! and caches rank and size information
      type(MPI_Comm), intent(in) :: mpi_comm_in
      type(comm_t) :: comm
      integer(int32) :: ierr

      ! The only communicator there is.
      comm%m_rank = 0
      comm%m_size = 1
      comm%is_valid = .true.
   end function create_comm_from_mpi

   function create_world_comm() result(comm)
   !! Creates a new communicator that duplicates MPI_COMM_WORLD.
   !! This is the standard way to obtain a communicator for application use.
      type(comm_t) :: comm
      type(MPI_Comm) :: dup_comm
      integer(int32) :: ierr

      ! The only communicator there is.
      comm%m_rank = 0
      comm%m_size = 1
      comm%is_valid = .true.
   end function create_world_comm

   function create_null_comm() result(comm)
   !! Creates an invalid/null communicator object that can be used
   !! for initialization or to represent absence of a communicator.
      type(comm_t) :: comm

      ! Explicitly initialize to null/invalid state
      ! The only communicator there is.
      comm%m_rank = 0
      comm%m_size = 1
      comm%is_valid = .false.
   end function create_null_comm

   pure function comm_rank(this) result(rank)
   !!
   !! Returns the 0-indexed rank of the calling process
      class(comm_t), intent(in) :: this
      integer :: rank
      rank = 0
   end function comm_rank

   pure function m_size_func(this) result(size)
   !! Returns the number of processes in the communicator
      class(comm_t), intent(in) :: this
      integer :: size
      size = 1
   end function m_size_func

   pure function comm_leader(this) result(is_leader)
   !! Returns true if the calling process has rank 0
      class(comm_t), intent(in) :: this
      logical :: is_leader
      is_leader = .true.
   end function comm_leader

   pure function comm_is_null(this) result(is_null)
      class(comm_t), intent(in) :: this
      logical :: is_null
      is_null = .not. this%is_valid
   end function comm_is_null

   function comm_get(this) result(mpi_comm_out)
      class(comm_t), intent(in) :: this
      type(MPI_Comm) :: mpi_comm_out

      mpi_comm_out = this%m_comm
   end function comm_get

   subroutine comm_barrier(this)
   !! Blocks until all processes in the communicator have called barrier
      class(comm_t), intent(in) :: this
      integer(int32) :: ierr
      ! Nothing to synchronise or release on one rank.
   end subroutine comm_barrier

   function comm_split_shared(this) result(new_comm)
   !! Creates a new communicator containing only processes that share
   !! memory with each other (typically processes on the same node)
      class(comm_t), intent(in) :: this
      type(comm_t) :: new_comm
      type(MPI_Comm) :: mpi_comm_new
      integer(int32) :: ierr

      ! Every split of one rank is that rank.
      new_comm = this
   end function comm_split_shared

   function comm_split_by_color(this, color) result(new_comm)
   !! Partitions the communicator into disjoint subgroups based on color.
   !! Processes with the same color end up in the same new communicator.
      class(comm_t), intent(in) :: this
      integer, intent(in) :: color
      type(comm_t) :: new_comm
      type(MPI_Comm) :: mpi_comm_new
      integer(int32) :: ierr

      ! Every split of one rank is that rank.
      new_comm = this
   end function comm_split_by_color

   function comm_discard_leader(this) result(new_comm)
      class(comm_t), intent(in) :: this
      type(comm_t) :: new_comm
      integer :: color

      ! Every split of one rank is that rank.
      new_comm = this
   end function comm_discard_leader

   function comm_discard_to(this, num_ranks) result(new_comm)
      class(comm_t), intent(in) :: this
      integer, intent(in) :: num_ranks
      type(comm_t) :: new_comm
      integer :: color

      ! Every split of one rank is that rank.
      new_comm = this
   end function comm_discard_to

   function comm_duplicate(this) result(new_comm)
      class(comm_t), intent(in) :: this
      type(comm_t) :: new_comm
      type(MPI_Comm) :: mpi_comm_new
      integer(int32) :: ierr

      ! Every split of one rank is that rank.
      new_comm = this
   end function comm_duplicate

   subroutine comm_send_integer(comm, data, dest, tag)
   !! Blocking send of an integer to specified destination
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: data
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_send_integer needs a peer rank"
   end subroutine comm_send_integer

   subroutine comm_send_integer_array(comm, data, dest, tag)
   !! Blocking send of an integer array to specified destination
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: data(:)
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_send_integer_array needs a peer rank"
   end subroutine comm_send_integer_array

   subroutine comm_send_integer64(comm, data, dest, tag)
   !! Blocking send of an integer64 to specified destination
      type(comm_t), intent(in) :: comm
      integer(int64), intent(in) :: data
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_send_integer64 needs a peer rank"
   end subroutine comm_send_integer64

   subroutine comm_send_integer64_array(comm, data, dest, tag)
   !! Blocking send of an integer64 array to specified destination
      type(comm_t), intent(in) :: comm
      integer(int64), intent(in) :: data(:)
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_send_integer64_array needs a peer rank"
   end subroutine comm_send_integer64_array

   subroutine comm_send_real_dp(comm, data, dest, tag)
      !! Blocking send of a single double precision real to specified destination
      type(comm_t), intent(in) :: comm
      real(dp), intent(in) :: data
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_send_real_dp needs a peer rank"
   end subroutine comm_send_real_dp

   subroutine comm_send_real_dp_array(comm, data, dest, tag)
      !! Blocking send of a double precision real array to specified destination
      type(comm_t), intent(in) :: comm
      real(dp), intent(in) :: data(:)
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_send_real_dp_array needs a peer rank"
   end subroutine comm_send_real_dp_array

   subroutine comm_send_real_dp_array_2d(comm, data, dest, tag)
   !! Blocking send of a 2D double precision real array to specified destination
      type(comm_t), intent(in) :: comm
      real(dp), intent(in) :: data(:, :)
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr, dim1, dim2

      ! Send dimensions first
      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_send_real_dp_array_2d needs a peer rank"
   end subroutine comm_send_real_dp_array_2d

   subroutine comm_send_real_dp_array_2d_n(comm, data, count, dest, tag)
      !! Blocking send of a contiguous 2D double-precision array using
      !! an explicit count.  Unlike `comm_send_real_dp_array_2d`, this
      !! variant does NOT prefix the message with the dimensions —
      !! caller and receiver agree on the shape via the protocol.  The
      !! receiver uses `comm_recv_real_dp_array_2d_n` (or any of the
      !! `_array_n` overloads) with the same `count`.
      type(comm_t), intent(in) :: comm
      real(dp), intent(in) :: data(:, :)
      integer(int32), intent(in) :: count
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_send_real_dp_array_2d_n needs a peer rank"
   end subroutine comm_send_real_dp_array_2d_n

   subroutine comm_send_integer_array_2d_n(comm, data, count, dest, tag)
      !! Blocking send of a contiguous 2D int32 array using an explicit
      !! count, no dim-prefix.  Pairs with `comm_recv_integer_array_2d_n`.
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: data(:, :)
      integer(int32), intent(in) :: count
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_send_integer_array_2d_n needs a peer rank"
   end subroutine comm_send_integer_array_2d_n

   subroutine comm_send_real_sp(comm, data, dest, tag)
      !! Blocking send of a single single-precision real to specified destination
      type(comm_t), intent(in) :: comm
      real(sp), intent(in) :: data
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_send_real_sp needs a peer rank"
   end subroutine comm_send_real_sp

   subroutine comm_send_real_sp_array(comm, data, dest, tag)
      !! Blocking send of a single-precision real array to specified destination
      type(comm_t), intent(in) :: comm
      real(sp), intent(in) :: data(:)
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_send_real_sp_array needs a peer rank"
   end subroutine comm_send_real_sp_array

   subroutine comm_send_real_sp_array_2d(comm, data, dest, tag)
   !! Blocking send of a 2D single-precision real array to specified destination
      type(comm_t), intent(in) :: comm
      real(sp), intent(in) :: data(:, :)
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr, dim1, dim2

      ! Send dimensions first
      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_send_real_sp_array_2d needs a peer rank"
   end subroutine comm_send_real_sp_array_2d

   subroutine comm_send_logical(comm, data, dest, tag)
   !! Blocking send of a logical value to specified destination
      type(comm_t), intent(in) :: comm
      logical, intent(in) :: data
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_send_logical needs a peer rank"
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

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_recv_integer needs a peer rank"
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
      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_recv_integer_array needs a peer rank"
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

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_recv_integer64 needs a peer rank"
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
      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_recv_integer64_array needs a peer rank"
   end subroutine comm_recv_integer64_array

   subroutine comm_recv_real_dp(comm, data, source, tag, status)
      !! Blocking receive of a single double precision real from specified source.
      type(comm_t), intent(in) :: comm
      real(dp), intent(out) :: data
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr
      type(MPI_Status), intent(out), optional :: status

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_recv_real_dp needs a peer rank"
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
      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_recv_real_dp_array needs a peer rank"
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
      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_recv_real_dp_array_2d needs a peer rank"
   end subroutine comm_recv_real_dp_array_2d

   subroutine comm_recv_real_dp_array_2d_n(comm, data, count, source, tag, status)
      !! Blocking receive of a contiguous 2D double-precision array
      !! with explicit `count`.  No dim-prefix protocol — the caller
      !! has already shaped `data`.  Mirrors `comm_send_real_dp_array_2d_n`.
      type(comm_t), intent(in)    :: comm
      real(dp), intent(out)   :: data(:, :)
      integer(int32), intent(in)    :: count
      integer(int32), intent(in)    :: source
      integer(int32), intent(in)    :: tag
      type(MPI_Status), intent(out), optional :: status
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_recv_real_dp_array_2d_n needs a peer rank"
   end subroutine comm_recv_real_dp_array_2d_n

   subroutine comm_recv_integer_array_2d_n(comm, data, count, source, tag, status)
      !! Blocking recv of a contiguous 2D int32 array with explicit count.
      type(comm_t), intent(in)    :: comm
      integer(int32), intent(out)   :: data(:, :)
      integer(int32), intent(in)    :: count
      integer(int32), intent(in)    :: source
      integer(int32), intent(in)    :: tag
      type(MPI_Status), intent(out), optional :: status
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_recv_integer_array_2d_n needs a peer rank"
   end subroutine comm_recv_integer_array_2d_n

   subroutine comm_recv_real_sp(comm, data, source, tag, status)
      !! Blocking receive of a single single-precision real from specified source.
      type(comm_t), intent(in) :: comm
      real(sp), intent(out) :: data
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr
      type(MPI_Status), intent(out), optional :: status

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_recv_real_sp needs a peer rank"
   end subroutine comm_recv_real_sp

   subroutine comm_recv_real_sp_array(comm, data, source, tag, status)
      !! Blocking receive of a single-precision real array from specified source.
      type(comm_t), intent(in) :: comm
      real(sp), allocatable, intent(out) :: data(:)
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(MPI_Status) :: status
      integer(int32) :: count
      integer(int32) :: ierr

      ! First probe to get message size
      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_recv_real_sp_array needs a peer rank"
   end subroutine comm_recv_real_sp_array

   subroutine comm_recv_real_sp_array_2d(comm, data, source, tag, status)
   !! Blocking receive of a 2D allocatable single-precision real array
      type(comm_t), intent(in) :: comm
      real(sp), intent(inout), allocatable :: data(:, :)
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(MPI_Status), intent(out) :: status
      integer(int32) :: ierr, count, dim1, dim2

      ! Receive dimensions first
      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_recv_real_sp_array_2d needs a peer rank"
   end subroutine comm_recv_real_sp_array_2d

   subroutine comm_recv_logical(comm, data, source, tag, status)
   !! Blocking receive of a logical value from specified source
      type(comm_t), intent(in) :: comm
      logical, intent(out) :: data
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(MPI_Status), intent(out), optional :: status
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_recv_logical needs a peer rank"
   end subroutine comm_recv_logical

   subroutine comm_recv_real_dp_array_n(comm, data, count, source, tag, status)
      !! Blocking receive into a pre-allocated double-precision array.
      type(comm_t), intent(in) :: comm
      real(dp), intent(out) :: data(:)
      integer(int32), intent(in) :: count
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(MPI_Status), intent(out), optional :: status
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_recv_real_dp_array_n needs a peer rank"
   end subroutine comm_recv_real_dp_array_n

   subroutine comm_recv_real_sp_array_n(comm, data, count, source, tag, status)
      !! Blocking receive into a pre-allocated single-precision array.
      type(comm_t), intent(in) :: comm
      real(sp), intent(out) :: data(:)
      integer(int32), intent(in) :: count
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(MPI_Status), intent(out), optional :: status
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_recv_real_sp_array_n needs a peer rank"
   end subroutine comm_recv_real_sp_array_n

   subroutine comm_recv_integer_array_n(comm, data, count, source, tag, status)
      !! Blocking receive into a pre-allocated int32 array.
      type(comm_t), intent(in) :: comm
      integer(int32), intent(out) :: data(:)
      integer(int32), intent(in) :: count
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(MPI_Status), intent(out), optional :: status
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_recv_integer_array_n needs a peer rank"
   end subroutine comm_recv_integer_array_n

   subroutine comm_recv_integer64_array_n(comm, data, count, source, tag, status)
      !! Blocking receive into a pre-allocated int64 array.
      type(comm_t), intent(in) :: comm
      integer(int64), intent(out) :: data(:)
      integer(int32), intent(in) :: count
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(MPI_Status), intent(out), optional :: status
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_recv_integer64_array_n needs a peer rank"
   end subroutine comm_recv_integer64_array_n

   subroutine comm_iprobe(comm, source, tag, message_pending, status)
      !! Non-blocking probe for incoming messages
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      logical, intent(out) :: message_pending
      type(MPI_Status), intent(out) :: status
      integer(int32) :: ierr

      ! No peer exists, so nothing is ever pending.
      ! Reporting false keeps a polling loop terminating
      ! rather than spinning forever.
      message_pending = .false.
   end subroutine comm_iprobe

   subroutine comm_probe(comm, source, tag, status)
      !! Blocking probe for incoming messages.  Returns once a message
      !! matching `(source, tag)` is queued at the receiver — caller
      !! reads `status%MPI_SOURCE` and `status%MPI_TAG` to decide who's
      !! talking and what kind of message it is, then issues the
      !! matching `recv`.  Use `MPI_ANY_SOURCE` / `MPI_ANY_TAG` to
      !! dispatch on whichever rank speaks first.
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(MPI_Status), intent(out) :: status
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_probe needs a peer rank"
   end subroutine comm_probe

   subroutine comm_finalize(this)
      !! Frees the MPI communicator resources
      class(comm_t), intent(inout) :: this
      integer(int32) :: ierr

      ! Nothing to synchronise or release on one rank.
   end subroutine comm_finalize

   subroutine abort_comm(comm, errorcode)
      !! Aborts all processes in the communicator with the given error code
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: errorcode
      integer(int32) :: ierr

      ! No ranks to bring down with us.
      error stop "pic_mpi_serial: abort_comm"
   end subroutine abort_comm

   subroutine comm_allgather_integer(comm, sendbuf, recvbuf)
      !! Gathers integer values from all processes in the communicator
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: sendbuf
      integer(int32), intent(out) :: recvbuf(:)
      integer(int32) :: ierr

      ! One rank: the gathered/reduced result is this
      ! rank's own contribution.
      recvbuf = sendbuf
   end subroutine comm_allgather_integer

   subroutine comm_bcast_integer(comm, buffer, count, root)
      !! Broadcasts integer data from root process to all processes in communicator
      type(comm_t), intent(in) :: comm
      integer(int32), intent(inout) :: buffer
      integer(int32), intent(in) :: count
      integer(int32), intent(in) :: root
      integer(int32) :: ierr

      ! One rank: the data is already where it needs to be.
   end subroutine comm_bcast_integer

   subroutine comm_bcast_integer64(comm, buffer, count, root)
      !! Broadcasts integer64 data from root process to all processes in communicator
      type(comm_t), intent(in) :: comm
      integer(int64), intent(inout) :: buffer
      integer(int32), intent(in) :: count
      integer(int32), intent(in) :: root
      integer(int32) :: ierr

      ! One rank: the data is already where it needs to be.
   end subroutine comm_bcast_integer64

   subroutine comm_bcast_real_dp(comm, buffer, count, root)
      !! Broadcasts double precision data from root process to all processes in communicator
      type(comm_t), intent(in) :: comm
      real(dp), intent(inout) :: buffer
      integer(int32), intent(in) :: count
      integer(int32), intent(in) :: root
      integer(int32) :: ierr

      ! One rank: the data is already where it needs to be.
   end subroutine comm_bcast_real_dp

   subroutine comm_bcast_real_dp_array(comm, buffer, count, root)
      !! Broadcasts double precision array from root process to all processes in communicator
      type(comm_t), intent(in) :: comm
      real(dp), intent(inout) :: buffer(:)
      integer(int32), intent(in) :: count
      integer(int32), intent(in) :: root
      integer(int32) :: ierr

      ! One rank: the data is already where it needs to be.
   end subroutine comm_bcast_real_dp_array

   subroutine comm_bcast_real_sp(comm, buffer, count, root)
      !! Broadcasts single-precision data from root process to all processes in communicator
      type(comm_t), intent(in) :: comm
      real(sp), intent(inout) :: buffer
      integer(int32), intent(in) :: count
      integer(int32), intent(in) :: root
      integer(int32) :: ierr

      ! One rank: the data is already where it needs to be.
   end subroutine comm_bcast_real_sp

   subroutine comm_bcast_real_sp_array(comm, buffer, count, root)
      !! Broadcasts single-precision array from root process to all processes in communicator
      type(comm_t), intent(in) :: comm
      real(sp), intent(inout) :: buffer(:)
      integer(int32), intent(in) :: count
      integer(int32), intent(in) :: root
      integer(int32) :: ierr

      ! One rank: the data is already where it needs to be.
   end subroutine comm_bcast_real_sp_array

   subroutine get_processor_name(name, namelen)
      !! Retrieves the name of the processor
      character(len=*), intent(inout) :: name
      integer(int32), intent(out) :: namelen
      integer(int32) :: ierr

      name = "serial"
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
      ! No MPI to start. A serial build is always as
      ! thread-safe as it is ever going to be.
      if (present(provided_thread_level)) provided_thread_level = MPI_THREAD_MULTIPLE
   end subroutine pic_mpi_init

   function pic_mpi_query_thread_level() result(thread_level)
      !! Query the current MPI thread support level
      integer(int32) :: thread_level
      integer(int32) :: ierr

      thread_level = MPI_THREAD_MULTIPLE
   end function pic_mpi_query_thread_level

   subroutine pic_mpi_finalize(ierr)
      !! Finalize MPI environment
      integer(int32), optional, intent(out) :: ierr
      integer(int32) :: ierr_local
      ! Nothing to synchronise or release on one rank.
   end subroutine pic_mpi_finalize

   pure function request_is_null(this) result(is_null)
      !! Checks if the request is null/invalid
      class(request_t), intent(in) :: this
      logical :: is_null
      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: request_is_null needs a peer rank"
   end function request_is_null

   function request_get(this) result(mpi_request_out)
      !! Retrieves the underlying MPI_Request handle
      class(request_t), intent(in) :: this
      type(MPI_Request) :: mpi_request_out

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: request_get needs a peer rank"
   end function request_get

   subroutine request_free(this)
      !! Frees the MPI request resources
      class(request_t), intent(inout) :: this
      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: request_free needs a peer rank"
   end subroutine request_free

   subroutine comm_isend_integer(comm, data, dest, tag, request)
   !! Initiates a non-blocking send operation. The request must be
   !! waited on using wait() or test() before the buffer can be reused.
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: data
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_isend_integer needs a peer rank"
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

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_isend_integer_array needs a peer rank"
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

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_isend_integer64 needs a peer rank"
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

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_isend_integer64_array needs a peer rank"
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

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_isend_real_dp needs a peer rank"
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

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_isend_real_dp_array needs a peer rank"
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
      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_isend_real_dp_array_2d needs a peer rank"
   end subroutine comm_isend_real_dp_array_2d

   subroutine comm_isend_real_dp_array_2d_n(comm, data, count, dest, tag, request)
      !! Non-blocking send of a contiguous 2D double-precision array
      !! with explicit count — no dim-prefix protocol.  Pairs with
      !! `comm_recv_real_dp_array_2d_n` on the receive side.  The
      !! caller must keep `data` valid until `wait`/`waitall` on the
      !! returned request completes.
      type(comm_t), intent(in)  :: comm
      real(dp), intent(in)  :: data(:, :)
      integer(int32), intent(in)  :: count
      integer(int32), intent(in)  :: dest
      integer(int32), intent(in)  :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_isend_real_dp_array_2d_n needs a peer rank"
   end subroutine comm_isend_real_dp_array_2d_n

   subroutine comm_isend_logical(comm, data, dest, tag, request)
   !! Non-blocking send of a logical value
      type(comm_t), intent(in) :: comm
      logical, intent(in) :: data
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_isend_logical needs a peer rank"
   end subroutine comm_isend_logical

   subroutine comm_isend_real_dp_array_n(comm, data, count, dest, tag, request)
      !! Non-blocking send with explicit count (for device pointers in host_data blocks)
      type(comm_t), intent(in) :: comm
      real(dp) :: data(:)
      integer(int32), intent(in) :: count
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_isend_real_dp_array_n needs a peer rank"
   end subroutine comm_isend_real_dp_array_n

   subroutine comm_isend_real_sp(comm, data, dest, tag, request)
   !! Non-blocking send of a single single-precision real
      type(comm_t), intent(in) :: comm
      real(sp), intent(in) :: data
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_isend_real_sp needs a peer rank"
   end subroutine comm_isend_real_sp

   subroutine comm_isend_real_sp_array(comm, data, dest, tag, request)
      !! Non-blocking send of a single-precision real array
      type(comm_t), intent(in) :: comm
      real(sp), intent(in) :: data(:)
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_isend_real_sp_array needs a peer rank"
   end subroutine comm_isend_real_sp_array

   subroutine comm_isend_real_sp_array_2d(comm, data, dest, tag, request)
   !! Non-blocking send of a 2D single-precision real array
      type(comm_t), intent(in) :: comm
      real(sp), intent(in) :: data(:, :)
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr, dim1, dim2

      ! Send dimensions first (blocking - simple approach)
      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_isend_real_sp_array_2d needs a peer rank"
   end subroutine comm_isend_real_sp_array_2d

   subroutine comm_isend_real_sp_array_n(comm, data, count, dest, tag, request)
      !! Non-blocking send with explicit count for single-precision (for device pointers in host_data blocks)
      type(comm_t), intent(in) :: comm
      real(sp) :: data(:)
      integer(int32), intent(in) :: count
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_isend_real_sp_array_n needs a peer rank"
   end subroutine comm_isend_real_sp_array_n

   subroutine comm_irecv_integer(comm, data, source, tag, request)
      !! Initiates a non-blocking receive operation. The request must be
      !! waited on using wait() or test() before the buffer can be used.
      type(comm_t), intent(in) :: comm
      integer(int32), intent(out) :: data
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_irecv_integer needs a peer rank"
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

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_irecv_integer_array needs a peer rank"
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

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_irecv_integer64 needs a peer rank"
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

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_irecv_integer64_array needs a peer rank"
   end subroutine comm_irecv_integer64_array

   subroutine comm_irecv_real_dp(comm, data, source, tag, request)
      type(comm_t), intent(in) :: comm
      real(dp), intent(out) :: data
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_irecv_real_dp needs a peer rank"
   end subroutine comm_irecv_real_dp

   subroutine comm_irecv_real_dp_array(comm, data, source, tag, request)
      type(comm_t), intent(in) :: comm
      real(dp), intent(out) :: data(:)
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_irecv_real_dp_array needs a peer rank"
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
      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_irecv_real_dp_array_2d needs a peer rank"
   end subroutine comm_irecv_real_dp_array_2d

   subroutine comm_irecv_logical(comm, data, source, tag, request)
   !! Non-blocking receive of a logical value
      type(comm_t), intent(in) :: comm
      logical, intent(out) :: data
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_irecv_logical needs a peer rank"
   end subroutine comm_irecv_logical

   subroutine comm_irecv_real_dp_array_n(comm, data, count, source, tag, request)
      !! Non-blocking receive with explicit count (for device pointers in host_data blocks)
      type(comm_t), intent(in) :: comm
      real(dp) :: data(:)
      integer(int32), intent(in) :: count
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_irecv_real_dp_array_n needs a peer rank"
   end subroutine comm_irecv_real_dp_array_n

   subroutine comm_irecv_real_sp(comm, data, source, tag, request)
      !! Non-blocking receive of a single single-precision real
      type(comm_t), intent(in) :: comm
      real(sp), intent(out) :: data
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_irecv_real_sp needs a peer rank"
   end subroutine comm_irecv_real_sp

   subroutine comm_irecv_real_sp_array(comm, data, source, tag, request)
      !! Non-blocking receive of a single-precision real array
      type(comm_t), intent(in) :: comm
      real(sp), intent(out) :: data(:)
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_irecv_real_sp_array needs a peer rank"
   end subroutine comm_irecv_real_sp_array

   subroutine comm_irecv_real_sp_array_2d(comm, data, source, tag, request)
   !! Non-blocking receive of a 2D allocatable single-precision real array
      type(comm_t), intent(in) :: comm
      real(sp), intent(inout), allocatable :: data(:, :)
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr, dim1, dim2

      ! Receive dimensions first (blocking - needed to allocate)
      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_irecv_real_sp_array_2d needs a peer rank"
   end subroutine comm_irecv_real_sp_array_2d

   subroutine comm_irecv_real_sp_array_n(comm, data, count, source, tag, request)
      !! Non-blocking receive with explicit count for single-precision (for device pointers in host_data blocks)
      type(comm_t), intent(in) :: comm
      real(sp) :: data(:)
      integer(int32), intent(in) :: count
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: comm_irecv_real_sp_array_n needs a peer rank"
   end subroutine comm_irecv_real_sp_array_n

   subroutine request_wait(request, status)
   !! Blocks until the operation associated with the request completes.
   !! The request is freed after completion.
      type(request_t), intent(inout) :: request
      type(MPI_Status), intent(out), optional :: status
      integer(int32) :: ierr

      ! Nothing was ever started, so nothing is outstanding.
   end subroutine request_wait

   subroutine request_waitall(requests, statuses)
   !! Blocks until all operations in the request array complete.
   !! All requests are freed after completion.
      type(request_t), intent(inout) :: requests(:)
      type(MPI_Status), intent(out), optional :: statuses(:)
      type(MPI_Request), allocatable :: mpi_requests(:)
      type(MPI_Status), allocatable :: temp_statuses(:)
      integer(int32) :: i, count, ierr

      ! Nothing was ever started, so nothing is outstanding.
   end subroutine request_waitall

   subroutine request_test(request, flag, status)
      type(request_t), intent(inout) :: request
      logical, intent(out) :: flag
      type(MPI_Status), intent(out), optional :: status
      integer(int32) :: ierr

      ! Nothing was started, so it is complete.
      flag = .true.
   end subroutine request_test

   function create_win_dp_array(comm, base, win_size) result(win)
      !! Create MPI window for RMA operations
      !!
      !! Creates a window exposing local memory to remote RMA operations.
      !! Used for DDI distributed arrays.
      type(comm_t), intent(in) :: comm
      real(dp), target :: base(:)
      integer(MPI_ADDRESS_KIND), intent(in) :: win_size
      type(win_t) :: win
      integer(int32) :: ierr
      integer(int32) :: disp_unit

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: create_win_dp_array needs a peer rank"
   end function create_win_dp_array

   function create_win_dynamic(comm) result(win)
      !! Create dynamic MPI window
      !!
      !! For windows where memory will be attached later.
      !! Useful for load balancing counters.
      type(comm_t), intent(in) :: comm
      type(win_t) :: win
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: create_win_dynamic needs a peer rank"
   end function create_win_dynamic

   subroutine create_win_allocate_dp_1d(comm, length, baseptr, win)
      use iso_c_binding, only: c_ptr, c_f_pointer
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: length
      real(dp), pointer, intent(out) :: baseptr(:)
      type(win_t), intent(out) :: win
      type(c_ptr) :: cptr
      integer(MPI_ADDRESS_KIND) :: win_size
      integer(int32) :: disp_unit, ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: create_win_allocate_dp_1d needs a peer rank"
   end subroutine create_win_allocate_dp_1d

   subroutine create_win_allocate_dp_2d(comm, dim1, dim2, baseptr, win)
      use iso_c_binding, only: c_ptr, c_f_pointer
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: dim1, dim2
      real(dp), pointer, intent(out) :: baseptr(:, :)
      type(win_t), intent(out) :: win
      type(c_ptr) :: cptr
      integer(MPI_ADDRESS_KIND) :: win_size
      integer(int32) :: disp_unit, ierr, total_size

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: create_win_allocate_dp_2d needs a peer rank"
   end subroutine create_win_allocate_dp_2d

   subroutine create_win_allocate_i64_1d(comm, length, baseptr, win)
      use iso_c_binding, only: c_ptr, c_f_pointer
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: length
      integer(int64), pointer, intent(out) :: baseptr(:)
      type(win_t), intent(out) :: win
      type(c_ptr) :: cptr
      integer(MPI_ADDRESS_KIND) :: win_size
      integer(int32) :: disp_unit, ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: create_win_allocate_i64_1d needs a peer rank"
   end subroutine create_win_allocate_i64_1d

   subroutine create_win_allocate_sp_1d(comm, length, baseptr, win)
      use iso_c_binding, only: c_ptr, c_f_pointer
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: length
      real(sp), pointer, intent(out) :: baseptr(:)
      type(win_t), intent(out) :: win
      type(c_ptr) :: cptr
      integer(MPI_ADDRESS_KIND) :: win_size
      integer(int32) :: disp_unit, ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: create_win_allocate_sp_1d needs a peer rank"
   end subroutine create_win_allocate_sp_1d

   subroutine create_win_allocate_i32_1d(comm, length, baseptr, win)
      use iso_c_binding, only: c_ptr, c_f_pointer
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: length
      integer(int32), pointer, intent(out) :: baseptr(:)
      type(win_t), intent(out) :: win
      type(c_ptr) :: cptr
      integer(MPI_ADDRESS_KIND) :: win_size
      integer(int32) :: disp_unit, ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: create_win_allocate_i32_1d needs a peer rank"
   end subroutine create_win_allocate_i32_1d

   pure function win_is_null(this) result(is_null)
      class(win_t), intent(in) :: this
      logical :: is_null
      is_null = .true.
   end function win_is_null

   function win_get_handle(this) result(mpi_win_out)
      class(win_t), intent(in) :: this
      type(MPI_Win) :: mpi_win_out

      mpi_win_out = MPI_WIN_NULL
   end function win_get_handle

   subroutine win_fence(this, assert)
      class(win_t), intent(in) :: this
      integer(int32), intent(in), optional :: assert
      integer(int32) :: ierr, assert_val

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: win_fence needs a peer rank"
   end subroutine win_fence

   subroutine win_lock(this, rank, lock_type)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: rank
      integer(int32), intent(in), optional :: lock_type
      integer(int32) :: ierr, ltype

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: win_lock needs a peer rank"
   end subroutine win_lock

   subroutine win_unlock(this, rank)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: rank
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: win_unlock needs a peer rank"
   end subroutine win_unlock

   subroutine win_lock_all(this, assert)
      class(win_t), intent(in) :: this
      integer(int32), intent(in), optional :: assert
      integer(int32) :: ierr, assert_val

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: win_lock_all needs a peer rank"
   end subroutine win_lock_all

   subroutine win_unlock_all(this)
      class(win_t), intent(in) :: this
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: win_unlock_all needs a peer rank"
   end subroutine win_unlock_all

   subroutine win_flush(this, rank)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: rank
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: win_flush needs a peer rank"
   end subroutine win_flush

   subroutine win_flush_all(this)
      class(win_t), intent(in) :: this
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: win_flush_all needs a peer rank"
   end subroutine win_flush_all

   subroutine win_get_dp(this, target_rank, target_disp, count, buffer)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      real(dp), intent(out) :: buffer(*)
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: win_get_dp needs a peer rank"
   end subroutine win_get_dp

   subroutine win_put_dp(this, target_rank, target_disp, count, buffer)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      real(dp), intent(in) :: buffer(*)
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: win_put_dp needs a peer rank"
   end subroutine win_put_dp

   subroutine win_rget_dp(this, target_rank, target_disp, count, buffer, request)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      real(dp), intent(out) :: buffer(*)
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: win_rget_dp needs a peer rank"
   end subroutine win_rget_dp

   subroutine win_rput_dp(this, target_rank, target_disp, count, buffer, request)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      real(dp), intent(in) :: buffer(*)
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: win_rput_dp needs a peer rank"
   end subroutine win_rput_dp

   subroutine win_accumulate_dp(this, target_rank, target_disp, count, buffer, op)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      real(dp), intent(in) :: buffer(*)
      type(MPI_Op), intent(in), optional :: op
      type(MPI_Op) :: op_used
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: win_accumulate_dp needs a peer rank"
   end subroutine win_accumulate_dp

   subroutine win_get_sp(this, target_rank, target_disp, count, buffer)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      real(sp), intent(out) :: buffer(*)
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: win_get_sp needs a peer rank"
   end subroutine win_get_sp

   subroutine win_put_sp(this, target_rank, target_disp, count, buffer)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      real(sp), intent(in) :: buffer(*)
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: win_put_sp needs a peer rank"
   end subroutine win_put_sp

   subroutine win_rget_sp(this, target_rank, target_disp, count, buffer, request)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      real(sp), intent(out) :: buffer(*)
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: win_rget_sp needs a peer rank"
   end subroutine win_rget_sp

   subroutine win_rput_sp(this, target_rank, target_disp, count, buffer, request)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      real(sp), intent(in) :: buffer(*)
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: win_rput_sp needs a peer rank"
   end subroutine win_rput_sp

   subroutine win_accumulate_sp(this, target_rank, target_disp, count, buffer, op)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      real(sp), intent(in) :: buffer(*)
      type(MPI_Op), intent(in), optional :: op
      type(MPI_Op) :: op_used
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: win_accumulate_sp needs a peer rank"
   end subroutine win_accumulate_sp

   subroutine win_get_i32(this, target_rank, target_disp, count, buffer)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      integer(int32), intent(out) :: buffer(*)
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: win_get_i32 needs a peer rank"
   end subroutine win_get_i32

   subroutine win_put_i32(this, target_rank, target_disp, count, buffer)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      integer(int32), intent(in) :: buffer(*)
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: win_put_i32 needs a peer rank"
   end subroutine win_put_i32

   subroutine win_rget_i32(this, target_rank, target_disp, count, buffer, request)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      integer(int32), intent(out) :: buffer(*)
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: win_rget_i32 needs a peer rank"
   end subroutine win_rget_i32

   subroutine win_rput_i32(this, target_rank, target_disp, count, buffer, request)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      integer(int32), intent(in) :: buffer(*)
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: win_rput_i32 needs a peer rank"
   end subroutine win_rput_i32

   subroutine win_accumulate_i32(this, target_rank, target_disp, count, buffer, op)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      integer(int32), intent(in) :: buffer(*)
      type(MPI_Op), intent(in), optional :: op
      type(MPI_Op) :: op_used
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: win_accumulate_i32 needs a peer rank"
   end subroutine win_accumulate_i32

   subroutine win_get_i64(this, target_rank, target_disp, count, buffer)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      integer(int64), intent(out) :: buffer(*)
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: win_get_i64 needs a peer rank"
   end subroutine win_get_i64

   subroutine win_put_i64(this, target_rank, target_disp, count, buffer)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      integer(int64), intent(in) :: buffer(*)
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: win_put_i64 needs a peer rank"
   end subroutine win_put_i64

   subroutine win_rget_i64(this, target_rank, target_disp, count, buffer, request)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      integer(int64), intent(out) :: buffer(*)
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: win_rget_i64 needs a peer rank"
   end subroutine win_rget_i64

   subroutine win_rput_i64(this, target_rank, target_disp, count, buffer, request)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      integer(int64), intent(in) :: buffer(*)
      type(request_t), intent(out) :: request
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: win_rput_i64 needs a peer rank"
   end subroutine win_rput_i64

   subroutine win_accumulate_i64(this, target_rank, target_disp, count, buffer, op)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int32), intent(in) :: count
      integer(int64), intent(in) :: buffer(*)
      type(MPI_Op), intent(in), optional :: op
      type(MPI_Op) :: op_used
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: win_accumulate_i64 needs a peer rank"
   end subroutine win_accumulate_i64

   subroutine win_fetch_and_add_i64(this, target_rank, target_disp, value, result)
      class(win_t), intent(in) :: this
      integer(int32), intent(in) :: target_rank
      integer(MPI_ADDRESS_KIND), intent(in) :: target_disp
      integer(int64), intent(in) :: value
      integer(int64), intent(out) :: result
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: win_fetch_and_add_i64 needs a peer rank"
   end subroutine win_fetch_and_add_i64

   subroutine win_finalize(this)
      class(win_t), intent(inout) :: this
      integer(int32) :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: win_finalize needs a peer rank"
   end subroutine win_finalize

   subroutine allreduce_dp(comm, buffer, op)
      type(comm_t), intent(in) :: comm
      real(dp), intent(inout) :: buffer
      type(MPI_Op), intent(in), optional :: op
      type(MPI_Op) :: op_used
      integer :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: allreduce_dp needs a peer rank"
   end subroutine allreduce_dp

   subroutine allreduce_dp_array(comm, buffer, count, op)
      type(comm_t), intent(in) :: comm
      real(dp), intent(inout) :: buffer(:)
      integer, intent(in), optional :: count
      type(MPI_Op), intent(in), optional :: op
      type(MPI_Op) :: op_used
      integer :: ierr, n

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: allreduce_dp_array needs a peer rank"
   end subroutine allreduce_dp_array

   subroutine allreduce_i32(comm, buffer, op)
      type(comm_t), intent(in) :: comm
      integer(int32), intent(inout) :: buffer
      type(MPI_Op), intent(in), optional :: op
      type(MPI_Op) :: op_used
      integer :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: allreduce_i32 needs a peer rank"
   end subroutine allreduce_i32

   subroutine allreduce_i32_array(comm, buffer, count, op)
      type(comm_t), intent(in) :: comm
      integer(int32), intent(inout) :: buffer(:)
      integer, intent(in), optional :: count
      type(MPI_Op), intent(in), optional :: op
      type(MPI_Op) :: op_used
      integer :: ierr, n

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: allreduce_i32_array needs a peer rank"
   end subroutine allreduce_i32_array

   subroutine allreduce_dp_to(comm, sendbuf, recvbuf, op)
      type(comm_t), intent(in) :: comm
      real(dp), intent(in) :: sendbuf
      real(dp), intent(out) :: recvbuf
      type(MPI_Op), intent(in), optional :: op
      type(MPI_Op) :: op_used
      integer :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: allreduce_dp_to needs a peer rank"
   end subroutine allreduce_dp_to

   subroutine allreduce_dp_array_to(comm, sendbuf, recvbuf, count, op)
      type(comm_t), intent(in) :: comm
      real(dp), intent(in) :: sendbuf(:)
      real(dp), intent(out) :: recvbuf(:)
      integer, intent(in), optional :: count
      type(MPI_Op), intent(in), optional :: op
      type(MPI_Op) :: op_used
      integer :: ierr, n

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: allreduce_dp_array_to needs a peer rank"
   end subroutine allreduce_dp_array_to

   subroutine allreduce_sp(comm, buffer, op)
      type(comm_t), intent(in) :: comm
      real(sp), intent(inout) :: buffer
      type(MPI_Op), intent(in), optional :: op
      type(MPI_Op) :: op_used
      integer :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: allreduce_sp needs a peer rank"
   end subroutine allreduce_sp

   subroutine allreduce_sp_array(comm, buffer, count, op)
      type(comm_t), intent(in) :: comm
      real(sp), intent(inout) :: buffer(:)
      integer, intent(in), optional :: count
      type(MPI_Op), intent(in), optional :: op
      type(MPI_Op) :: op_used
      integer :: ierr, n

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: allreduce_sp_array needs a peer rank"
   end subroutine allreduce_sp_array

   subroutine allreduce_sp_to(comm, sendbuf, recvbuf, op)
      type(comm_t), intent(in) :: comm
      real(sp), intent(in) :: sendbuf
      real(sp), intent(out) :: recvbuf
      type(MPI_Op), intent(in), optional :: op
      type(MPI_Op) :: op_used
      integer :: ierr

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: allreduce_sp_to needs a peer rank"
   end subroutine allreduce_sp_to

   subroutine allreduce_sp_array_to(comm, sendbuf, recvbuf, count, op)
      type(comm_t), intent(in) :: comm
      real(sp), intent(in) :: sendbuf(:)
      real(sp), intent(out) :: recvbuf(:)
      integer, intent(in), optional :: count
      type(MPI_Op), intent(in), optional :: op
      type(MPI_Op) :: op_used
      integer :: ierr, n

      ! A single rank has no peer. Reaching this means the caller
      ! did not take its `size() == 1` path -- which is a bug worth
      ! stopping for, not one to hide by returning quietly.
      error stop "pic_mpi_serial: allreduce_sp_array_to needs a peer rank"
   end subroutine allreduce_sp_array_to

end module pic_mpi_serial
