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
                      MPI_Bcast, MPI_Init, MPI_Finalize, &
                      MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_MAX_PROCESSOR_NAME, &
                      operator(==), operator(/=), MPI_DOUBLE_PRECISION
   implicit none
   private

   public :: comm_t, comm_world, comm_null
   public :: send, recv, isend, irecv
   public :: request_t, wait, waitall, test
   public :: iprobe, abort_comm, allgather, get_processor_name, bcast
   public :: pic_mpi_init, pic_mpi_finalize

   ! Export MPI types and constants needed by applications
   public :: MPI_Status, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_MAX_PROCESSOR_NAME

   ! Request type for non-blocking operations
   type :: request_t
      private
      type(MPI_Request) :: m_request = MPI_REQUEST_NULL
      logical :: is_valid = .false.
   contains
      procedure :: is_null => request_is_null
      procedure :: get => request_get
      procedure :: free => request_free
   end type request_t

   type :: comm_t
      private
      type(MPI_Comm) :: m_comm = MPI_COMM_NULL
      integer(int32) :: m_rank = -1
      integer(int32) :: m_size = -1
      logical :: is_valid = .false.
   contains
      procedure :: rank => comm_rank
      procedure :: size => m_size_func
      procedure :: leader => comm_leader
      procedure :: is_null => comm_is_null
      procedure :: get => comm_get

      procedure :: barrier => comm_barrier

      procedure :: split => comm_split_shared
      procedure :: split_by => comm_split_by_color
      procedure :: discard_leader => comm_discard_leader
      procedure :: discard_to => comm_discard_to
      procedure :: duplicate => comm_duplicate

      procedure :: finalize => comm_finalize
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

   function create_world_comm() result(comm)
      type(comm_t) :: comm
      type(MPI_Comm) :: dup_comm
      integer(int32) :: ierr

      call MPI_Comm_dup(MPI_COMM_WORLD, dup_comm, ierr)
      comm = create_comm_from_mpi(dup_comm)

   end function create_world_comm

   function create_null_comm() result(comm)
      type(comm_t) :: comm

      ! Explicitly initialize to null/invalid state
      comm%m_comm = MPI_COMM_NULL
      comm%m_rank = -1
      comm%m_size = -1
      comm%is_valid = .false.
   end function create_null_comm

   pure function comm_rank(this) result(rank)
      class(comm_t), intent(in) :: this
      integer :: rank
      rank = this%m_rank
   end function comm_rank

   pure function m_size_func(this) result(size)
      class(comm_t), intent(in) :: this
      integer :: size
      size = this%m_size
   end function m_size_func

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

   subroutine comm_barrier(this)
      class(comm_t), intent(in) :: this
      integer(int32) :: ierr
      call MPI_Barrier(this%m_comm, ierr)
   end subroutine comm_barrier

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

   subroutine comm_send_integer(comm, data, dest, tag)
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: data
      integer(int32), intent(in) :: dest
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr

      call MPI_Send(data, 1, MPI_INTEGER, dest, tag, comm%m_comm, ierr)
   end subroutine comm_send_integer

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

   subroutine pic_mpi_init(ierr)
      !! Initialize MPI environment
      integer(int32), optional, intent(out) :: ierr
      integer(int32) :: ierr_local
      call MPI_Init(ierr_local)
      if (present(ierr)) ierr = ierr_local
   end subroutine pic_mpi_init

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
