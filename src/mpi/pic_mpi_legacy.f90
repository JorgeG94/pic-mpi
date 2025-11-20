module pic_legacy_mpi
   use mpi
   use pic_types, only: int32
   implicit none
   private

   public :: comm_t, comm_world, comm_null
   public :: send, recv
   public ::  iprobe

   type :: comm_t
      private
      integer :: m_comm = MPI_COMM_NULL
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
   end interface send

   interface recv
      module procedure :: comm_recv_integer
      module procedure :: comm_recv_integer_array
   end interface recv

   interface iprobe
      module procedure :: comm_iprobe
   end interface iprobe

contains

   function create_comm_from_mpi(mpi_comm_in) result(comm)
      integer, intent(in) :: mpi_comm_in
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
      integer :: dup_comm
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
      integer :: mpi_comm_out

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
      integer :: mpi_comm_new
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
      integer :: mpi_comm_new
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
      integer :: mpi_comm_new
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

   subroutine comm_recv_integer(comm, data, source, tag, status)
      type(comm_t), intent(in) :: comm
      integer(int32), intent(out) :: data
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      integer(int32) :: ierr
      integer, intent(out), optional :: status(MPI_STATUS_SIZE)
      integer :: stat(MPI_STATUS_SIZE)

      if (present(status)) then
         call MPI_Recv(data, 1, MPI_INTEGER, source, tag, comm%m_comm, status, ierr)
      else
         call MPI_Recv(data, 1, MPI_INTEGER, source, tag, comm%m_comm, stat, ierr)
      end if
   end subroutine comm_recv_integer

   subroutine comm_recv_integer_array(comm, data, source, tag, status)
      type(comm_t), intent(in) :: comm
      integer(int32), allocatable, intent(out) :: data(:)
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      integer :: status(MPI_STATUS_SIZE)
      integer(int32) :: count
      integer(int32) :: ierr

      ! First probe to get message size
      call MPI_Probe(source, tag, comm%m_comm, status, ierr)
      call MPI_Get_count(status, MPI_INTEGER, count, ierr)

      ! Allocate and receive
      allocate (data(count))
      call MPI_Recv(data, count, MPI_INTEGER, source, tag, comm%m_comm, status, ierr)
   end subroutine comm_recv_integer_array

   subroutine comm_iprobe(comm, source, tag, message_pending, status)
      type(comm_t), intent(in) :: comm
      integer(int32), intent(in) :: source
      integer(int32), intent(in) :: tag
      logical, intent(out) :: message_pending
      integer, intent(out) :: status(MPI_STATUS_SIZE)
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

end module pic_legacy_mpi