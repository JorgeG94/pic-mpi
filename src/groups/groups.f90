!> Groups and scope management for DDI compatibility
!!
!! Provides DDI-compatible group creation and communicator scope switching.
!! Maintains a registry of communicators that can be switched for operations.
!!
!! Usage:
!!   call groups_init(world_comm)
!!   call ddi_group_create(ngroups, DDI_COMM_WORLD, DDI_COMM_GROUP, DDI_COMM_MASTER)
!!   call ddi_scope(DDI_COMM_GROUP)
!!   ! Operations now use group communicator
!!   call ddi_scope(DDI_COMM_WORLD)
!!   call groups_finalize()
module groups
   use pic_types, only: int32
   use pic_mpi_lib, only: comm_t
   use mpi_f08, only: MPI_UNDEFINED
   use groups_types, only: ddi_group_t, DDI_COMM_WORLD, DDI_COMM_GROUP, DDI_COMM_MASTER
   implicit none
   private

   ! Re-export types and constants
   public :: ddi_group_t
   public :: DDI_COMM_WORLD, DDI_COMM_GROUP, DDI_COMM_MASTER

   ! Initialization and finalization
   public :: groups_init, groups_finalize

   ! Group operations
   public :: ddi_group_create

   ! Scope operations
   public :: ddi_scope, ddi_ascope, get_working_comm

   ! Communicator registry access
   public :: groups_register_comm, groups_get_comm, groups_get_working_comm_id

   ! Maximum number of registered communicators
   integer(int32), parameter :: MAX_COMMS = 32

   ! Communicator registry
   type(comm_t), save :: comm_registry(MAX_COMMS)
   logical, save :: comm_valid(MAX_COMMS) = .false.
   integer(int32), save :: working_comm_id = DDI_COMM_WORLD

   ! Group state
   type(ddi_group_t), save :: current_group
   logical, save :: groups_initialized = .false.
   logical, save :: groups_created = .false.

contains

   !> Initialize the groups module
   !!
   !! Sets up the communicator registry with the world communicator.
   !! @param comm World communicator
   subroutine groups_init(comm)
      type(comm_t), intent(in) :: comm
      integer(int32) :: i

      ! Clear registry
      do i = 1, MAX_COMMS
         comm_valid(i) = .false.
      end do

      ! Register world communicator at slot 1 (DDI_COMM_WORLD)
      comm_registry(DDI_COMM_WORLD) = comm%duplicate()
      comm_valid(DDI_COMM_WORLD) = .true.
      working_comm_id = DDI_COMM_WORLD

      ! Initialize group state
      current_group%world_comm = comm_registry(DDI_COMM_WORLD)
      current_group%group_id = 0
      current_group%ngroups = 1
      current_group%is_master = comm_registry(DDI_COMM_WORLD)%leader()

      groups_initialized = .true.
      groups_created = .false.
   end subroutine groups_init

   !> Finalize the groups module
   !!
   !! Frees all registered communicators.
   subroutine groups_finalize()
      integer(int32) :: i

      if (.not. groups_initialized) return

      ! Free all valid communicators
      do i = 1, MAX_COMMS
         if (comm_valid(i)) then
            call comm_registry(i)%finalize()
            comm_valid(i) = .false.
         end if
      end do

      groups_initialized = .false.
      groups_created = .false.
      working_comm_id = DDI_COMM_WORLD
   end subroutine groups_finalize

   !> Create DDI-style groups
   !!
   !! Splits the world communicator into ngroups groups with equal ranks.
   !! Creates a masters-only communicator for rank 0 of each group.
   !!
   !! @param ngroups Number of groups to create
   !! @param world_id Output: ID for world communicator (always DDI_COMM_WORLD)
   !! @param group_id Output: ID for group communicator (always DDI_COMM_GROUP)
   !! @param master_id Output: ID for masters communicator (always DDI_COMM_MASTER)
   subroutine ddi_group_create(ngroups, world_id, group_id, master_id)
      integer(int32), intent(in) :: ngroups
      integer(int32), intent(out) :: world_id, group_id, master_id

      type(comm_t) :: world_comm, group_comm, master_comm
      integer(int32) :: world_size, world_rank
      integer(int32) :: ranks_per_group, my_group, color

      if (.not. groups_initialized) then
         error stop "groups: groups_init must be called before ddi_group_create"
      end if

      world_comm = comm_registry(DDI_COMM_WORLD)
      world_size = world_comm%size()
      world_rank = world_comm%rank()

      ! Calculate group assignment
      if (ngroups <= 0) then
         error stop "groups: ngroups must be positive"
      end if

      if (ngroups > world_size) then
         error stop "groups: ngroups cannot exceed number of MPI ranks"
      end if

      ranks_per_group = world_size/ngroups
      my_group = world_rank/ranks_per_group
      if (my_group >= ngroups) my_group = ngroups - 1  ! Handle uneven division

      ! Create group communicator (split by group assignment)
      group_comm = world_comm%split_by(my_group)
      comm_registry(DDI_COMM_GROUP) = group_comm
      comm_valid(DDI_COMM_GROUP) = .true.

      ! Create masters communicator (only rank 0 of each group participates)
      if (group_comm%leader()) then
         color = 0  ! Masters join
      else
         color = MPI_UNDEFINED  ! Non-masters excluded
      end if
      master_comm = world_comm%split_by(color)
      comm_registry(DDI_COMM_MASTER) = master_comm
      comm_valid(DDI_COMM_MASTER) = .true.

      ! Update group state
      current_group%world_comm = world_comm
      current_group%group_comm = group_comm
      current_group%master_comm = master_comm
      current_group%group_id = my_group
      current_group%ngroups = ngroups
      current_group%is_master = group_comm%leader()

      groups_created = .true.

      ! Return communicator IDs
      world_id = DDI_COMM_WORLD
      group_id = DDI_COMM_GROUP
      master_id = DDI_COMM_MASTER
   end subroutine ddi_group_create

   !> Switch the active working communicator (synchronous)
   !!
   !! Changes which communicator is used for subsequent operations.
   !! Performs a barrier on both old and new communicators for synchronization.
   !!
   !! @param comm_id Communicator ID to switch to
   subroutine ddi_scope(comm_id)
      integer(int32), intent(in) :: comm_id

      if (.not. groups_initialized) then
         error stop "groups: groups_init must be called before ddi_scope"
      end if

      if (comm_id < 1 .or. comm_id > MAX_COMMS) then
         error stop "groups: invalid communicator ID in ddi_scope"
      end if

      if (.not. comm_valid(comm_id)) then
         error stop "groups: communicator not registered in ddi_scope"
      end if

      ! Skip if no change
      if (comm_id == working_comm_id) return

      ! Synchronize on old communicator before switching
      if (comm_valid(working_comm_id)) then
         call comm_registry(working_comm_id)%barrier()
      end if

      ! Switch to new communicator
      working_comm_id = comm_id

      ! Synchronize on new communicator
      call comm_registry(working_comm_id)%barrier()
   end subroutine ddi_scope

   !> Switch the active working communicator (asynchronous)
   !!
   !! Changes which communicator is used for subsequent operations.
   !! Does NOT perform any synchronization (caller is responsible).
   !!
   !! @param comm_id Communicator ID to switch to
   subroutine ddi_ascope(comm_id)
      integer(int32), intent(in) :: comm_id

      if (.not. groups_initialized) then
         error stop "groups: groups_init must be called before ddi_ascope"
      end if

      if (comm_id < 1 .or. comm_id > MAX_COMMS) then
         error stop "groups: invalid communicator ID in ddi_ascope"
      end if

      if (.not. comm_valid(comm_id)) then
         error stop "groups: communicator not registered in ddi_ascope"
      end if

      ! Just switch, no synchronization
      working_comm_id = comm_id
   end subroutine ddi_ascope

   !> Get the current working communicator
   !!
   !! @return Current working communicator
   function get_working_comm() result(comm)
      type(comm_t) :: comm

      if (.not. groups_initialized) then
         error stop "groups: groups_init must be called before get_working_comm"
      end if

      comm = comm_registry(working_comm_id)
   end function get_working_comm

   !> Get current working communicator ID
   function groups_get_working_comm_id() result(comm_id)
      integer(int32) :: comm_id
      comm_id = working_comm_id
   end function groups_get_working_comm_id

   !> Register a communicator in the registry
   !!
   !! @param comm Communicator to register
   !! @param comm_id Output: assigned ID
   subroutine groups_register_comm(comm, comm_id)
      type(comm_t), intent(in) :: comm
      integer(int32), intent(out) :: comm_id
      integer(int32) :: i

      if (.not. groups_initialized) then
         error stop "groups: groups_init must be called before groups_register_comm"
      end if

      ! Find free slot (starting from 4 since 1-3 are reserved)
      comm_id = -1
      do i = 4, MAX_COMMS
         if (.not. comm_valid(i)) then
            comm_id = i
            exit
         end if
      end do

      if (comm_id < 0) then
         error stop "groups: no free slots in communicator registry"
      end if

      comm_registry(comm_id) = comm%duplicate()
      comm_valid(comm_id) = .true.
   end subroutine groups_register_comm

   !> Get a communicator by ID
   !!
   !! @param comm_id Communicator ID
   !! @return The communicator
   function groups_get_comm(comm_id) result(comm)
      integer(int32), intent(in) :: comm_id
      type(comm_t) :: comm

      if (.not. groups_initialized) then
         error stop "groups: groups_init must be called before groups_get_comm"
      end if

      if (comm_id < 1 .or. comm_id > MAX_COMMS) then
         error stop "groups: invalid communicator ID"
      end if

      if (.not. comm_valid(comm_id)) then
         error stop "groups: communicator not registered"
      end if

      comm = comm_registry(comm_id)
   end function groups_get_comm

end module groups
