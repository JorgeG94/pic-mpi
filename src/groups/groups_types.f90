!> Type definitions for DDI groups and scope management
!!
!! Provides types for managing multiple communicator groups and scope switching
!! for DDI compatibility. Supports creating subgroups from a world communicator
!! and switching the active working communicator for array operations.
module groups_types
   use pic_types, only: int32
   use pic_mpi_lib, only: comm_t
   implicit none
   private

   public :: ddi_group_t
   public :: DDI_COMM_WORLD, DDI_COMM_GROUP, DDI_COMM_MASTER

   ! Standard DDI communicator IDs
   integer(int32), parameter :: DDI_COMM_WORLD = 1   !! World communicator
   integer(int32), parameter :: DDI_COMM_GROUP = 2   !! Current group communicator
   integer(int32), parameter :: DDI_COMM_MASTER = 3  !! Masters-only communicator

   !> DDI group descriptor
   !!
   !! Represents a group configuration with world, group, and master communicators.
   !! Each rank belongs to exactly one group. Group masters (rank 0 in each group)
   !! can communicate via the master communicator.
   type :: ddi_group_t
      type(comm_t) :: world_comm     !! Parent communicator (all ranks)
      type(comm_t) :: group_comm     !! Group communicator (same-group ranks)
      type(comm_t) :: master_comm    !! Masters-only communicator
      integer(int32) :: group_id     !! Which group this rank belongs to (0-indexed)
      integer(int32) :: ngroups      !! Total number of groups
      logical :: is_master           !! True if rank 0 in group
   end type ddi_group_t

end module groups_types
