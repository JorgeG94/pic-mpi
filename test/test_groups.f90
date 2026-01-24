!> Test groups and scope management module
!!
!! Tests group creation, scope switching, and communicator management.
!! Run with: mpirun -np 4 ./test_groups
program test_groups
   use pic_mpi_lib
   use pic_types, only: int32, int64, dp
   use groups
   use darrays
   implicit none

   type(comm_t) :: world_comm
   integer :: n_passed, n_failed

   n_passed = 0
   n_failed = 0

   call pic_mpi_init()
   world_comm = comm_world()

   if (world_comm%size() < 4) then
      if (world_comm%leader()) then
         print *, "ERROR: This test requires at least 4 MPI ranks"
      end if
      call pic_mpi_finalize()
      stop 1
   end if

   if (world_comm%leader()) then
      print *, "========================================"
      print *, "Testing Groups and Scope Management"
      print *, "========================================"
      print *, "Number of ranks:", world_comm%size()
   end if

   ! Initialize groups
   call groups_init(world_comm)
   call darrays_init(world_comm)

   call world_comm%barrier()

   call test_init_state()
   call test_group_create()
   call test_scope_switch()
   call test_scope_array_isolation()

   call world_comm%barrier()

   ! Summary
   if (world_comm%leader()) then
      print *, "========================================"
      print *, "Results: ", n_passed, " passed, ", n_failed, " failed"
      print *, "========================================"
   end if

   call darrays_finalize()
   call groups_finalize()
   call world_comm%finalize()
   call pic_mpi_finalize()

   if (n_failed > 0) stop 1

contains

   subroutine test_init_state()
      type(comm_t) :: work_comm
      integer(int32) :: np, me
      logical :: test_ok

      if (world_comm%leader()) print *, "Test 1: Initial state..."

      work_comm = get_working_comm()
      np = work_comm%size()
      me = work_comm%rank()

      ! Initially working comm should be world comm
      test_ok = (np == world_comm%size()) .and. (me == world_comm%rank())

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Initial state is world comm"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Initial state is not world comm"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_init_state

   subroutine test_group_create()
      integer(int32) :: world_id, group_id, master_id
      integer(int32) :: ngroups
      type(comm_t) :: group_comm, master_comm
      logical :: test_ok
      integer(int32) :: expected_group_size

      if (world_comm%leader()) print *, "Test 2: Group creation..."

      ngroups = 2
      call ddi_group_create(ngroups, world_id, group_id, master_id)

      test_ok = .true.

      ! Check returned IDs match constants
      if (world_id /= DDI_COMM_WORLD) test_ok = .false.
      if (group_id /= DDI_COMM_GROUP) test_ok = .false.
      if (master_id /= DDI_COMM_MASTER) test_ok = .false.

      ! Check group communicator size (should be world_size / ngroups)
      group_comm = groups_get_comm(DDI_COMM_GROUP)
      expected_group_size = world_comm%size()/ngroups
      if (group_comm%size() /= expected_group_size) then
         test_ok = .false.
      end if

      ! Check masters communicator (only group leaders should have valid comm)
      master_comm = groups_get_comm(DDI_COMM_MASTER)
      if (group_comm%leader()) then
         ! Masters should see ngroups participants
         if (master_comm%size() /= ngroups) test_ok = .false.
      else
         ! Non-masters should have null/invalid comm
         if (.not. master_comm%is_null()) test_ok = .false.
      end if

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Group creation"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Group creation"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_group_create

   subroutine test_scope_switch()
      type(comm_t) :: work_comm
      integer(int32) :: np_world, np_group
      logical :: test_ok

      if (world_comm%leader()) print *, "Test 3: Scope switching..."

      test_ok = .true.

      ! Get world scope size
      call ddi_scope(DDI_COMM_WORLD)
      work_comm = get_working_comm()
      np_world = work_comm%size()

      ! Switch to group scope
      call ddi_scope(DDI_COMM_GROUP)
      work_comm = get_working_comm()
      np_group = work_comm%size()

      ! Group should be smaller than world
      if (np_group >= np_world) test_ok = .false.
      if (np_group /= world_comm%size()/2) test_ok = .false.

      ! Switch back to world
      call ddi_scope(DDI_COMM_WORLD)
      work_comm = get_working_comm()
      if (work_comm%size() /= np_world) test_ok = .false.

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Scope switching"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Scope switching"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_scope_switch

   subroutine test_scope_array_isolation()
      integer(int32) :: handle1, handle2
      integer(int32) :: ilo, ihi, jlo, jhi
      integer(int32) :: nrows, ncols, i, j
      real(dp), allocatable :: buffer(:), verify(:)
      type(comm_t) :: group_comm
      logical :: test_ok
      real(dp) :: expected_val

      if (world_comm%leader()) print *, "Test 4: Array isolation by scope..."

      nrows = 10
      ncols = 4

      ! Create array in world scope
      call ddi_scope(DDI_COMM_WORLD)
      call darray_create(nrows, ncols, handle1, 0.0_dp)

      ! Fill with world rank values
      call darray_distrib(handle1, world_comm%rank(), ilo, ihi, jlo, jhi)
      if (jhi >= jlo) then
         allocate (buffer(nrows*(jhi - jlo + 1)))
         buffer = real(world_comm%rank(), dp)
         call darray_put(handle1, 1, nrows, jlo, jhi, buffer)
         deallocate (buffer)
      end if

      call world_comm%barrier()

      ! Switch to group scope and create a new array
      call ddi_scope(DDI_COMM_GROUP)
      group_comm = get_working_comm()

      call darray_create(nrows, ncols, handle2, 0.0_dp)

      ! Fill with group rank values (should be different from world rank)
      call darray_distrib(handle2, group_comm%rank(), ilo, ihi, jlo, jhi)
      if (jhi >= jlo) then
         allocate (buffer(nrows*(jhi - jlo + 1)))
         buffer = real(group_comm%rank() + 100, dp)  ! Offset to distinguish
         call darray_put(handle2, 1, nrows, jlo, jhi, buffer)
         deallocate (buffer)
      end if

      call group_comm%barrier()

      ! Verify arrays are independent
      test_ok = .true.

      ! Check world array still has world rank values
      call ddi_scope(DDI_COMM_WORLD)
      if (world_comm%leader()) then
         allocate (verify(nrows*ncols))
         call darray_get(handle1, 1, nrows, 1, ncols, verify)

         ! First column should have rank 0's value
         if (abs(verify(1) - 0.0_dp) > 1.0e-10_dp) then
            test_ok = .false.
         end if
         deallocate (verify)
      end if

      ! Check group array has group rank values
      call ddi_scope(DDI_COMM_GROUP)
      if (group_comm%leader()) then
         allocate (verify(nrows*ncols))
         call darray_get(handle2, 1, nrows, 1, ncols, verify)

         ! First element should be 100 (group rank 0 + 100)
         if (abs(verify(1) - 100.0_dp) > 1.0e-10_dp) then
            test_ok = .false.
         end if
         deallocate (verify)
      end if

      ! Cleanup
      call ddi_scope(DDI_COMM_GROUP)
      call darray_destroy(handle2)
      call ddi_scope(DDI_COMM_WORLD)
      call darray_destroy(handle1)

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Array isolation by scope"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Array isolation by scope"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_scope_array_isolation

end program test_groups
