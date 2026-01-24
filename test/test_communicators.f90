!> Test MPI communicator operations
!!
!! Tests split, split_by, duplicate, discard_leader, discard_to
!! Run with: mpirun -np 4 ./test_communicators
program test_communicators
   use pic_mpi_lib
   use pic_types, only: dp, int32
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
      print *, "Testing MPI Communicator Operations"
      print *, "========================================"
      print *, "Number of ranks:", world_comm%size()
   end if

   call world_comm%barrier()

   call test_duplicate()
   call test_split_by_color()
   call test_split_shared()
   call test_discard_leader()
   call test_discard_to()

   call world_comm%barrier()

   if (world_comm%leader()) then
      print *, "========================================"
      print *, "Results: ", n_passed, " passed, ", n_failed, " failed"
      print *, "========================================"
   end if

   call world_comm%finalize()
   call pic_mpi_finalize()

contains

   subroutine test_duplicate()
      type(comm_t) :: dup_comm
      logical :: test_ok

      if (world_comm%leader()) print *, "Test 1: Communicator duplicate..."

      test_ok = .true.

      dup_comm = world_comm%duplicate()

      ! Verify properties match
      if (dup_comm%is_null()) then
         test_ok = .false.
         print *, "  Rank", world_comm%rank(), ": duplicate returned null"
      else if (dup_comm%size() /= world_comm%size()) then
         test_ok = .false.
         print *, "  Rank", world_comm%rank(), ": size mismatch"
      else if (dup_comm%rank() /= world_comm%rank()) then
         test_ok = .false.
         print *, "  Rank", world_comm%rank(), ": rank mismatch"
      end if

      call dup_comm%finalize()

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Duplicate"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Duplicate"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_duplicate

   subroutine test_split_by_color()
      type(comm_t) :: split_comm
      integer :: color
      logical :: test_ok

      if (world_comm%leader()) print *, "Test 2: Split by color..."

      test_ok = .true.

      ! Split into two groups: even and odd ranks
      if (mod(world_comm%rank(), 2) == 0) then
         color = 0  ! Even ranks
      else
         color = 1  ! Odd ranks
      end if

      split_comm = world_comm%split_by(color)

      if (split_comm%is_null()) then
         test_ok = .false.
         print *, "  Rank", world_comm%rank(), ": split returned null"
      else
         ! Verify size is approximately half
         ! For N ranks: evens get ceil(N/2), odds get floor(N/2)
         if (color == 0) then
            if (split_comm%size() /= (world_comm%size() + 1)/2) then
               test_ok = .false.
               print *, "  Rank", world_comm%rank(), ": even group size wrong"
            end if
         else
            if (split_comm%size() /= world_comm%size()/2) then
               test_ok = .false.
               print *, "  Rank", world_comm%rank(), ": odd group size wrong"
            end if
         end if
      end if

      if (.not. split_comm%is_null()) call split_comm%finalize()

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Split by color"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Split by color"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_split_by_color

   subroutine test_split_shared()
      type(comm_t) :: shared_comm
      logical :: test_ok

      if (world_comm%leader()) print *, "Test 3: Split shared memory..."

      test_ok = .true.

      ! Split into shared memory groups (same node)
      shared_comm = world_comm%split()

      if (shared_comm%is_null()) then
         test_ok = .false.
         print *, "  Rank", world_comm%rank(), ": shared split returned null"
      else
         ! On a single node, all ranks should be in same communicator
         ! Verify rank is valid and within size
         if (shared_comm%rank() < 0 .or. shared_comm%rank() >= shared_comm%size()) then
            test_ok = .false.
            print *, "  Rank", world_comm%rank(), ": invalid local rank"
         end if
      end if

      if (.not. shared_comm%is_null()) call shared_comm%finalize()

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Split shared memory"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Split shared memory"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_split_shared

   subroutine test_discard_leader()
      type(comm_t) :: no_leader_comm
      logical :: test_ok

      if (world_comm%leader()) print *, "Test 4: Discard leader..."

      test_ok = .true.

      no_leader_comm = world_comm%discard_leader()

      if (world_comm%rank() == 0) then
         ! Leader should get null communicator
         if (.not. no_leader_comm%is_null()) then
            test_ok = .false.
            print *, "  Leader got non-null communicator"
         end if
      else
         ! Non-leaders should get valid communicator
         if (no_leader_comm%is_null()) then
            test_ok = .false.
            print *, "  Rank", world_comm%rank(), ": got null communicator"
         else if (no_leader_comm%size() /= world_comm%size() - 1) then
            test_ok = .false.
            print *, "  Rank", world_comm%rank(), ": wrong size", &
               no_leader_comm%size(), "expected", world_comm%size() - 1
         end if
      end if

      if (.not. no_leader_comm%is_null()) call no_leader_comm%finalize()

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Discard leader"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Discard leader"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_discard_leader

   subroutine test_discard_to()
      type(comm_t) :: subset_comm
      integer, parameter :: KEEP_RANKS = 2
      logical :: test_ok

      if (world_comm%leader()) print *, "Test 5: Discard to (keep first 2)..."

      test_ok = .true.

      subset_comm = world_comm%discard_to(KEEP_RANKS)

      if (world_comm%rank() < KEEP_RANKS) then
         ! First N ranks should get valid communicator
         if (subset_comm%is_null()) then
            test_ok = .false.
            print *, "  Rank", world_comm%rank(), ": got null (should be valid)"
         else if (subset_comm%size() /= KEEP_RANKS) then
            test_ok = .false.
            print *, "  Rank", world_comm%rank(), ": wrong size"
         end if
      else
         ! Other ranks should get null
         if (.not. subset_comm%is_null()) then
            test_ok = .false.
            print *, "  Rank", world_comm%rank(), ": got valid (should be null)"
         end if
      end if

      if (.not. subset_comm%is_null()) call subset_comm%finalize()

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Discard to"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Discard to"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_discard_to

end program test_communicators
