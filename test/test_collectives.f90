!> Test MPI collective operations
!!
!! Tests allreduce, allgather, bcast, barrier
!! Run with: mpirun -np 4 ./test_collectives
program test_collectives
   use pic_mpi_lib
   use pic_types, only: dp, int32, int64
   implicit none

   type(comm_t) :: world_comm
   integer :: n_passed, n_failed

   n_passed = 0
   n_failed = 0

   call pic_mpi_init()
   world_comm = comm_world()

   if (world_comm%size() < 2) then
      if (world_comm%leader()) then
         print *, "ERROR: This test requires at least 2 MPI ranks"
      end if
      call pic_mpi_finalize()
      stop 1
   end if

   if (world_comm%leader()) then
      print *, "========================================"
      print *, "Testing MPI Collective Operations"
      print *, "========================================"
      print *, "Number of ranks:", world_comm%size()
   end if

   call world_comm%barrier()

   call test_barrier()
   call test_allreduce_scalar()
   call test_allreduce_array()
   call test_allgather()
   call test_bcast()

   call world_comm%barrier()

   if (world_comm%leader()) then
      print *, "========================================"
      print *, "Results: ", n_passed, " passed, ", n_failed, " failed"
      print *, "========================================"
   end if

   call world_comm%finalize()
   call pic_mpi_finalize()

contains

   subroutine test_barrier()
      integer :: i
      logical :: test_ok

      if (world_comm%leader()) print *, "Test 1: Barrier synchronization..."

      test_ok = .true.

      ! Simple barrier test - all ranks should reach this point
      do i = 1, 10
         call world_comm%barrier()
      end do

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Barrier"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Barrier"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_barrier

   subroutine test_allreduce_scalar()
      real(dp) :: value, expected
      integer(int32) :: ivalue, iexpected
      logical :: test_ok

      if (world_comm%leader()) print *, "Test 2: Allreduce scalar..."

      test_ok = .true.

      ! Double precision sum
      value = real(world_comm%rank() + 1, dp)
      call allreduce(world_comm, value)

      ! Expected: 1 + 2 + ... + N = N*(N+1)/2
      expected = real(world_comm%size()*(world_comm%size() + 1)/2, dp)
      if (abs(value - expected) > 1.0e-10_dp) then
         test_ok = .false.
         print *, "  Rank", world_comm%rank(), "DP mismatch: expected", expected, "got", value
      end if

      ! Integer sum
      ivalue = world_comm%rank() + 1
      call allreduce(world_comm, ivalue)

      iexpected = world_comm%size()*(world_comm%size() + 1)/2
      if (ivalue /= iexpected) then
         test_ok = .false.
         print *, "  Rank", world_comm%rank(), "INT mismatch: expected", iexpected, "got", ivalue
      end if

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Allreduce scalar (sum =", expected, ")"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Allreduce scalar"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_allreduce_scalar

   subroutine test_allreduce_array()
      integer, parameter :: LENGTH = 100
      real(dp) :: values(LENGTH), expected
      integer :: i
      logical :: test_ok

      if (world_comm%leader()) print *, "Test 3: Allreduce array..."

      test_ok = .true.

      ! Each rank contributes rank+1 to each element
      values = real(world_comm%rank() + 1, dp)

      call allreduce(world_comm, values)

      ! Expected: sum of all ranks' contributions
      expected = real(world_comm%size()*(world_comm%size() + 1)/2, dp)

      do i = 1, LENGTH
         if (abs(values(i) - expected) > 1.0e-10_dp) then
            test_ok = .false.
            print *, "  Rank", world_comm%rank(), "mismatch at", i
            exit
         end if
      end do

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Allreduce array"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Allreduce array"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_allreduce_array

   subroutine test_allgather()
      integer :: send_val
      integer, allocatable :: recv_vals(:)
      integer :: i
      logical :: test_ok

      if (world_comm%leader()) print *, "Test 4: Allgather..."

      test_ok = .true.

      send_val = world_comm%rank()*10
      allocate (recv_vals(world_comm%size()))

      call allgather(world_comm, send_val, recv_vals)

      ! Verify: recv_vals(i) should be (i-1)*10
      do i = 1, world_comm%size()
         if (recv_vals(i) /= (i - 1)*10) then
            test_ok = .false.
            print *, "  Rank", world_comm%rank(), "mismatch at", i, &
               ": expected", (i - 1)*10, "got", recv_vals(i)
            exit
         end if
      end do

      deallocate (recv_vals)

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Allgather"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Allgather"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_allgather

   subroutine test_bcast()
      integer :: value
      integer(int64) :: value64
      logical :: test_ok

      if (world_comm%leader()) print *, "Test 5: Broadcast..."

      test_ok = .true.

      ! Integer broadcast from root
      if (world_comm%rank() == 0) then
         value = 42
      else
         value = 0
      end if

      call bcast(world_comm, value, 1, 0)

      if (value /= 42) then
         test_ok = .false.
         print *, "  Rank", world_comm%rank(), "INT32 mismatch: expected 42, got", value
      end if

      ! Integer64 broadcast
      if (world_comm%rank() == 0) then
         value64 = 123456789012_int64
      else
         value64 = 0_int64
      end if

      call bcast(world_comm, value64, 1, 0)

      if (value64 /= 123456789012_int64) then
         test_ok = .false.
         print *, "  Rank", world_comm%rank(), "INT64 mismatch"
      end if

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Broadcast"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Broadcast"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_bcast

end program test_collectives
