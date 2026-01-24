!> Test MPI RMA Accumulate operation
!!
!! Tests atomic accumulate for concurrent updates
!! This pattern is useful for Fock matrix accumulation
!! Run with: mpirun -np 4 ./test_rma_accumulate
program test_rma_accumulate
   use pic_mpi_lib
   use pic_types, only: dp, int32
   use mpi_f08, only: MPI_ADDRESS_KIND
   implicit none

   type(comm_t) :: world_comm
   integer :: n_passed, n_failed

   n_passed = 0
   n_failed = 0

   ! Initialize MPI (no args = default MPI_THREAD_FUNNELED)
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
      print *, "Testing RMA Accumulate Operations"
      print *, "========================================"
      print *, "Number of ranks:", world_comm%size()
   end if

   call world_comm%barrier()

   ! Test 1: Simple accumulate (all workers add to rank 0)
   call test_simple_accumulate()

   ! Test 2: Multiple accumulates from same source
   call test_multiple_accumulates()

   ! Test 3: Concurrent accumulates from all ranks
   call test_concurrent_accumulates()

   call world_comm%barrier()

   ! Summary
   if (world_comm%leader()) then
      print *, "========================================"
      print *, "Results: ", n_passed, " passed, ", n_failed, " failed"
      print *, "========================================"
   end if

   call world_comm%finalize()
   call pic_mpi_finalize()

contains

   subroutine test_simple_accumulate()
      integer, parameter :: LENGTH = 10
      real(dp), pointer :: buffer(:)
      real(dp) :: contrib(LENGTH)
      type(win_t) :: win
      integer(MPI_ADDRESS_KIND) :: disp
      integer :: i
      logical :: test_ok
      real(dp) :: expected

      if (world_comm%leader()) print *, "Test 1: Simple accumulate (workers -> rank 0)..."

      call win_allocate(world_comm, LENGTH, buffer, win)

      ! Initialize rank 0's buffer to zeros, workers prepare contributions
      if (world_comm%rank() == 0) then
         buffer = 0.0_dp
      else
         buffer = 0.0_dp
         contrib = [(real(world_comm%rank(), dp), i=1, LENGTH)]
      end if

      ! Synchronize with fence
      call win%fence()

      ! All workers accumulate to rank 0
      if (world_comm%rank() /= 0) then
         disp = 0_MPI_ADDRESS_KIND
         call win%accumulate_dp(0, disp, LENGTH, contrib)
      end if

      ! Synchronize after accumulate
      call win%fence()

      ! Verify on rank 0
      ! Expected: sum of ranks 1 to N-1 for each element
      test_ok = .true.
      if (world_comm%rank() == 0) then
         expected = 0.0_dp
         do i = 1, world_comm%size() - 1
            expected = expected + real(i, dp)
         end do

         do i = 1, LENGTH
            if (abs(buffer(i) - expected) > 1.0e-10_dp) then
               test_ok = .false.
               print *, "  Mismatch at", i, ": expected", expected, "got", buffer(i)
               exit
            end if
         end do
      end if

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Simple accumulate (sum =", expected, ")"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Simple accumulate"
            n_failed = n_failed + 1
         end if
      end if

      call win%finalize()
   end subroutine test_simple_accumulate

   subroutine test_multiple_accumulates()
      integer, parameter :: LENGTH = 5
      real(dp), pointer :: buffer(:)
      real(dp) :: contrib(LENGTH)
      type(win_t) :: win
      integer(MPI_ADDRESS_KIND) :: disp
      integer :: i, iter
      integer, parameter :: NUM_ITERS = 10
      logical :: test_ok
      real(dp) :: expected

      if (world_comm%leader()) print *, "Test 2: Multiple accumulates from rank 1..."

      call win_allocate(world_comm, LENGTH, buffer, win)

      buffer = 0.0_dp
      contrib = 1.0_dp  ! Each accumulate adds 1.0

      call win%fence()

      ! Rank 1 does multiple accumulates to rank 0
      if (world_comm%rank() == 1) then
         disp = 0_MPI_ADDRESS_KIND
         do iter = 1, NUM_ITERS
            call win%accumulate_dp(0, disp, LENGTH, contrib)
         end do
      end if

      call win%fence()

      ! Verify on rank 0: each element should be NUM_ITERS
      test_ok = .true.
      if (world_comm%rank() == 0) then
         expected = real(NUM_ITERS, dp)
         do i = 1, LENGTH
            if (abs(buffer(i) - expected) > 1.0e-10_dp) then
               test_ok = .false.
               print *, "  Mismatch at", i, ": expected", expected, "got", buffer(i)
               exit
            end if
         end do
      end if

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Multiple accumulates"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Multiple accumulates"
            n_failed = n_failed + 1
         end if
      end if

      call win%finalize()
   end subroutine test_multiple_accumulates

   subroutine test_concurrent_accumulates()
      integer, parameter :: LENGTH = 100
      real(dp), pointer :: buffer(:)
      real(dp) :: contrib(LENGTH)
      type(win_t) :: win
      integer(MPI_ADDRESS_KIND) :: disp
      integer :: i
      logical :: test_ok
      real(dp) :: expected, sum_ranks

      if (world_comm%leader()) print *, "Test 3: Concurrent accumulates (all -> rank 0)..."

      call win_allocate(world_comm, LENGTH, buffer, win)

      ! Initialize
      buffer = 0.0_dp
      ! Each rank contributes its rank number multiplied by element index
      contrib = [(real(world_comm%rank()*i, dp), i=1, LENGTH)]

      call win%fence()

      ! All ranks (including rank 0) accumulate to rank 0
      disp = 0_MPI_ADDRESS_KIND
      call win%accumulate_dp(0, disp, LENGTH, contrib)

      call win%fence()

      ! Verify on rank 0
      ! Expected for element i: i * sum(0, 1, 2, ..., N-1) = i * N*(N-1)/2
      test_ok = .true.
      if (world_comm%rank() == 0) then
         sum_ranks = real(world_comm%size()*(world_comm%size() - 1)/2, dp)

         do i = 1, LENGTH
            expected = real(i, dp)*sum_ranks
            if (abs(buffer(i) - expected) > 1.0e-8_dp) then
               test_ok = .false.
               print *, "  Mismatch at", i, ": expected", expected, "got", buffer(i)
               exit
            end if
         end do
      end if

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Concurrent accumulates"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Concurrent accumulates"
            n_failed = n_failed + 1
         end if
      end if

      call win%finalize()
   end subroutine test_concurrent_accumulates

end program test_rma_accumulate
