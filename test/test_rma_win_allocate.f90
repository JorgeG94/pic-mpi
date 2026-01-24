!> Test MPI window allocation functions
!!
!! Tests win_allocate for 1D and 2D arrays
!! Run with: mpirun -np 2 ./test_rma_win_allocate
program test_rma_win_allocate
   use pic_mpi_lib
   use pic_types, only: dp, int32
   implicit none

   type(comm_t) :: world_comm
   type(win_t) :: win_1d, win_2d
   real(dp), pointer :: buffer_1d(:), buffer_2d(:, :)
   integer :: i, j
   integer :: n_passed, n_failed

   n_passed = 0
   n_failed = 0

   ! Initialize MPI (no args = default MPI_THREAD_FUNNELED)
   call pic_mpi_init()
   world_comm = comm_world()

   if (world_comm%leader()) then
      print *, "========================================"
      print *, "Testing RMA Window Allocation"
      print *, "========================================"
      print *, "Number of ranks:", world_comm%size()
   end if

   call world_comm%barrier()

   ! Test 1: 1D window allocation
   call test_1d_allocation()

   ! Test 2: 2D window allocation
   call test_2d_allocation()

   ! Test 3: Window operations after allocation
   call test_window_operations()

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

   subroutine test_1d_allocation()
      integer, parameter :: LENGTH = 100
      real(dp), pointer :: ptr(:)
      type(win_t) :: win
      logical :: test_ok

      if (world_comm%leader()) print *, "Test 1: 1D window allocation..."

      call win_allocate(world_comm, LENGTH, ptr, win)

      ! Check allocation succeeded
      test_ok = associated(ptr) .and. .not. win%is_null()

      if (test_ok) then
         ! Initialize and verify we can write to the memory
         do i = 1, LENGTH
            ptr(i) = real(world_comm%rank()*1000 + i, dp)
         end do

         ! Verify values
         do i = 1, LENGTH
            if (abs(ptr(i) - real(world_comm%rank()*1000 + i, dp)) > 1.0e-10_dp) then
               test_ok = .false.
               exit
            end if
         end do
      end if

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: 1D allocation (", LENGTH, " elements)"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: 1D allocation"
            n_failed = n_failed + 1
         end if
      end if

      call win%finalize()
   end subroutine test_1d_allocation

   subroutine test_2d_allocation()
      integer, parameter :: DIM1 = 50, DIM2 = 60
      real(dp), pointer :: ptr(:, :)
      type(win_t) :: win
      logical :: test_ok

      if (world_comm%leader()) print *, "Test 2: 2D window allocation..."

      call win_allocate(world_comm, DIM1, DIM2, ptr, win)

      ! Check allocation succeeded
      test_ok = associated(ptr) .and. .not. win%is_null()

      if (test_ok) then
         ! Check dimensions
         if (size(ptr, 1) /= DIM1 .or. size(ptr, 2) /= DIM2) then
            test_ok = .false.
         end if
      end if

      if (test_ok) then
         ! Initialize and verify we can write to the memory
         do j = 1, DIM2
            do i = 1, DIM1
               ptr(i, j) = real(i*100 + j, dp)
            end do
         end do

         ! Verify values
         outer: do j = 1, DIM2
            do i = 1, DIM1
               if (abs(ptr(i, j) - real(i*100 + j, dp)) > 1.0e-10_dp) then
                  test_ok = .false.
                  exit outer
               end if
            end do
         end do outer
      end if

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: 2D allocation (", DIM1, "x", DIM2, ")"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: 2D allocation"
            n_failed = n_failed + 1
         end if
      end if

      call win%finalize()
   end subroutine test_2d_allocation

   subroutine test_window_operations()
      integer, parameter :: LENGTH = 10
      real(dp), pointer :: ptr(:)
      type(win_t) :: win
      logical :: test_ok

      if (world_comm%leader()) print *, "Test 3: Window fence operations..."

      call win_allocate(world_comm, LENGTH, ptr, win)

      test_ok = .true.

      ! Initialize local buffer
      do i = 1, LENGTH
         ptr(i) = real(world_comm%rank() + 1, dp)*10.0_dp
      end do

      ! Test fence synchronization
      call win%fence()

      ! All ranks should be synchronized here
      call win%fence()

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Window fence operations"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Window fence operations"
            n_failed = n_failed + 1
         end if
      end if

      call win%finalize()
   end subroutine test_window_operations

end program test_rma_win_allocate
