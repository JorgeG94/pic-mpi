!> Test MPI RMA Non-blocking operations (Rget/Rput)
!!
!! Tests rget_dp, rput_dp with lock_all/unlock_all/flush_all
!! This is the pattern needed for metalquicha hessian transfers
!! Run with: mpirun -np 4 ./test_rma_nonblocking
program test_rma_nonblocking
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
      print *, "Testing RMA Non-blocking Operations"
      print *, "========================================"
      print *, "Number of ranks:", world_comm%size()
   end if

   call world_comm%barrier()

   ! Test 1: Non-blocking get with lock/unlock
   call test_rget_lock()

   ! Test 2: Non-blocking put with lock/unlock
   call test_rput_lock()

   ! Test 3: lock_all/flush_all pattern (metalquicha pattern)
   call test_lock_all_pattern()

   ! Test 4: Multiple concurrent rgets
   call test_multiple_rgets()

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

   subroutine test_rget_lock()
      integer, parameter :: LENGTH = 100
      real(dp), pointer :: buffer(:)
      real(dp) :: recv_data(LENGTH)
      type(win_t) :: win
      type(request_t) :: req
      integer(MPI_ADDRESS_KIND) :: disp
      integer :: i
      logical :: test_ok

      if (world_comm%leader()) print *, "Test 1: Non-blocking Rget with lock/unlock..."

      call win_allocate(world_comm, LENGTH, buffer, win)

      ! Rank 1 initializes data
      if (world_comm%rank() == 1) then
         buffer = [(real(i*7, dp), i=1, LENGTH)]
      else
         buffer = 0.0_dp
      end if
      recv_data = -1.0_dp

      call world_comm%barrier()

      ! Rank 0 does non-blocking get from rank 1
      if (world_comm%rank() == 0) then
         disp = 0_MPI_ADDRESS_KIND
         call win%lock(1)  ! Lock rank 1's window
         call win%rget_dp(1, disp, LENGTH, recv_data, req)
         call wait(req)    ! Wait for completion
         call win%unlock(1)
      end if

      call world_comm%barrier()

      ! Verify
      test_ok = .true.
      if (world_comm%rank() == 0) then
         do i = 1, LENGTH
            if (abs(recv_data(i) - real(i*7, dp)) > 1.0e-10_dp) then
               test_ok = .false.
               print *, "  Mismatch at", i, ": expected", real(i*7, dp), "got", recv_data(i)
               exit
            end if
         end do
      end if

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Non-blocking Rget"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Non-blocking Rget"
            n_failed = n_failed + 1
         end if
      end if

      call win%finalize()
   end subroutine test_rget_lock

   subroutine test_rput_lock()
      integer, parameter :: LENGTH = 100
      real(dp), pointer :: buffer(:)
      real(dp) :: send_data(LENGTH)
      type(win_t) :: win
      type(request_t) :: req
      integer(MPI_ADDRESS_KIND) :: disp
      integer :: i
      logical :: test_ok

      if (world_comm%leader()) print *, "Test 2: Non-blocking Rput with lock/unlock..."

      call win_allocate(world_comm, LENGTH, buffer, win)

      buffer = 0.0_dp
      if (world_comm%rank() == 0) then
         send_data = [(real(i*3, dp), i=1, LENGTH)]
      end if

      call world_comm%barrier()

      ! Rank 0 does non-blocking put to rank 1
      if (world_comm%rank() == 0) then
         disp = 0_MPI_ADDRESS_KIND
         call win%lock(1)
         call win%rput_dp(1, disp, LENGTH, send_data, req)
         call wait(req)
         call win%unlock(1)
      end if

      call world_comm%barrier()

      ! Verify on rank 1
      test_ok = .true.
      if (world_comm%rank() == 1) then
         do i = 1, LENGTH
            if (abs(buffer(i) - real(i*3, dp)) > 1.0e-10_dp) then
               test_ok = .false.
               print *, "  Mismatch at", i, ": expected", real(i*3, dp), "got", buffer(i)
               exit
            end if
         end do
      end if

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Non-blocking Rput"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Non-blocking Rput"
            n_failed = n_failed + 1
         end if
      end if

      call win%finalize()
   end subroutine test_rput_lock

   subroutine test_lock_all_pattern()
      integer, parameter :: LENGTH = 50
      real(dp), pointer :: buffer(:)
      real(dp), allocatable :: recv_data(:, :)
      type(win_t) :: win
      type(request_t), allocatable :: requests(:)
      integer(MPI_ADDRESS_KIND) :: disp
      integer :: i, j, num_workers
      logical :: test_ok
      real(dp) :: expected

      if (world_comm%leader()) print *, "Test 3: lock_all/flush_all pattern..."

      call win_allocate(world_comm, LENGTH, buffer, win)

      ! Each rank initializes its buffer with rank-specific values
      buffer = [(real(world_comm%rank()*1000 + i, dp), i=1, LENGTH)]

      call world_comm%barrier()

      ! Coordinator (rank 0) fetches from all workers using lock_all
      test_ok = .true.
      if (world_comm%rank() == 0) then
         num_workers = world_comm%size() - 1
         if (num_workers > 0) then
            allocate (recv_data(LENGTH, num_workers))
            allocate (requests(num_workers))
            recv_data = -1.0_dp

            ! Lock all windows
            call win%lock_all()

            ! Issue non-blocking gets to all workers
            do i = 1, num_workers
               disp = 0_MPI_ADDRESS_KIND
               call win%rget_dp(i, disp, LENGTH, recv_data(:, i), requests(i))
            end do

            ! Flush to ensure completion
            call win%flush_all()

            ! Wait for all requests
            call waitall(requests)

            ! Unlock all
            call win%unlock_all()

            ! Verify data from each worker
            do i = 1, num_workers
               do j = 1, LENGTH
                  expected = real(i*1000 + j, dp)
                  if (abs(recv_data(j, i) - expected) > 1.0e-10_dp) then
                     test_ok = .false.
                     print *, "  Mismatch from worker", i, "at", j, &
                        ": expected", expected, "got", recv_data(j, i)
                     exit
                  end if
               end do
               if (.not. test_ok) exit
            end do

            deallocate (recv_data, requests)
         end if
      end if

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: lock_all/flush_all pattern"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: lock_all/flush_all pattern"
            n_failed = n_failed + 1
         end if
      end if

      call win%finalize()
   end subroutine test_lock_all_pattern

   subroutine test_multiple_rgets()
      integer, parameter :: LENGTH = 25
      real(dp), pointer :: buffer(:)
      real(dp) :: recv1(LENGTH), recv2(LENGTH)
      type(win_t) :: win
      type(request_t) :: req1, req2
      integer(MPI_ADDRESS_KIND) :: disp
      integer :: i
      logical :: test_ok

      if (world_comm%size() < 3) then
         if (world_comm%leader()) then
            print *, "Test 4: Skipped (needs 3+ ranks)"
         end if
         return
      end if

      if (world_comm%leader()) print *, "Test 4: Multiple concurrent Rgets..."

      call win_allocate(world_comm, LENGTH, buffer, win)

      ! Ranks 1 and 2 have different data
      if (world_comm%rank() == 1) then
         buffer = [(real(i*11, dp), i=1, LENGTH)]
      else if (world_comm%rank() == 2) then
         buffer = [(real(i*13, dp), i=1, LENGTH)]
      else
         buffer = 0.0_dp
      end if

      recv1 = -1.0_dp
      recv2 = -1.0_dp

      call world_comm%barrier()

      ! Rank 0 fetches from ranks 1 and 2 concurrently
      test_ok = .true.
      if (world_comm%rank() == 0) then
         call win%lock_all()

         disp = 0_MPI_ADDRESS_KIND
         call win%rget_dp(1, disp, LENGTH, recv1, req1)
         call win%rget_dp(2, disp, LENGTH, recv2, req2)

         ! Wait for both
         call wait(req1)
         call wait(req2)

         call win%unlock_all()

         ! Verify data from rank 1
         do i = 1, LENGTH
            if (abs(recv1(i) - real(i*11, dp)) > 1.0e-10_dp) then
               test_ok = .false.
               print *, "  Mismatch from rank 1 at", i
               exit
            end if
         end do

         ! Verify data from rank 2
         if (test_ok) then
            do i = 1, LENGTH
               if (abs(recv2(i) - real(i*13, dp)) > 1.0e-10_dp) then
                  test_ok = .false.
                  print *, "  Mismatch from rank 2 at", i
                  exit
               end if
            end do
         end if
      end if

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Multiple concurrent Rgets"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Multiple concurrent Rgets"
            n_failed = n_failed + 1
         end if
      end if

      call win%finalize()
   end subroutine test_multiple_rgets

end program test_rma_nonblocking
