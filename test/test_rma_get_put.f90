!> Test MPI RMA Get/Put operations
!!
!! Tests blocking get_dp and put_dp with fence synchronization
!! Run with: mpirun -np 2 ./test_rma_get_put
program test_rma_get_put
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
      print *, "Testing RMA Get/Put Operations"
      print *, "========================================"
      print *, "Number of ranks:", world_comm%size()
   end if

   call world_comm%barrier()

   ! Test 1: Simple put from rank 0 to rank 1
   call test_simple_put()

   ! Test 2: Simple get from rank 1 to rank 0
   call test_simple_get()

   ! Test 3: Bidirectional exchange
   call test_bidirectional()

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

   subroutine test_simple_put()
      integer, parameter :: LENGTH = 10
      real(dp), pointer :: buffer(:)
      real(dp) :: send_data(LENGTH)
      type(win_t) :: win
      integer(MPI_ADDRESS_KIND) :: disp
      integer :: i
      logical :: test_ok

      if (world_comm%leader()) print *, "Test 1: Simple Put (rank 0 -> rank 1)..."

      ! Allocate window on all ranks
      call win_allocate(world_comm, LENGTH, buffer, win)

      ! Initialize
      buffer = 0.0_dp
      if (world_comm%rank() == 0) then
         send_data = [(real(i*10, dp), i=1, LENGTH)]
      end if

      ! Synchronize before RMA
      call win%fence()

      ! Rank 0 puts data into rank 1's window
      if (world_comm%rank() == 0) then
         disp = 0_MPI_ADDRESS_KIND
         call win%put_dp(1, disp, LENGTH, send_data)
      end if

      ! Synchronize after RMA
      call win%fence()

      ! Verify on rank 1
      test_ok = .true.
      if (world_comm%rank() == 1) then
         do i = 1, LENGTH
            if (abs(buffer(i) - real(i*10, dp)) > 1.0e-10_dp) then
               test_ok = .false.
               print *, "  Mismatch at ", i, ": expected", real(i*10, dp), "got", buffer(i)
               exit
            end if
         end do
      end if

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Simple Put"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Simple Put"
            n_failed = n_failed + 1
         end if
      end if

      call win%finalize()
   end subroutine test_simple_put

   subroutine test_simple_get()
      integer, parameter :: LENGTH = 10
      real(dp), pointer :: buffer(:)
      real(dp) :: recv_data(LENGTH)
      type(win_t) :: win
      integer(MPI_ADDRESS_KIND) :: disp
      integer :: i
      logical :: test_ok

      if (world_comm%leader()) print *, "Test 2: Simple Get (rank 0 <- rank 1)..."

      ! Allocate window on all ranks
      call win_allocate(world_comm, LENGTH, buffer, win)

      ! Initialize - rank 1 has the source data
      if (world_comm%rank() == 1) then
         buffer = [(real(i*5, dp), i=1, LENGTH)]
      else
         buffer = 0.0_dp
      end if
      recv_data = 0.0_dp

      ! Synchronize before RMA
      call win%fence()

      ! Rank 0 gets data from rank 1's window
      if (world_comm%rank() == 0) then
         disp = 0_MPI_ADDRESS_KIND
         call win%get_dp(1, disp, LENGTH, recv_data)
      end if

      ! Synchronize after RMA
      call win%fence()

      ! Verify on rank 0
      test_ok = .true.
      if (world_comm%rank() == 0) then
         do i = 1, LENGTH
            if (abs(recv_data(i) - real(i*5, dp)) > 1.0e-10_dp) then
               test_ok = .false.
               print *, "  Mismatch at ", i, ": expected", real(i*5, dp), "got", recv_data(i)
               exit
            end if
         end do
      end if

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Simple Get"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Simple Get"
            n_failed = n_failed + 1
         end if
      end if

      call win%finalize()
   end subroutine test_simple_get

   subroutine test_bidirectional()
      integer, parameter :: LENGTH = 20
      real(dp), pointer :: buffer(:)
      real(dp) :: recv_data(LENGTH)
      type(win_t) :: win
      integer(MPI_ADDRESS_KIND) :: disp
      integer :: i, partner
      logical :: test_ok
      real(dp) :: expected

      if (world_comm%leader()) print *, "Test 3: Bidirectional exchange..."

      ! Allocate window on all ranks
      call win_allocate(world_comm, LENGTH, buffer, win)

      ! Initialize with rank-specific data
      buffer = [(real(world_comm%rank()*100 + i, dp), i=1, LENGTH)]
      recv_data = 0.0_dp

      ! Synchronize before RMA
      call win%fence()

      ! Each rank gets data from its partner
      if (world_comm%rank() == 0) then
         partner = 1
      else if (world_comm%rank() == 1) then
         partner = 0
      else
         partner = -1
      end if

      if (partner >= 0) then
         disp = 0_MPI_ADDRESS_KIND
         call win%get_dp(partner, disp, LENGTH, recv_data)
      end if

      ! Synchronize after RMA
      call win%fence()

      ! Verify
      test_ok = .true.
      if (partner >= 0) then
         do i = 1, LENGTH
            expected = real(partner*100 + i, dp)
            if (abs(recv_data(i) - expected) > 1.0e-10_dp) then
               test_ok = .false.
               print *, "  Rank", world_comm%rank(), "mismatch at", i, &
                  ": expected", expected, "got", recv_data(i)
               exit
            end if
         end do
      end if

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Bidirectional exchange"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Bidirectional exchange"
            n_failed = n_failed + 1
         end if
      end if

      call win%finalize()
   end subroutine test_bidirectional

end program test_rma_get_put
