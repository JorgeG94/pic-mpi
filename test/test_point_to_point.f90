!> Test MPI point-to-point operations
!!
!! Tests send/recv, isend/irecv for various data types
!! Run with: mpirun -np 2 ./test_point_to_point
program test_point_to_point
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
      print *, "Testing MPI Point-to-Point Operations"
      print *, "========================================"
      print *, "Number of ranks:", world_comm%size()
   end if

   call world_comm%barrier()

   call test_send_recv_int32()
   call test_send_recv_int64()
   call test_send_recv_real_dp()
   call test_send_recv_logical()
   call test_send_recv_arrays()
   call test_isend_irecv()
   call test_iprobe()

   call world_comm%barrier()

   if (world_comm%leader()) then
      print *, "========================================"
      print *, "Results: ", n_passed, " passed, ", n_failed, " failed"
      print *, "========================================"
   end if

   call world_comm%finalize()
   call pic_mpi_finalize()

contains

   subroutine test_send_recv_int32()
      integer(int32) :: send_val, recv_val
      logical :: test_ok

      if (world_comm%leader()) print *, "Test 1: Send/Recv int32..."

      test_ok = .true.

      if (world_comm%rank() == 0) then
         send_val = 12345
         call send(world_comm, send_val, 1, 100)
         call recv(world_comm, recv_val, 1, 101)
         if (recv_val /= 54321) test_ok = .false.
      else if (world_comm%rank() == 1) then
         call recv(world_comm, recv_val, 0, 100)
         if (recv_val /= 12345) test_ok = .false.
         send_val = 54321
         call send(world_comm, send_val, 0, 101)
      end if

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Send/Recv int32"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Send/Recv int32"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_send_recv_int32

   subroutine test_send_recv_int64()
      integer(int64) :: send_val, recv_val
      logical :: test_ok

      if (world_comm%leader()) print *, "Test 2: Send/Recv int64..."

      test_ok = .true.

      if (world_comm%rank() == 0) then
         send_val = 9876543210_int64
         call send(world_comm, send_val, 1, 200)
         call recv(world_comm, recv_val, 1, 201)
         if (recv_val /= 1234567890_int64) test_ok = .false.
      else if (world_comm%rank() == 1) then
         call recv(world_comm, recv_val, 0, 200)
         if (recv_val /= 9876543210_int64) test_ok = .false.
         send_val = 1234567890_int64
         call send(world_comm, send_val, 0, 201)
      end if

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Send/Recv int64"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Send/Recv int64"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_send_recv_int64

   subroutine test_send_recv_real_dp()
      real(dp) :: send_val, recv_val
      logical :: test_ok

      if (world_comm%leader()) print *, "Test 3: Send/Recv real(dp)..."

      test_ok = .true.

      if (world_comm%rank() == 0) then
         send_val = 3.14159265358979_dp
         call send(world_comm, send_val, 1, 300)
         call recv(world_comm, recv_val, 1, 301)
         if (abs(recv_val - 2.71828182845905_dp) > 1.0e-14_dp) test_ok = .false.
      else if (world_comm%rank() == 1) then
         call recv(world_comm, recv_val, 0, 300)
         if (abs(recv_val - 3.14159265358979_dp) > 1.0e-14_dp) test_ok = .false.
         send_val = 2.71828182845905_dp
         call send(world_comm, send_val, 0, 301)
      end if

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Send/Recv real(dp)"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Send/Recv real(dp)"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_send_recv_real_dp

   subroutine test_send_recv_logical()
      logical :: send_val, recv_val
      logical :: test_ok

      if (world_comm%leader()) print *, "Test 4: Send/Recv logical..."

      test_ok = .true.

      if (world_comm%rank() == 0) then
         send_val = .true.
         call send(world_comm, send_val, 1, 400)
         call recv(world_comm, recv_val, 1, 401)
         if (recv_val .neqv. .false.) test_ok = .false.
      else if (world_comm%rank() == 1) then
         call recv(world_comm, recv_val, 0, 400)
         if (recv_val .neqv. .true.) test_ok = .false.
         send_val = .false.
         call send(world_comm, send_val, 0, 401)
      end if

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Send/Recv logical"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Send/Recv logical"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_send_recv_logical

   subroutine test_send_recv_arrays()
      integer, parameter :: LENGTH = 50
      integer(int32) :: send_arr(LENGTH)
      integer(int32), allocatable :: recv_arr(:)
      real(dp) :: send_dp(LENGTH), recv_dp(LENGTH)
      type(MPI_Status) :: status
      integer :: i
      logical :: test_ok

      if (world_comm%leader()) print *, "Test 5: Send/Recv arrays..."

      test_ok = .true.

      if (world_comm%rank() == 0) then
         send_arr = [(i*2, i=1, LENGTH)]
         call send(world_comm, send_arr, 1, 500)

         call recv(world_comm, recv_arr, 1, 501, status)
         if (size(recv_arr) /= LENGTH) then
            test_ok = .false.
         else
            do i = 1, LENGTH
               if (recv_arr(i) /= i*3) then
                  test_ok = .false.
                  exit
               end if
            end do
         end if
      else if (world_comm%rank() == 1) then
         call recv(world_comm, recv_arr, 0, 500, status)
         if (size(recv_arr) /= LENGTH) then
            test_ok = .false.
         else
            do i = 1, LENGTH
               if (recv_arr(i) /= i*2) then
                  test_ok = .false.
                  exit
               end if
            end do
         end if

         send_arr = [(i*3, i=1, LENGTH)]
         call send(world_comm, send_arr, 0, 501)
      end if

      if (allocated(recv_arr)) deallocate (recv_arr)

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Send/Recv arrays"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Send/Recv arrays"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_send_recv_arrays

   subroutine test_isend_irecv()
      integer(int32) :: send_val, recv_val
      type(request_t) :: send_req, recv_req
      logical :: test_ok

      if (world_comm%leader()) print *, "Test 6: Isend/Irecv..."

      test_ok = .true.

      if (world_comm%rank() == 0) then
         send_val = 999
         recv_val = 0
         call isend(world_comm, send_val, 1, 600, send_req)
         call irecv(world_comm, recv_val, 1, 601, recv_req)
         call wait(send_req)
         call wait(recv_req)
         if (recv_val /= 888) test_ok = .false.
      else if (world_comm%rank() == 1) then
         send_val = 888
         recv_val = 0
         call irecv(world_comm, recv_val, 0, 600, recv_req)
         call isend(world_comm, send_val, 0, 601, send_req)
         call wait(recv_req)
         call wait(send_req)
         if (recv_val /= 999) test_ok = .false.
      end if

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Isend/Irecv"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Isend/Irecv"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_isend_irecv

   subroutine test_iprobe()
      integer(int32) :: send_val, recv_val
      type(MPI_Status) :: status
      logical :: has_message, test_ok
      integer :: attempts

      if (world_comm%leader()) print *, "Test 7: Iprobe..."

      test_ok = .true.

      if (world_comm%rank() == 0) then
         ! Wait for message from rank 1
         has_message = .false.
         attempts = 0
         do while (.not. has_message .and. attempts < 10000)
            call iprobe(world_comm, 1, 700, has_message, status)
            attempts = attempts + 1
         end do

         if (has_message) then
            call recv(world_comm, recv_val, 1, 700)
            if (recv_val /= 777) test_ok = .false.
         else
            test_ok = .false.
            print *, "  Iprobe did not detect message"
         end if
      else if (world_comm%rank() == 1) then
         send_val = 777
         call send(world_comm, send_val, 0, 700)
      end if

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Iprobe"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Iprobe"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_iprobe

end program test_point_to_point
