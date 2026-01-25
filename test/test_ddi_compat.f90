!> Test DDI compatibility wrapper
!!
!! Tests DDI API functions for GAMESS compatibility.
!! Run with: mpirun -np 4 ./test_ddi_compat
program test_ddi_compat
   use pic_types, only: int32, int64, dp
   use ddi_compat
   implicit none

   integer :: n_passed, n_failed
   integer(int32) :: np, me

   n_passed = 0
   n_failed = 0

   call ddi_init()
   call ddi_nproc(np, me)

   if (np < 2) then
      if (me == 0) then
         print *, "ERROR: This test requires at least 2 MPI ranks"
      end if
      call ddi_finalize()
      stop 1
   end if

   if (me == 0) then
      print *, "========================================"
      print *, "Testing DDI Compatibility Layer"
      print *, "========================================"
      print *, "Number of ranks:", np
   end if

   call ddi_sync(0_int64)

   call test_nproc()
   call test_create_destroy()
   call test_get_put()
   call test_accumulate()
   call test_zero()
   call test_gsumf()
   call test_gsumi()
   call test_dlb()
   call test_groups()
   call test_nnode()
   call test_ngroup()
   call test_async_comm()
   call test_comm_ops()
   call test_comm_create_destroy()
   call test_arr_fill()

   call ddi_sync(0_int64)

   ! Summary
   if (me == 0) then
      print *, "========================================"
      print *, "Results: ", n_passed, " passed, ", n_failed, " failed"
      print *, "========================================"
   end if

   call ddi_finalize()

   if (n_failed > 0) stop 1

contains

   subroutine test_nproc()
      integer(int32) :: np_test, me_test
      logical :: test_ok

      if (me == 0) print *, "Test 1: ddi_nproc..."

      call ddi_nproc(np_test, me_test)
      test_ok = (np_test == np) .and. (me_test == me)

      call ddi_sync(0_int64)

      if (me == 0) then
         if (test_ok) then
            print *, "  PASS: ddi_nproc"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: ddi_nproc"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_nproc

   subroutine test_create_destroy()
      integer(int32) :: handle, ilo, ihi, jlo, jhi
      logical :: test_ok

      if (me == 0) print *, "Test 2: ddi_create/ddi_destroy..."

      call ddi_create(100, 50, handle, 0.0_dp)

      ! Check distribution
      call ddi_distrib(handle, me, ilo, ihi, jlo, jhi)

      test_ok = (ilo == 1) .and. (ihi == 100) .and. (jlo >= 1) .and. (jhi <= 50)

      call ddi_destroy(handle)

      call ddi_sync(0_int64)

      if (me == 0) then
         if (test_ok) then
            print *, "  PASS: ddi_create/ddi_destroy"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: ddi_create/ddi_destroy"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_create_destroy

   subroutine test_get_put()
      integer(int32) :: handle, ilo, ihi, jlo, jhi
      integer(int32) :: nrows, ncols, i, j
      real(dp), allocatable :: buffer(:), verify(:)
      logical :: test_ok

      if (me == 0) print *, "Test 3: ddi_get/ddi_put..."

      nrows = 10
      ncols = 8
      call ddi_create(nrows, ncols, handle, 0.0_dp)

      call ddi_distrib(handle, me, ilo, ihi, jlo, jhi)

      ! Put local data
      if (jhi >= jlo) then
         allocate (buffer(nrows*(jhi - jlo + 1)))
         do j = jlo, jhi
            do i = 1, nrows
               buffer(i + (j - jlo)*nrows) = real(i*100 + j, dp)
            end do
         end do
         call ddi_put(handle, 1, nrows, jlo, jhi, buffer)
         deallocate (buffer)
      end if

      call ddi_sync(0_int64)

      ! Rank 0 verifies
      test_ok = .true.
      if (me == 0) then
         allocate (verify(nrows*ncols))
         call ddi_get(handle, 1, nrows, 1, ncols, verify)

         do j = 1, ncols
            do i = 1, nrows
               if (abs(verify(i + (j - 1)*nrows) - real(i*100 + j, dp)) > 1.0e-10_dp) then
                  test_ok = .false.
                  exit
               end if
            end do
            if (.not. test_ok) exit
         end do
         deallocate (verify)
      end if

      call ddi_destroy(handle)
      call ddi_sync(0_int64)

      if (me == 0) then
         if (test_ok) then
            print *, "  PASS: ddi_get/ddi_put"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: ddi_get/ddi_put"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_get_put

   subroutine test_accumulate()
      integer(int32) :: handle, nrows, ncols, i
      real(dp), allocatable :: buffer(:), verify(:)
      real(dp) :: expected
      logical :: test_ok

      if (me == 0) print *, "Test 4: ddi_acc..."

      nrows = 5
      ncols = 4
      call ddi_create(nrows, ncols, handle, 0.0_dp)

      call ddi_sync(0_int64)

      ! All ranks accumulate 1.0 to all elements
      allocate (buffer(nrows*ncols))
      buffer = 1.0_dp
      call ddi_acc(handle, 1, nrows, 1, ncols, buffer)
      deallocate (buffer)

      call ddi_sync(0_int64)

      ! Rank 0 verifies sum equals number of ranks
      test_ok = .true.
      if (me == 0) then
         allocate (verify(nrows*ncols))
         call ddi_get(handle, 1, nrows, 1, ncols, verify)

         expected = real(np, dp)
         do i = 1, nrows*ncols
            if (abs(verify(i) - expected) > 1.0e-10_dp) then
               test_ok = .false.
               exit
            end if
         end do
         deallocate (verify)
      end if

      call ddi_destroy(handle)
      call ddi_sync(0_int64)

      if (me == 0) then
         if (test_ok) then
            print *, "  PASS: ddi_acc"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: ddi_acc"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_accumulate

   subroutine test_gsumf()
      real(dp) :: buffer(10), expected
      integer(int32) :: i
      logical :: test_ok

      if (me == 0) print *, "Test 5: ddi_gsumf..."

      ! Each rank contributes its rank value
      do i = 1, 10
         buffer(i) = real(me + i, dp)
      end do

      call ddi_gsumf(0, buffer, 10)

      ! Expected: sum of (rank + i) for all ranks
      test_ok = .true.
      do i = 1, 10
         ! Sum of (0+i) + (1+i) + ... + (np-1+i) = np*i + (0+1+...+np-1) = np*i + np*(np-1)/2
         expected = real(np*i + np*(np - 1)/2, dp)
         if (abs(buffer(i) - expected) > 1.0e-10_dp) then
            test_ok = .false.
            exit
         end if
      end do

      call ddi_sync(0_int64)

      if (me == 0) then
         if (test_ok) then
            print *, "  PASS: ddi_gsumf"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: ddi_gsumf"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_gsumf

   subroutine test_gsumi()
      integer(int64) :: buffer(10), expected  ! int64 for GAMESS -i8 compatibility
      integer(int32) :: i
      logical :: test_ok

      if (me == 0) print *, "Test 6: ddi_gsumi..."

      ! Each rank contributes its rank value
      do i = 1, 10
         buffer(i) = int(me + i, int64)
      end do

      call ddi_gsumi(0_int64, buffer, 10_int64)

      ! Expected: sum of (rank + i) for all ranks
      test_ok = .true.
      do i = 1, 10
         expected = int(np*i + np*(np - 1)/2, int64)
         if (buffer(i) /= expected) then
            test_ok = .false.
            exit
         end if
      end do

      call ddi_sync(0_int64)

      if (me == 0) then
         if (test_ok) then
            print *, "  PASS: ddi_gsumi"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: ddi_gsumi"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_gsumi

   subroutine test_dlb()
      integer(int64) :: counter, counters(100)
      integer(int64) :: total_work(1)  ! int64 for ddi_gsumi
      integer(int32) :: i, nwork
      logical :: test_ok

      if (me == 0) print *, "Test 7: ddi_dlbreset/ddi_dlbnext..."

      nwork = 20

      call ddi_dlbreset()

      ! Initialize counters
      do i = 1, 100
         counters(i) = -1_int64
      end do

      ! Each rank gets work items
      i = 1
      do
         call ddi_dlbnext(counter)
         if (counter >= nwork) exit
         counters(i) = counter
         i = i + 1
      end do

      call ddi_sync(0_int64)

      ! Verify total work done
      total_work(1) = int(i - 1, int64)
      call ddi_gsumi(0_int64, total_work, 1_int64)

      test_ok = .true.
      if (me == 0) then
         if (total_work(1) /= nwork) then
            test_ok = .false.
         end if
      end if

      call ddi_sync(0_int64)

      if (me == 0) then
         if (test_ok) then
            print *, "  PASS: ddi_dlbreset/ddi_dlbnext"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: ddi_dlbreset/ddi_dlbnext"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_dlb

   subroutine test_groups()
      integer(int32) :: world_id, group_id, master_id
      integer(int32) :: np_world, me_world, np_group, me_group
      logical :: test_ok

      if (me == 0) print *, "Test 8: ddi_group_create/ddi_scope..."

      ! Save world state
      call ddi_nproc(np_world, me_world)

      ! Create 2 groups (requires at least 4 ranks)
      if (np >= 4) then
         call ddi_group_create(2, world_id, group_id, master_id)

         test_ok = .true.

         ! Switch to group scope
         call ddi_scope(group_id)
         call ddi_nproc(np_group, me_group)

         ! Group should be half the size of world
         if (np_group /= np_world/2) test_ok = .false.

         ! Switch back to world
         call ddi_scope(world_id)
         call ddi_nproc(np_group, me_group)

         if (np_group /= np_world) test_ok = .false.
         if (me_group /= me_world) test_ok = .false.
      else
         test_ok = .true.  ! Skip if not enough ranks
         if (me == 0) print *, "  (Skipped - need 4 ranks for group test)"
      end if

      call ddi_sync(0_int64)

      if (me == 0) then
         if (test_ok) then
            print *, "  PASS: ddi_group_create/ddi_scope"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: ddi_group_create/ddi_scope"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_groups

   subroutine test_zero()
      integer(int32) :: handle, nrows, ncols, i
      real(dp), allocatable :: buffer(:), verify(:)
      logical :: test_ok

      if (me == 0) print *, "Test 9: ddi_zero..."

      nrows = 10
      ncols = 8
      call ddi_create(nrows, ncols, handle, 99.0_dp)

      ! Zero the array
      call ddi_zero(handle)

      call ddi_sync(0_int64)

      ! Rank 0 verifies all zeros
      test_ok = .true.
      if (me == 0) then
         allocate (verify(nrows*ncols))
         call ddi_get(handle, 1, nrows, 1, ncols, verify)

         do i = 1, nrows*ncols
            if (abs(verify(i)) > 1.0e-10_dp) then
               test_ok = .false.
               exit
            end if
         end do
         deallocate (verify)
      end if

      call ddi_destroy(handle)
      call ddi_sync(0_int64)

      if (me == 0) then
         if (test_ok) then
            print *, "  PASS: ddi_zero"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: ddi_zero"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_zero

   subroutine test_nnode()
      integer(int32) :: nnodes, mynode
      logical :: test_ok

      if (me == 0) print *, "Test 10: ddi_nnode..."

      call ddi_nnode(nnodes, mynode)

      test_ok = .true.
      ! Sanity checks: at least 1 node, mynode in valid range
      if (nnodes < 1) test_ok = .false.
      if (mynode < 0 .or. mynode >= nnodes) test_ok = .false.

      call ddi_sync(0_int64)

      if (me == 0) then
         if (test_ok) then
            print *, "  PASS: ddi_nnode (", nnodes, " nodes)"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: ddi_nnode"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_nnode

   subroutine test_ngroup()
      integer(int32) :: ngroups, mygroup
      logical :: test_ok

      if (me == 0) print *, "Test 11: ddi_ngroup..."

      call ddi_ngroup(ngroups, mygroup)

      test_ok = .true.
      ! Before group creation, should be 1 group with ID 0
      if (ngroups < 1) test_ok = .false.
      if (mygroup < 0 .or. mygroup >= ngroups) test_ok = .false.

      call ddi_sync(0_int64)

      if (me == 0) then
         if (test_ok) then
            print *, "  PASS: ddi_ngroup"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: ddi_ngroup"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_ngroup

   subroutine test_async_comm()
      real(dp) :: send_buf(10), recv_buf(10)
      integer(int32) :: send_req, recv_req, i
      logical :: test_ok

      if (me == 0) print *, "Test 12: ddi_isend/ddi_irecv/ddi_wait..."

      test_ok = .true.

      if (np >= 2) then
         if (me == 0) then
            ! Rank 0 sends to rank 1
            do i = 1, 10
               send_buf(i) = real(i*100, dp)
            end do
            call ddi_isend(send_buf, 80, 1, send_req)
            call ddi_wait(send_req)
         else if (me == 1) then
            ! Rank 1 receives from rank 0
            call ddi_irecv(recv_buf, 80, 0, recv_req)
            call ddi_wait(recv_req)

            ! Verify received data
            do i = 1, 10
               if (abs(recv_buf(i) - real(i*100, dp)) > 1.0e-10_dp) then
                  test_ok = .false.
                  exit
               end if
            end do
         end if
      end if

      call ddi_sync(0_int64)

      ! Gather test results
      if (me == 1) then
         send_buf(1) = 0.0_dp
         if (.not. test_ok) send_buf(1) = 1.0_dp
         call ddi_send(send_buf, 8, 0)
      else if (me == 0 .and. np >= 2) then
         call ddi_recv(recv_buf, 8, 1)
         if (recv_buf(1) > 0.5_dp) test_ok = .false.
      end if

      call ddi_sync(0_int64)

      if (me == 0) then
         if (test_ok) then
            print *, "  PASS: ddi_isend/ddi_irecv/ddi_wait"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: ddi_isend/ddi_irecv/ddi_wait"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_async_comm

   subroutine test_comm_ops()
      integer(int32) :: np_comm, me_comm
      real(dp) :: buffer_f(5), expected
      integer(int64) :: buffer_i(5), expected_i  ! int64 for GAMESS -i8 compatibility
      integer(int32) :: i
      logical :: test_ok

      if (me == 0) print *, "Test 13: ddi_*_comm functions..."

      test_ok = .true.

      ! Test ddi_nproc_comm on DDI_COMM_WORLD
      call ddi_nproc_comm(DDI_COMM_WORLD, np_comm, me_comm)
      if (np_comm /= np .or. me_comm /= me) test_ok = .false.

      ! Test ddi_sync_comm
      call ddi_sync_comm(DDI_COMM_WORLD)

      ! Test ddi_gsumf_comm
      do i = 1, 5
         buffer_f(i) = real(me + i, dp)
      end do
      call ddi_gsumf_comm(DDI_COMM_WORLD, buffer_f, 5)

      ! Verify: sum of (rank + i) = np*i + np*(np-1)/2
      do i = 1, 5
         expected = real(np*i + np*(np - 1)/2, dp)
         if (abs(buffer_f(i) - expected) > 1.0e-10_dp) then
            test_ok = .false.
            exit
         end if
      end do

      ! Test ddi_gsumi_comm
      do i = 1, 5
         buffer_i(i) = int(me + i, int64)
      end do
      call ddi_gsumi_comm(DDI_COMM_WORLD, buffer_i, 5)

      ! Verify: sum of (rank + i) = np*i + np*(np-1)/2
      do i = 1, 5
         expected_i = int(np*i + np*(np - 1)/2, int64)
         if (buffer_i(i) /= expected_i) then
            test_ok = .false.
            exit
         end if
      end do

      ! Test ddi_bcast_comm (double precision)
      buffer_f = 0.0_dp
      if (me == 0) then
         do i = 1, 5
            buffer_f(i) = real(i*10, dp)
         end do
      end if
      call ddi_bcast_comm(DDI_COMM_WORLD, 'F', buffer_f, 5, 0)

      ! Verify broadcast values
      do i = 1, 5
         if (abs(buffer_f(i) - real(i*10, dp)) > 1.0e-10_dp) then
            test_ok = .false.
            exit
         end if
      end do

      call ddi_sync(0_int64)

      if (me == 0) then
         if (test_ok) then
            print *, "  PASS: ddi_*_comm functions"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: ddi_*_comm functions"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_comm_ops

   subroutine test_comm_create_destroy()
      integer(int32) :: new_comm_id, np_new, me_new, color
      logical :: test_ok

      if (me == 0) print *, "Test 14: ddi_comm_create/ddi_comm_destroy..."

      test_ok = .true.

      ! Split into 2 groups: even ranks (color=0) and odd ranks (color=1)
      if (mod(me, 2) == 0) then
         color = 0
      else
         color = 1
      end if

      call ddi_comm_create(DDI_COMM_WORLD, color, new_comm_id)

      ! Verify we got a valid comm_id (should be >= 4)
      if (new_comm_id < 4) test_ok = .false.

      ! Check the new communicator size
      call ddi_nproc_comm(new_comm_id, np_new, me_new)

      ! Even ranks: should have (np+1)/2 ranks
      ! Odd ranks: should have np/2 ranks
      if (color == 0) then
         if (np_new /= (np + 1)/2) test_ok = .false.
      else
         if (np_new /= np/2) test_ok = .false.
      end if

      ! Test sync on the new communicator
      call ddi_sync_comm(new_comm_id)

      ! Destroy the communicator
      call ddi_comm_destroy(new_comm_id)

      call ddi_sync(0_int64)

      if (me == 0) then
         if (test_ok) then
            print *, "  PASS: ddi_comm_create/ddi_comm_destroy"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: ddi_comm_create/ddi_comm_destroy"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_comm_create_destroy

   subroutine test_arr_fill()
      integer(int32) :: handle, nrows, ncols, i
      real(dp), allocatable :: verify(:)
      real(dp) :: fill_value
      logical :: test_ok

      if (me == 0) print *, "Test 15: ddi_arr_fill..."

      nrows = 10
      ncols = 8
      fill_value = 42.0_dp
      call ddi_create(nrows, ncols, handle, 0.0_dp)

      ! Fill the array with a constant
      call ddi_arr_fill(handle, fill_value)

      call ddi_sync(0_int64)

      ! Rank 0 verifies all values
      test_ok = .true.
      if (me == 0) then
         allocate (verify(nrows*ncols))
         call ddi_get(handle, 1, nrows, 1, ncols, verify)

         do i = 1, nrows*ncols
            if (abs(verify(i) - fill_value) > 1.0e-10_dp) then
               test_ok = .false.
               exit
            end if
         end do
         deallocate (verify)
      end if

      call ddi_destroy(handle)
      call ddi_sync(0_int64)

      if (me == 0) then
         if (test_ok) then
            print *, "  PASS: ddi_arr_fill"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: ddi_arr_fill"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_arr_fill

end program test_ddi_compat
