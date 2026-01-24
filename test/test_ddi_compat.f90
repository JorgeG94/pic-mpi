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

   call ddi_sync(0)

   call test_nproc()
   call test_create_destroy()
   call test_get_put()
   call test_accumulate()
   call test_gsumf()
   call test_gsumi()
   call test_dlb()
   call test_groups()

   call ddi_sync(0)

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

      call ddi_sync(0)

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

      call ddi_sync(0)

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

      call ddi_sync(0)

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
      call ddi_sync(0)

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

      call ddi_sync(0)

      ! All ranks accumulate 1.0 to all elements
      allocate (buffer(nrows*ncols))
      buffer = 1.0_dp
      call ddi_acc(handle, 1, nrows, 1, ncols, buffer)
      deallocate (buffer)

      call ddi_sync(0)

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
      call ddi_sync(0)

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

      call ddi_sync(0)

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
      integer(int32) :: buffer(10), expected, i
      logical :: test_ok

      if (me == 0) print *, "Test 6: ddi_gsumi..."

      ! Each rank contributes its rank value
      do i = 1, 10
         buffer(i) = me + i
      end do

      call ddi_gsumi(0, buffer, 10)

      ! Expected: sum of (rank + i) for all ranks
      test_ok = .true.
      do i = 1, 10
         expected = np*i + np*(np - 1)/2
         if (buffer(i) /= expected) then
            test_ok = .false.
            exit
         end if
      end do

      call ddi_sync(0)

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
      integer(int32) :: i, nwork, total_work(1)
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

      call ddi_sync(0)

      ! Verify total work done
      total_work(1) = i - 1
      call ddi_gsumi(0, total_work, 1)

      test_ok = .true.
      if (me == 0) then
         if (total_work(1) /= nwork) then
            test_ok = .false.
         end if
      end if

      call ddi_sync(0)

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

      call ddi_sync(0)

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

end program test_ddi_compat
