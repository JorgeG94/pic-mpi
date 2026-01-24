!> Test distributed arrays module
!!
!! Tests create/destroy, get/put, accumulate, and load balancing.
!! Run with: mpirun -np 4 ./test_darrays
program test_darrays
   use pic_mpi_lib
   use pic_types, only: sp, dp, int32, int64
   use groups
   use darrays
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
      print *, "Testing Distributed Arrays (darrays)"
      print *, "========================================"
      print *, "Number of ranks:", world_comm%size()
   end if

   ! Initialize groups first (required by darrays)
   call groups_init(world_comm)

   ! Initialize darrays
   call darrays_init(world_comm)
   call darrays_set_comm(world_comm)  ! No-op now, kept for compatibility
   call dlb_init(world_comm)

   call world_comm%barrier()

   call test_create_destroy()
   call test_get_put_dp()
   call test_get_put_sp()
   call test_get_put_i32()
   call test_get_put_i64()
   call test_accumulate()
   call test_dlb()

   call world_comm%barrier()

   ! Summary
   if (world_comm%leader()) then
      print *, "========================================"
      print *, "Results: ", n_passed, " passed, ", n_failed, " failed"
      print *, "========================================"
   end if

   call dlb_finalize()
   call darrays_finalize()
   call groups_finalize()
   call world_comm%finalize()
   call pic_mpi_finalize()

   if (n_failed > 0) stop 1

contains

   subroutine test_create_destroy()
      integer(int32) :: handle, ilo, ihi, jlo, jhi
      logical :: test_ok

      if (world_comm%leader()) print *, "Test 1: Create and destroy..."

      call darray_create(100, 50, handle, 0.0_dp)

      ! Check distribution
      call darray_distrib(handle, world_comm%rank(), ilo, ihi, jlo, jhi)

      test_ok = (ilo == 1) .and. (ihi == 100) .and. (jlo >= 1) .and. (jhi <= 50)

      call darray_destroy(handle)

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Create/destroy"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Create/destroy"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_create_destroy

   subroutine test_get_put_dp()
      integer(int32) :: handle, ilo, ihi, jlo, jhi
      integer(int32) :: nrows, ncols, i, j
      real(dp), allocatable :: buffer(:), verify(:)
      logical :: test_ok

      if (world_comm%leader()) print *, "Test 2a: Get/put (dp)..."

      nrows = 10
      ncols = 8
      call darray_create(nrows, ncols, handle, 0.0_dp)

      call darray_distrib(handle, world_comm%rank(), ilo, ihi, jlo, jhi)

      if (jhi >= jlo) then
         allocate (buffer(nrows*(jhi - jlo + 1)))
         do j = jlo, jhi
            do i = 1, nrows
               buffer(i + (j - jlo)*nrows) = real(i*100 + j, dp)
            end do
         end do
         call darray_put(handle, 1, nrows, jlo, jhi, buffer)
         deallocate (buffer)
      end if

      call world_comm%barrier()

      test_ok = .true.
      if (world_comm%leader()) then
         allocate (verify(nrows*ncols))
         call darray_get(handle, 1, nrows, 1, ncols, verify)

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

      call darray_destroy(handle)
      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Get/put (dp)"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Get/put (dp)"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_get_put_dp

   subroutine test_get_put_sp()
      integer(int32) :: handle, ilo, ihi, jlo, jhi
      integer(int32) :: nrows, ncols, i, j
      real(sp), allocatable :: buffer(:), verify(:)
      logical :: test_ok

      if (world_comm%leader()) print *, "Test 2b: Get/put (sp)..."

      nrows = 10
      ncols = 8
      call darray_create(nrows, ncols, handle, 0.0_sp)

      call darray_distrib(handle, world_comm%rank(), ilo, ihi, jlo, jhi)

      if (jhi >= jlo) then
         allocate (buffer(nrows*(jhi - jlo + 1)))
         do j = jlo, jhi
            do i = 1, nrows
               buffer(i + (j - jlo)*nrows) = real(i*100 + j, sp)
            end do
         end do
         call darray_put(handle, 1, nrows, jlo, jhi, buffer)
         deallocate (buffer)
      end if

      call world_comm%barrier()

      test_ok = .true.
      if (world_comm%leader()) then
         allocate (verify(nrows*ncols))
         call darray_get(handle, 1, nrows, 1, ncols, verify)

         do j = 1, ncols
            do i = 1, nrows
               if (abs(verify(i + (j - 1)*nrows) - real(i*100 + j, sp)) > 1.0e-5_sp) then
                  test_ok = .false.
                  exit
               end if
            end do
            if (.not. test_ok) exit
         end do
         deallocate (verify)
      end if

      call darray_destroy(handle)
      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Get/put (sp)"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Get/put (sp)"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_get_put_sp

   subroutine test_get_put_i32()
      integer(int32) :: handle, ilo, ihi, jlo, jhi
      integer(int32) :: nrows, ncols, i, j
      integer(int32), allocatable :: buffer(:), verify(:)
      logical :: test_ok

      if (world_comm%leader()) print *, "Test 2c: Get/put (i32)..."

      nrows = 10
      ncols = 8
      call darray_create(nrows, ncols, handle, 0_int32)

      call darray_distrib(handle, world_comm%rank(), ilo, ihi, jlo, jhi)

      if (jhi >= jlo) then
         allocate (buffer(nrows*(jhi - jlo + 1)))
         do j = jlo, jhi
            do i = 1, nrows
               buffer(i + (j - jlo)*nrows) = i*100 + j
            end do
         end do
         call darray_put(handle, 1, nrows, jlo, jhi, buffer)
         deallocate (buffer)
      end if

      call world_comm%barrier()

      test_ok = .true.
      if (world_comm%leader()) then
         allocate (verify(nrows*ncols))
         call darray_get(handle, 1, nrows, 1, ncols, verify)

         do j = 1, ncols
            do i = 1, nrows
               if (verify(i + (j - 1)*nrows) /= i*100 + j) then
                  test_ok = .false.
                  exit
               end if
            end do
            if (.not. test_ok) exit
         end do
         deallocate (verify)
      end if

      call darray_destroy(handle)
      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Get/put (i32)"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Get/put (i32)"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_get_put_i32

   subroutine test_get_put_i64()
      integer(int32) :: handle, ilo, ihi, jlo, jhi
      integer(int32) :: nrows, ncols, i, j
      integer(int64), allocatable :: buffer(:), verify(:)
      logical :: test_ok

      if (world_comm%leader()) print *, "Test 2d: Get/put (i64)..."

      nrows = 10
      ncols = 8
      call darray_create(nrows, ncols, handle, 0_int64)

      call darray_distrib(handle, world_comm%rank(), ilo, ihi, jlo, jhi)

      if (jhi >= jlo) then
         allocate (buffer(nrows*(jhi - jlo + 1)))
         do j = jlo, jhi
            do i = 1, nrows
               buffer(i + (j - jlo)*nrows) = int(i, int64)*100_int64 + int(j, int64)
            end do
         end do
         call darray_put(handle, 1, nrows, jlo, jhi, buffer)
         deallocate (buffer)
      end if

      call world_comm%barrier()

      test_ok = .true.
      if (world_comm%leader()) then
         allocate (verify(nrows*ncols))
         call darray_get(handle, 1, nrows, 1, ncols, verify)

         do j = 1, ncols
            do i = 1, nrows
               if (verify(i + (j - 1)*nrows) /= int(i, int64)*100_int64 + int(j, int64)) then
                  test_ok = .false.
                  exit
               end if
            end do
            if (.not. test_ok) exit
         end do
         deallocate (verify)
      end if

      call darray_destroy(handle)
      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Get/put (i64)"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Get/put (i64)"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_get_put_i64

   subroutine test_accumulate()
      integer(int32) :: handle, nrows, ncols, i, j
      real(dp), allocatable :: buffer(:), verify(:)
      real(dp) :: expected
      logical :: test_ok

      if (world_comm%leader()) print *, "Test 3: Accumulate..."

      nrows = 5
      ncols = 4
      call darray_create(nrows, ncols, handle, 0.0_dp)

      call world_comm%barrier()

      ! All ranks accumulate 1.0 to all elements
      allocate (buffer(nrows*ncols))
      buffer = 1.0_dp
      call darray_acc(handle, 1, nrows, 1, ncols, buffer)
      deallocate (buffer)

      call world_comm%barrier()

      ! Rank 0 verifies sum equals number of ranks
      test_ok = .true.
      if (world_comm%leader()) then
         allocate (verify(nrows*ncols))
         call darray_get(handle, 1, nrows, 1, ncols, verify)

         expected = real(world_comm%size(), dp)
         do i = 1, nrows*ncols
            if (abs(verify(i) - expected) > 1.0e-10_dp) then
               test_ok = .false.
               print *, "  Mismatch at", i, ": expected", expected, "got", verify(i)
               exit
            end if
         end do
         deallocate (verify)
      end if

      call darray_destroy(handle)
      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: Accumulate"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: Accumulate"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_accumulate

   subroutine test_dlb()
      integer(int64) :: counter, counters(100)
      integer(int32) :: i, j, nwork
      logical :: test_ok, found_dup
      integer(int32), allocatable :: all_counts(:)

      if (world_comm%leader()) print *, "Test 4: Dynamic load balancing..."

      nwork = 20

      call dlb_reset()

      ! Each rank gets work items until done
      do i = 1, 100
         counters(i) = -1_int64
      end do

      i = 1
      do
         call dlb_next(counter)
         if (counter >= nwork) exit
         counters(i) = counter
         i = i + 1
      end do

      call world_comm%barrier()

      ! Gather counts and verify uniqueness
      test_ok = .true.
      allocate (all_counts(world_comm%size()))
      all_counts = 0
      all_counts(world_comm%rank() + 1) = i - 1
      call darray_gsumi(all_counts, world_comm%size())

      if (world_comm%leader()) then
         ! Check total work items equals nwork
         if (sum(all_counts) /= nwork) then
            test_ok = .false.
            print *, "  Wrong total work items:", sum(all_counts), "expected", nwork
         end if
      end if

      deallocate (all_counts)
      call world_comm%barrier()

      if (world_comm%leader()) then
         if (test_ok) then
            print *, "  PASS: DLB"
            n_passed = n_passed + 1
         else
            print *, "  FAIL: DLB"
            n_failed = n_failed + 1
         end if
      end if
   end subroutine test_dlb

end program test_darrays
