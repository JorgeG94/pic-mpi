!> Test building an owned comm_t from a caller's Fortran communicator handle
!!
!! This is the path a host application takes when it hands pic-mpi a
!! communicator it already owns -- GAMESS/DDI passing its compute
!! communicator, say -- rather than letting pic-mpi duplicate
!! MPI_COMM_WORLD.
!!
!! The raw MPI calls here go through `mpi_f08`, and the integer handles are
!! taken from its MPI_VAL components. The older `mpi` module would read more
!! naturally, since its communicators already ARE those integers -- but it
!! cannot be used: a PIC_USE_VAPAA build links MPI::MPI_C only and vapaa
!! supplies no legacy `mpi` module, so every `mpi`-module symbol is
!! undefined at link time there. mpi_f08 exists in all three MPI
!! configurations, vapaa's included.
!!
!! What has to hold:
!!   - the result is a DIFFERENT communicator (congruent, not identical),
!!     so pic-mpi's messages cannot match the caller's
!!   - rank and size agree with the original
!!   - collectives work on it
!!   - it works for a sub-communicator, not just world
!!   - finalize() frees pic-mpi's copy and leaves the caller's alone
!!   - a null handle yields an invalid comm_t
!!
!! Works at any size from 2 ranks up; CI runs it on 2.
!!
!! Run with: mpirun -np 2 ./test_comm_from_handle
program test_comm_from_handle
   use pic_mpi_lib, only: comm_t, comm_world, comm_from_handle, allreduce, &
                          pic_mpi_init, pic_mpi_finalize
   use mpi_f08, only: MPI_Comm, mpi_world => MPI_COMM_WORLD, &
                      mpi_null => MPI_COMM_NULL, &
                      MPI_Comm_split, MPI_Comm_compare, MPI_Comm_free, &
                      MPI_Comm_rank, MPI_Comm_size, MPI_Barrier, &
                      MPI_IDENT, MPI_CONGRUENT
   implicit none

   type(comm_t) :: world_comm
   integer :: n_passed, n_failed

   n_passed = 0
   n_failed = 0

   call pic_mpi_init()
   world_comm = comm_world()

   ! Two ranks is enough for every check here, and is what CI runners can
   ! give. Nothing below assumes a particular size: the split test compares
   ! each half's reduction against that half's own membership, which
   ! distinguishes a correct split from a broken one at two ranks (1 vs 2)
   ! just as it does at four (2 vs 4).
   if (world_comm%size() < 2) then
      if (world_comm%leader()) then
         print *, "ERROR: This test requires at least 2 MPI ranks"
      end if
      call pic_mpi_finalize()
      stop 1
   end if

   if (world_comm%leader()) then
      print *, "========================================"
      print *, "Testing comm_from_handle"
      print *, "========================================"
      print *, "Number of ranks:", world_comm%size()
   end if

   call world_comm%barrier()

   call test_wraps_world()
   call test_wraps_split_halves()
   call test_finalize_spares_the_original()
   call test_null_handle()

   call world_comm%barrier()

   if (world_comm%leader()) then
      print *, "========================================"
      print *, "Results: ", n_passed, " passed, ", n_failed, " failed"
      print *, "========================================"
   end if

   call allreduce(world_comm, n_failed)

   call pic_mpi_finalize()

   if (n_failed > 0) stop 1

contains

   !> The comm_t's communicator as a plain integer handle.
   !!
   !! comm_t%get() returns whatever the active backend calls a
   !! communicator: type(MPI_Comm) under mpi_f08, a bare integer under the
   !! legacy `mpi` module. The tests below need the integer form, to hand to
   !! the `mpi`-module routines they check against.
   function handle_of(comm) result(fhandle)
      type(comm_t), intent(in) :: comm
      integer :: fhandle
#if !defined(USE_LEGACY)
      type(MPI_Comm) :: held
#endif

#if defined(USE_LEGACY)
      ! Under the legacy backend a communicator already is its handle.
      fhandle = comm%get()
#else
      ! A component cannot be referenced on a function result directly, so
      ! the result is held before MPI_VAL is read off it.
      held = comm%get()
      fhandle = held%MPI_VAL
#endif
   end function handle_of

   !> The same handle space seen as an mpi_f08 communicator.
   !!
   !! Needed to hand a pic-mpi communicator to MPI_Comm_compare, whose
   !! arguments are type(MPI_Comm) whichever backend is active.
   function as_mpi_comm(fhandle) result(c)
      integer, intent(in) :: fhandle
      type(MPI_Comm) :: c

      c%MPI_VAL = fhandle
   end function as_mpi_comm

   !> Report one test's outcome, agreed across all ranks.
   !!
   !! Every check here is per-rank, and a failure on any single rank is a
   !! failure of the test -- so the flag is reduced before the leader
   !! prints, rather than trusting rank 0's own view.
   subroutine report(name, ok_local)
      character(len=*), intent(in) :: name
      logical, intent(in) :: ok_local
      integer :: bad

      bad = 0
      if (.not. ok_local) bad = 1
      call allreduce(world_comm, bad)

      call world_comm%barrier()

      if (world_comm%leader()) then
         if (bad == 0) then
            print *, "  PASS: "//name
            n_passed = n_passed + 1
         else
            print *, "  FAIL: "//name
            n_failed = n_failed + 1
         end if
      end if
   end subroutine report

   !> Wrapping world's handle gives a congruent but distinct communicator.
   subroutine test_wraps_world()
      type(comm_t) :: wrapped
      integer :: ierr, cmp, total
      logical :: ok

      ok = .true.

      wrapped = comm_from_handle(mpi_world%MPI_VAL)

      if (wrapped%is_null()) then
         ok = .false.
         print *, "  Rank", world_comm%rank(), "wrapped world came back null"
      end if

      if (wrapped%rank() /= world_comm%rank()) then
         ok = .false.
         print *, "  Rank", world_comm%rank(), "rank mismatch:", wrapped%rank()
      end if

      if (wrapped%size() /= world_comm%size()) then
         ok = .false.
         print *, "  Rank", world_comm%rank(), "size mismatch:", wrapped%size()
      end if

      ! The point of duplicating: same group, different context. MPI_IDENT
      ! would mean the handle was merely wrapped, and pic-mpi's traffic
      ! could then match the caller's.
      call MPI_Comm_compare(mpi_world, as_mpi_comm(handle_of(wrapped)), cmp, ierr)
      if (cmp /= MPI_CONGRUENT) then
         ok = .false.
         if (cmp == MPI_IDENT) then
            print *, "  Rank", world_comm%rank(), "communicator was wrapped, not duplicated"
         else
            print *, "  Rank", world_comm%rank(), "unexpected compare result:", cmp
         end if
      end if

      total = 1
      call allreduce(wrapped, total)
      if (total /= world_comm%size()) then
         ok = .false.
         print *, "  Rank", world_comm%rank(), "allreduce on wrapped comm gave", total
      end if

      call wrapped%finalize()

      call report("wrapping MPI_COMM_WORLD's handle", ok)
   end subroutine test_wraps_world

   !> A sub-communicator wraps just as well as world, and stays separate.
   subroutine test_wraps_split_halves()
      type(comm_t) :: wrapped
      integer :: ierr, half, total
      type(MPI_Comm) :: half_comm
      integer :: my_rank, half_size
      logical :: ok

      ok = .true.

      ! Lower half and upper half, by world rank.
      half = 0
      if (world_comm%rank() >= world_comm%size()/2) half = 1

      call MPI_Comm_split(mpi_world, half, world_comm%rank(), half_comm, ierr)

      wrapped = comm_from_handle(half_comm%MPI_VAL)

      call MPI_Comm_rank(half_comm, my_rank, ierr)
      call MPI_Comm_size(half_comm, half_size, ierr)

      if (wrapped%rank() /= my_rank .or. wrapped%size() /= half_size) then
         ok = .false.
         print *, "  Rank", world_comm%rank(), "half rank/size mismatch:", &
            wrapped%rank(), my_rank, wrapped%size(), half_size
      end if

      ! Each half sums only its own members. Had the wrap leaked across
      ! halves -- or silently landed on world -- this would come out as the
      ! world size instead. The two differ at every supported size: 1 vs 2
      ! on two ranks, 2 vs 4 on four.
      total = 1
      call allreduce(wrapped, total)
      if (total /= half_size) then
         ok = .false.
         print *, "  Rank", world_comm%rank(), "half allreduce gave", total, "expected", half_size
      end if

      ! The check above is only worth something if the half really is
      ! smaller than world; assert that rather than assume it.
      if (half_size >= world_comm%size()) then
         ok = .false.
         print *, "  Rank", world_comm%rank(), "half is not smaller than world:", &
            half_size, world_comm%size()
      end if

      call wrapped%finalize()
      call MPI_Comm_free(half_comm, ierr)

      call report("wrapping a split sub-communicator", ok)
   end subroutine test_wraps_split_halves

   !> finalize() must free only pic-mpi's duplicate.
   subroutine test_finalize_spares_the_original()
      type(comm_t) :: wrapped
      integer :: ierr, cmp
      type(MPI_Comm) :: own_comm
      logical :: ok

      ok = .true.

      ! A communicator this test owns, so that freeing it afterwards is
      ! legitimate -- freeing MPI_COMM_WORLD would not be.
      call MPI_Comm_split(mpi_world, 0, world_comm%rank(), own_comm, ierr)

      wrapped = comm_from_handle(own_comm%MPI_VAL)
      call wrapped%finalize()

      ! The original must still be usable. A barrier on a freed communicator
      ! is an error, so completing one is the check.
      ierr = 0
      call MPI_Barrier(own_comm, ierr)
      if (ierr /= 0) then
         ok = .false.
         print *, "  Rank", world_comm%rank(), "barrier on the original failed:", ierr
      end if

      ! And it still compares as a real communicator against itself.
      call MPI_Comm_compare(own_comm, own_comm, cmp, ierr)
      if (cmp /= MPI_IDENT) then
         ok = .false.
         print *, "  Rank", world_comm%rank(), "original no longer identical to itself:", cmp
      end if

      ! Its owner can still free it, which is the other half of "untouched".
      call MPI_Comm_free(own_comm, ierr)
      if (ierr /= 0) then
         ok = .false.
         print *, "  Rank", world_comm%rank(), "freeing the original failed:", ierr
      end if

      call report("finalize leaves the caller's communicator alone", ok)
   end subroutine test_finalize_spares_the_original

   !> A null handle is not an error, just an invalid comm_t.
   subroutine test_null_handle()
      type(comm_t) :: wrapped
      logical :: ok

      ok = .true.

      wrapped = comm_from_handle(mpi_null%MPI_VAL)

      if (.not. wrapped%is_null()) then
         ok = .false.
         print *, "  Rank", world_comm%rank(), "null handle produced a usable comm_t"
      end if

      call report("a null handle gives an invalid comm_t", ok)
   end subroutine test_null_handle

end program test_comm_from_handle
