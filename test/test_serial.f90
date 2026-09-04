!> The single-rank backend, exercised rather than merely compiled.
!!
!! Everything here is a claim `pic_mpi_serial` makes about what one rank
!! means, and each is checked. A backend that builds but reports `size() == 0`
!! or leaves a broadcast buffer untouched-but-wrong would pass a compile and
!! fail a program.
!!
!! Deliberately NOT run under `mpirun`: the point of this build is that there
!! is no MPI to launch it with.
program test_serial
   use pic_types, only: int32, dp
   use pic_mpi_lib, only: comm_t, comm_world, pic_mpi_init, pic_mpi_finalize, &
                          bcast, allgather, allreduce, iprobe, MPI_Status, &
                          MPI_ANY_SOURCE, MPI_ANY_TAG
   implicit none

   type(comm_t) :: comm
   type(MPI_Status) :: status
   integer(int32) :: ivalue, mine, gathered(1)
   real(dp) :: dvalue, dsend(3), drecv(3)
   logical :: pending
   integer :: failures

   failures = 0

   call pic_mpi_init()
   comm = comm_world()

   ! --- what one rank is ---------------------------------------------------
   call expect(comm%rank() == 0, "the only rank is rank 0")
   call expect(comm%size() == 1, "the only communicator has size 1")
   call expect(comm%leader(), "the only rank is the leader")
   call expect(.not. comm%is_null(), "comm_world is not null")

   ! --- a broadcast has already arrived ------------------------------------
   !
   ! There is no root to send from and no peer to send to, so the buffer must
   ! come back holding exactly what it went in holding. Checking that catches
   ! a body that zeroed it "to be safe".
   ivalue = 42_int32
   call bcast(comm, ivalue, 1_int32, 0_int32)
   call expect(ivalue == 42_int32, "a broadcast leaves the buffer alone")

   ! --- a gather returns this rank's own contribution -----------------------
   !
   ! `allgather` takes a scalar per rank and fills one slot per rank, so on
   ! one rank the result is a single element holding this rank's value.
   mine = 7_int32
   gathered = 0_int32
   call allgather(comm, mine, gathered)
   call expect(gathered(1) == 7_int32, &
               "a gather over one rank returns that rank's own value")

   ! --- a reduction over one rank is that rank's own contribution ----------
   !
   ! Unlike a send, a reduction IS defined at one rank -- for every op the
   ! answer is the single contribution -- so these are identities, not errors.
   ! In place there is nothing to do; the `_to` forms copy.
   dvalue = 3.5_dp
   call allreduce(comm, dvalue)
   call expect(dvalue == 3.5_dp, &
               "an in-place reduction leaves the buffer alone")

   dsend = [1.0_dp, 2.0_dp, 3.0_dp]
   drecv = -1.0_dp
   call allreduce(comm, dsend, drecv)
   call expect(all(drecv == dsend), &
               "a reduction into a buffer copies every element")

   ! `count` has to bound the same range MPI_Allreduce would be given, so the
   ! tail past it keeps whatever the caller left there.
   drecv = -1.0_dp
   call allreduce(comm, dsend, drecv, 2)
   call expect(drecv(1) == 1.0_dp .and. drecv(2) == 2.0_dp &
               .and. drecv(3) == -1.0_dp, &
               "`count` bounds the reduced range, as it does under MPI")

   ! --- nothing is ever pending --------------------------------------------
   !
   ! This one matters more than it looks: a manager loop polls `iprobe` until
   ! something arrives, so a backend that reported `.true.` here would spin
   ! forever rather than fail.
   pending = .true.
   call iprobe(comm, MPI_ANY_SOURCE, MPI_ANY_TAG, pending, status)
   call expect(.not. pending, "no message is ever pending on one rank")

   call pic_mpi_finalize()

   if (failures == 0) then
      print '(a)', "test_serial: PASS"
   else
      print '(a,i0,a)', "test_serial: FAIL (", failures, " checks)"
      error stop 1
   end if

contains

   subroutine expect(condition, what)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: what
      if (condition) then
         print '(a,a)', "  ok    ", what
      else
         print '(a,a)', "  FAIL  ", what
         failures = failures + 1
      end if
   end subroutine expect

end program test_serial
