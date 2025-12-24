# Missing MPI Functions for DDI Replacement

Based on GAMESS usage analysis and your existing `pic_mpi` wrappers.

---

## What You Already Have ✅

Your `pic_mpi` already covers:
- ✅ `barrier` → DDI_SYNC (113 calls)
- ✅ `bcast` → DDI_BCAST (656 calls)
- ✅ `send/recv` → DDI_SEND/RECV (22 calls)
- ✅ `isend/irecv/wait` → DDI_ISEND/IRECV/WAIT (116 calls)
- ✅ `allgather` → DDI_MASTERS_BCAST uses this
- ✅ Communicator splitting → for node-level operations

**Your wrappers cover ~900 DDI calls already!**

---

## What You Need to Add to `pic_mpi`

### Priority 1: Allreduce (CRITICAL - 1,466 calls!)

```fortran
! Add to pic_mpi interface block:
interface allreduce
   module procedure :: comm_allreduce_real_dp
   module procedure :: comm_allreduce_real_dp_array
   module procedure :: comm_allreduce_integer
   module procedure :: comm_allreduce_integer_array
end interface allreduce

public :: allreduce

! Implementation:
subroutine comm_allreduce_real_dp_array(comm, sendbuf, recvbuf, op)
   use mpi_f08, only: MPI_SUM, MPI_Allreduce, MPI_IN_PLACE
   type(comm_t), intent(in) :: comm
   real(dp), intent(inout) :: sendbuf(:)  ! in-place operation
   real(dp), intent(out) :: recvbuf(:)
   integer, intent(in), optional :: op
   integer(int32) :: ierr
   integer :: mpi_op

   if (present(op)) then
      mpi_op = op
   else
      mpi_op = MPI_SUM
   end if

   call MPI_Allreduce(MPI_IN_PLACE, sendbuf, size(sendbuf), &
                      MPI_DOUBLE_PRECISION, mpi_op, comm%m_comm, ierr)
end subroutine comm_allreduce_real_dp_array

! Similar for integer, integer_array, real_dp scalar
```

**Covers**: DDI_GSUMF (1,301 calls) + DDI_GSUMI (165 calls)

---

### Priority 2: MPI-3 RMA Operations (560 calls)

Add to `pic_mpi`:

```fortran
use mpi_f08, only: MPI_Win, MPI_Win_create, MPI_Win_free, &
                   MPI_Get, MPI_Put, MPI_Accumulate, &
                   MPI_Win_fence, MPI_Win_lock, MPI_Win_unlock

! Window type
type :: win_t
   private
   type(MPI_Win) :: m_win = MPI_WIN_NULL
   logical :: is_valid = .false.
contains
   procedure :: fence => win_fence
   procedure :: lock => win_lock
   procedure :: unlock => win_unlock
   procedure :: finalize => win_finalize
end type win_t

public :: win_t, win_create

interface win_create
   module procedure create_win_real_dp
end interface

function create_win_real_dp(comm, base, size) result(win)
   type(comm_t), intent(in) :: comm
   real(dp), target, intent(in) :: base(:)
   integer(int64), intent(in) :: size
   type(win_t) :: win
   integer :: ierr

   call MPI_Win_create(base, size, 8_int64, MPI_INFO_NULL, &
                       comm%get(), win%m_win, ierr)
   win%is_valid = .true.
end function create_win_real_dp

! Get operation
subroutine win_get(win, target_rank, target_offset, count, buffer)
   class(win_t), intent(in) :: win
   integer, intent(in) :: target_rank, target_offset, count
   real(dp), intent(out) :: buffer(:)
   integer :: ierr

   call MPI_Get(buffer, count, MPI_DOUBLE_PRECISION, &
                target_rank, target_offset, count, &
                MPI_DOUBLE_PRECISION, win%m_win, ierr)
end subroutine win_get

! Put operation
subroutine win_put(win, target_rank, target_offset, count, buffer)
   class(win_t), intent(in) :: win
   integer, intent(in) :: target_rank, target_offset, count
   real(dp), intent(in) :: buffer(:)
   integer :: ierr

   call MPI_Put(buffer, count, MPI_DOUBLE_PRECISION, &
                target_rank, target_offset, count, &
                MPI_DOUBLE_PRECISION, win%m_win, ierr)
end subroutine win_put

! Accumulate operation (CRITICAL for DDI_ACC)
subroutine win_accumulate(win, target_rank, target_offset, count, buffer)
   use mpi_f08, only: MPI_SUM
   class(win_t), intent(in) :: win
   integer, intent(in) :: target_rank, target_offset, count
   real(dp), intent(in) :: buffer(:)
   integer :: ierr

   call MPI_Accumulate(buffer, count, MPI_DOUBLE_PRECISION, &
                       target_rank, target_offset, count, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, win%m_win, ierr)
end subroutine win_accumulate

! Fence for synchronization
subroutine win_fence(this)
   class(win_t), intent(in) :: this
   integer :: ierr
   call MPI_Win_fence(0, this%m_win, ierr)
end subroutine win_fence

! Lock/unlock for passive target RMA
subroutine win_lock(this, rank, lock_type)
   use mpi_f08, only: MPI_LOCK_SHARED, MPI_LOCK_EXCLUSIVE
   class(win_t), intent(in) :: this
   integer, intent(in) :: rank
   integer, intent(in), optional :: lock_type
   integer :: ierr, ltype

   if (present(lock_type)) then
      ltype = lock_type
   else
      ltype = MPI_LOCK_SHARED
   end if

   call MPI_Win_lock(ltype, rank, 0, this%m_win, ierr)
end subroutine win_lock

subroutine win_unlock(this, rank)
   class(win_t), intent(in) :: this
   integer, intent(in) :: rank
   integer :: ierr
   call MPI_Win_unlock(rank, this%m_win, ierr)
end subroutine win_unlock
```

**Covers**: DDI_GET (435 calls) + DDI_PUT (100 calls) + DDI_ACC (25 calls)

---

### Priority 3: Fetch-and-Add (352 calls)

```fortran
use mpi_f08, only: MPI_Fetch_and_op, MPI_SUM

subroutine win_fetch_and_add(win, target_rank, offset, value, result)
   class(win_t), intent(in) :: win
   integer, intent(in) :: target_rank
   integer(int64), intent(in) :: offset
   integer(int64), intent(in) :: value
   integer(int64), intent(out) :: result
   integer :: ierr

   call MPI_Fetch_and_op(value, result, MPI_INTEGER8, &
                         target_rank, offset, MPI_SUM, win%m_win, ierr)
end subroutine win_fetch_and_add
```

**Covers**: DDI_DLBNEXT (352 calls) + DDI_DLBRESET (292 calls)

---

## Summary: Lines of Code to Add

| Feature | Lines | Difficulty | GAMESS Calls Covered |
|---------|-------|-----------|---------------------|
| **Allreduce** | ~80 | Easy | 1,466 (31%) |
| **RMA (Get/Put/Acc)** | ~200 | Medium | 560 (12%) |
| **Fetch-and-add** | ~30 | Easy | 352 (8%) |
| **Total** | **~310 lines** | | **2,378 calls (51%)** |

Combined with what you already have, this gives you **~3,300 calls (70% of GAMESS)**.

---

## What You DON'T Need from DDI

Since you have modern MPI wrappers:
- ❌ Data server model (DDI's main complexity)
- ❌ TCP/IP sockets (you have MPI)
- ❌ LAPI support (obsolete)
- ❌ Custom collective tree algorithms (MPI's are better)
- ❌ Most of DDI's 15,000 lines

---

## DDI Functions You Need to Write (Beyond MPI Wrappers)

These are DDI-specific logic, not pure MPI wrappers:

### 1. Distributed Array Management (~200 lines)

```fortran
module ddi_arrays
   use pic_mpi
   use pic_types

   type :: ddi_array_t
      integer :: nrows, ncols
      integer, allocatable :: col_map(:)      ! Which rank owns each column
      integer :: first_col, ncols_local
      real(dp), allocatable :: local_data(:)
      type(win_t) :: win
      logical :: active = .false.
   end type ddi_array_t

contains

   subroutine ddi_create(comm, nrows, ncols, handle)
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: nrows, ncols
      integer, intent(out) :: handle

      ! Column-wise distribution
      ! Allocate local portion
      ! Create MPI window
      ! Store in global array registry
   end subroutine ddi_create

   subroutine ddi_get(handle, ilo, ihi, jlo, jhi, buffer)
      ! Use win%get() for remote columns
      ! Direct memcpy for local columns
   end subroutine ddi_get

   subroutine ddi_put(handle, ilo, ihi, jlo, jhi, buffer)
      ! Use win%put() for remote columns
   end subroutine ddi_put

   subroutine ddi_acc(handle, ilo, ihi, jlo, jhi, buffer)
      ! Use win%accumulate() for remote columns
      ! Local atomic add for local columns
   end subroutine ddi_acc

end module ddi_arrays
```

### 2. Load Balancing (~50 lines)

```fortran
module ddi_dlb
   use pic_mpi

   type(win_t) :: dlb_window
   integer(int64) :: dlb_counter = 0

contains

   subroutine ddi_dlb_reset(comm)
      ! Reset counter to 0 on rank 0
      call dlb_window%fence()
   end subroutine ddi_dlb_reset

   subroutine ddi_dlb_next(comm, counter)
      integer(int64), intent(out) :: counter
      integer(int64) :: one = 1

      ! Atomic fetch-and-add
      call dlb_window%lock(0)
      call dlb_window%fetch_and_add(0, 0_int64, one, counter)
      call dlb_window%unlock(0)
   end subroutine ddi_dlb_next

end module ddi_dlb
```

### 3. Column Distribution Helper (~30 lines)

```fortran
function build_column_map(ncols, nranks) result(map)
   integer, intent(in) :: ncols, nranks
   integer, allocatable :: map(:)
   integer :: col, rank, cols_per_rank, extra, my_cols

   allocate(map(ncols))
   cols_per_rank = ncols / nranks
   extra = mod(ncols, nranks)

   col = 0
   do rank = 0, nranks-1
      my_cols = cols_per_rank
      if (rank < extra) my_cols = my_cols + 1

      map(col:col+my_cols-1) = rank
      col = col + my_cols
   end do
end function build_column_map
```

---

## Complete Implementation Estimate

| Component | Lines | Where |
|-----------|-------|-------|
| **Add to pic_mpi** | 310 | RMA + allreduce + fetch-and-add |
| **DDI array logic** | 200 | New module |
| **DDI load balancing** | 50 | New module |
| **Helpers** | 50 | Column maps, etc. |
| **Fortran interface** | 100 | C-to-Fortran if GAMESS is Fortran |
| **Total** | **710 lines** | |

Compare to DDI's 15,000 lines!

---

## Testing Strategy

1. **Test MPI additions** (310 lines):
```fortran
! Test allreduce
call comm%allreduce(data)

! Test RMA window
win = win_create(comm, array, size)
call win%fence()
call win%get(target_rank, offset, count, buffer)
call win%fence()
```

2. **Test DDI arrays** (200 lines):
```fortran
call ddi_create(comm, 100, 100, handle)
call ddi_put(handle, 0, 9, 0, 9, my_data)
call ddi_get(handle, 0, 9, 0, 9, buffer)
! Verify buffer == my_data
```

3. **Integration with GAMESS**:
```fortran
! Replace DDI_GSUMF calls
call ddi_gsumf(array, n) → call comm%allreduce(array)

! Replace DDI_GET calls
call ddi_get(handle, ...) → your ddi_get wrapper
```

---

## Recommendation

1. **Add 310 lines to `pic_mpi`** (allreduce, RMA, fetch-and-add) - **1 week**
2. **Write 250 lines of DDI logic** (arrays, load balancing) - **1 week**
3. **Test with GAMESS** - **1 week**

**Total: 3 weeks, 560 lines of code**

vs. maintaining DDI's 15,000 lines of legacy cruft.

You're 80% there already!
