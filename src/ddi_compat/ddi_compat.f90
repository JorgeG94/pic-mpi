!> DDI-compatible API wrapper (backward compatibility module)
!!
!! This module provides the DDI API for use with 'use ddi_compat'.
!! For external linkage (GAMESS-style), the symbols in ddi_api.f90
!! are used instead.
!!
!! Usage:
!!   use ddi_compat
!!   call ddi_init()
!!   call ddi_create(idim, jdim, handle)
!!   call ddi_put(handle, ilo, ihi, jlo, jhi, buffer)
!!   call ddi_get(handle, ilo, ihi, jlo, jhi, buffer)
!!   call ddi_destroy(handle)
!!   call ddi_finalize()
!!
module ddi_compat
   use pic_types, only: int32, int64, sp, dp
   use ddi_impl, only: DDI_COMM_WORLD, DDI_COMM_GROUP, DDI_COMM_MASTER
   use ddi_impl, only: ddi_init_impl, ddi_finalize_impl, ddi_memory_impl, ddi_output_impl
   use ddi_impl, only: ddi_create_dp_impl, ddi_create_i32_impl, ddi_create_i64_impl
   use ddi_impl, only: ddi_destroy_impl, ddi_distrib_impl, ddi_zero_impl
   use ddi_impl, only: ddi_get_dp_impl, ddi_get_i32_impl, ddi_get_i64_impl
   use ddi_impl, only: ddi_put_dp_impl, ddi_put_i32_impl, ddi_put_i64_impl
   use ddi_impl, only: ddi_acc_dp_impl, ddi_acc_i32_impl, ddi_acc_i64_impl
   use ddi_impl, only: ddi_gsumf_impl, ddi_gsumi_impl, ddi_bcast_impl
   use ddi_impl, only: ddi_send_impl, ddi_recv_impl, ddi_recvany_impl
   use ddi_impl, only: ddi_isend_impl, ddi_irecv_impl, ddi_wait_impl
   use ddi_impl, only: ddi_sync_impl, ddi_nproc_impl, ddi_nnode_impl, ddi_ngroup_impl
   use ddi_impl, only: ddi_dlbreset_impl, ddi_dlbnext_impl
   use ddi_impl, only: ddi_gdlbreset_impl, ddi_gdlbnext_impl
   use ddi_impl, only: ddi_group_create_impl, ddi_scope_impl, ddi_ascope_impl
   use ddi_impl, only: ddi_timer_reset_impl, ddi_timer_output_impl
   use ddi_impl, only: ddi_comm_create_impl, ddi_comm_destroy_impl
   use ddi_impl, only: ddi_nproc_comm_impl, ddi_sync_comm_impl
   use ddi_impl, only: ddi_gsumf_comm_impl, ddi_gsumi_comm_impl, ddi_bcast_comm_impl
   use ddi_impl, only: ddi_arr_fill_impl, ddi_arr_scale_impl
   use ddi_impl, only: ddi_arr_min_impl, ddi_arr_max_impl
   use ddi_impl, only: ddi_arr_dot_impl, ddi_arr_add_impl, ddi_arr_acc_impl
   use ddi_impl, only: ddi_getacc_impl, ddi_group_create_custom_impl

   implicit none
   private

   ! Re-export DDI communicator constants
   public :: DDI_COMM_WORLD, DDI_COMM_GROUP, DDI_COMM_MASTER

   ! Initialization and finalization
   public :: ddi_init, ddi_finalize, ddi_memory, ddi_output

   ! Array operations (generic interfaces for dp, i32, i64)
   public :: ddi_create, ddi_destroy, ddi_distrib, ddi_zero
   public :: ddi_get, ddi_put, ddi_acc

   ! Communication
   public :: ddi_gsumf, ddi_gsumi, ddi_bcast
   public :: ddi_send, ddi_recv, ddi_recvany
   public :: ddi_isend, ddi_irecv, ddi_wait

   ! Synchronization
   public :: ddi_sync, ddi_nproc, ddi_nnode, ddi_ngroup

   ! Load balancing
   public :: ddi_dlbreset, ddi_dlbnext
   public :: ddi_gdlbreset, ddi_gdlbnext

   ! Groups and scope
   public :: ddi_group_create, ddi_scope, ddi_ascope

   ! Timer
   public :: ddi_timer_reset, ddi_timer_output

   ! Communicator-specific operations
   public :: ddi_comm_create, ddi_comm_destroy
   public :: ddi_nproc_comm, ddi_sync_comm
   public :: ddi_gsumf_comm, ddi_gsumi_comm, ddi_bcast_comm

   ! Array-wide operations
   public :: ddi_arr_fill, ddi_arr_scale
   public :: ddi_arr_min, ddi_arr_max
   public :: ddi_arr_dot, ddi_arr_add, ddi_arr_acc

   ! Combined operations
   public :: ddi_getacc

   ! Custom group creation
   public :: ddi_group_create_custom

   ! Generic interfaces for array operations
   interface ddi_create
      module procedure ddi_create_dp
      module procedure ddi_create_i32
      module procedure ddi_create_i64
   end interface ddi_create

   interface ddi_get
      module procedure ddi_get_dp
      module procedure ddi_get_i32
      module procedure ddi_get_i64
   end interface ddi_get

   interface ddi_put
      module procedure ddi_put_dp
      module procedure ddi_put_i32
      module procedure ddi_put_i64
   end interface ddi_put

   interface ddi_acc
      module procedure ddi_acc_dp
      module procedure ddi_acc_i32
      module procedure ddi_acc_i64
   end interface ddi_acc

contains

   ! ========================================================================
   ! Initialization / Finalization
   ! ========================================================================

   subroutine ddi_init()
      call ddi_init_impl()
   end subroutine ddi_init

   subroutine ddi_finalize()
      call ddi_finalize_impl()
   end subroutine ddi_finalize

   subroutine ddi_memory(memddi)
      integer(int64), intent(in) :: memddi
      call ddi_memory_impl(memddi)
   end subroutine ddi_memory

   subroutine ddi_output(level)
      integer(int32), intent(in) :: level
      call ddi_output_impl(level)
   end subroutine ddi_output

   ! ========================================================================
   ! Array Creation
   ! ========================================================================

   subroutine ddi_create_dp(idim, jdim, handle, init_val)
      integer(int32), intent(in) :: idim, jdim
      integer(int32), intent(out) :: handle
      real(dp), intent(in) :: init_val
      call ddi_create_dp_impl(idim, jdim, handle, init_val)
   end subroutine ddi_create_dp

   subroutine ddi_create_i32(idim, jdim, handle, init_val)
      integer(int32), intent(in) :: idim, jdim
      integer(int32), intent(out) :: handle
      integer(int32), intent(in) :: init_val
      call ddi_create_i32_impl(idim, jdim, handle, init_val)
   end subroutine ddi_create_i32

   subroutine ddi_create_i64(idim, jdim, handle, init_val)
      integer(int32), intent(in) :: idim, jdim
      integer(int32), intent(out) :: handle
      integer(int64), intent(in) :: init_val
      call ddi_create_i64_impl(idim, jdim, handle, init_val)
   end subroutine ddi_create_i64

   subroutine ddi_destroy(handle)
      integer(int32), intent(in) :: handle
      call ddi_destroy_impl(handle)
   end subroutine ddi_destroy

   subroutine ddi_distrib(handle, rank, ilo, ihi, jlo, jhi)
      integer(int32), intent(in) :: handle, rank
      integer(int32), intent(out) :: ilo, ihi, jlo, jhi
      call ddi_distrib_impl(handle, rank, ilo, ihi, jlo, jhi)
   end subroutine ddi_distrib

   subroutine ddi_zero(handle)
      integer(int32), intent(in) :: handle
      call ddi_zero_impl(handle)
   end subroutine ddi_zero

   ! ========================================================================
   ! Get Operations
   ! ========================================================================

   subroutine ddi_get_dp(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      real(dp), intent(out) :: buffer(*)
      call ddi_get_dp_impl(handle, ilo, ihi, jlo, jhi, buffer)
   end subroutine ddi_get_dp

   subroutine ddi_get_i32(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      integer(int32), intent(out) :: buffer(*)
      call ddi_get_i32_impl(handle, ilo, ihi, jlo, jhi, buffer)
   end subroutine ddi_get_i32

   subroutine ddi_get_i64(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      integer(int64), intent(out) :: buffer(*)
      call ddi_get_i64_impl(handle, ilo, ihi, jlo, jhi, buffer)
   end subroutine ddi_get_i64

   ! ========================================================================
   ! Put Operations
   ! ========================================================================

   subroutine ddi_put_dp(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      real(dp), intent(in) :: buffer(*)
      call ddi_put_dp_impl(handle, ilo, ihi, jlo, jhi, buffer)
   end subroutine ddi_put_dp

   subroutine ddi_put_i32(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      integer(int32), intent(in) :: buffer(*)
      call ddi_put_i32_impl(handle, ilo, ihi, jlo, jhi, buffer)
   end subroutine ddi_put_i32

   subroutine ddi_put_i64(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      integer(int64), intent(in) :: buffer(*)
      call ddi_put_i64_impl(handle, ilo, ihi, jlo, jhi, buffer)
   end subroutine ddi_put_i64

   ! ========================================================================
   ! Accumulate Operations
   ! ========================================================================

   subroutine ddi_acc_dp(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      real(dp), intent(in) :: buffer(*)
      call ddi_acc_dp_impl(handle, ilo, ihi, jlo, jhi, buffer)
   end subroutine ddi_acc_dp

   subroutine ddi_acc_i32(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      integer(int32), intent(in) :: buffer(*)
      call ddi_acc_i32_impl(handle, ilo, ihi, jlo, jhi, buffer)
   end subroutine ddi_acc_i32

   subroutine ddi_acc_i64(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      integer(int64), intent(in) :: buffer(*)
      call ddi_acc_i64_impl(handle, ilo, ihi, jlo, jhi, buffer)
   end subroutine ddi_acc_i64

   ! ========================================================================
   ! Communication
   ! ========================================================================

   subroutine ddi_gsumf(tag, buffer, n)
      integer(int32), intent(in) :: tag, n
      real(dp), intent(inout) :: buffer(*)
      call ddi_gsumf_impl(tag, buffer, n)
   end subroutine ddi_gsumf

   subroutine ddi_gsumi(tag, buffer, n)
      integer(int32), intent(in) :: tag, n
      integer(int64), intent(inout) :: buffer(*)  ! int64 for GAMESS -i8 compatibility
      call ddi_gsumi_impl(tag, buffer, n)
   end subroutine ddi_gsumi

   subroutine ddi_bcast(tag, dtype, buffer, n, root)
      integer(int32), intent(in) :: tag, n, root
      character(len=1), intent(in) :: dtype
      real(dp), intent(inout), target :: buffer(*)
      call ddi_bcast_impl(tag, dtype, buffer, n, root)
   end subroutine ddi_bcast

   subroutine ddi_send(buffer, size, to)
      real(dp), intent(in) :: buffer(*)
      integer(int32), intent(in) :: size, to
      call ddi_send_impl(buffer, size, to)
   end subroutine ddi_send

   subroutine ddi_recv(buffer, size, from)
      real(dp), intent(out) :: buffer(*)
      integer(int32), intent(in) :: size, from
      call ddi_recv_impl(buffer, size, from)
   end subroutine ddi_recv

   subroutine ddi_recvany(buffer, size, from)
      real(dp), intent(out) :: buffer(*)
      integer(int32), intent(in) :: size
      integer(int32), intent(out) :: from
      call ddi_recvany_impl(buffer, size, from)
   end subroutine ddi_recvany

   subroutine ddi_isend(buffer, size, to, req)
      real(dp), intent(in) :: buffer(*)
      integer(int32), intent(in) :: size, to
      integer(int32), intent(out) :: req
      call ddi_isend_impl(buffer, size, to, req)
   end subroutine ddi_isend

   subroutine ddi_irecv(buffer, size, from, req)
      real(dp), intent(out) :: buffer(*)
      integer(int32), intent(in) :: size, from
      integer(int32), intent(out) :: req
      call ddi_irecv_impl(buffer, size, from, req)
   end subroutine ddi_irecv

   subroutine ddi_wait(req)
      integer(int32), intent(in) :: req
      call ddi_wait_impl(req)
   end subroutine ddi_wait

   ! ========================================================================
   ! Synchronization
   ! ========================================================================

   subroutine ddi_sync(tag)
      integer(int32), intent(in) :: tag
      call ddi_sync_impl(tag)
   end subroutine ddi_sync

   subroutine ddi_nproc(np, me)
      integer(int32), intent(out) :: np, me
      call ddi_nproc_impl(np, me)
   end subroutine ddi_nproc

   subroutine ddi_nnode(nnodes, mynode)
      integer(int32), intent(out) :: nnodes, mynode
      call ddi_nnode_impl(nnodes, mynode)
   end subroutine ddi_nnode

   subroutine ddi_ngroup(ngroups, mygroup)
      integer(int32), intent(out) :: ngroups, mygroup
      call ddi_ngroup_impl(ngroups, mygroup)
   end subroutine ddi_ngroup

   ! ========================================================================
   ! Load Balancing
   ! ========================================================================

   subroutine ddi_dlbreset()
      call ddi_dlbreset_impl()
   end subroutine ddi_dlbreset

   subroutine ddi_dlbnext(counter)
      integer(int64), intent(out) :: counter
      call ddi_dlbnext_impl(counter)
   end subroutine ddi_dlbnext

   subroutine ddi_gdlbreset()
      call ddi_gdlbreset_impl()
   end subroutine ddi_gdlbreset

   subroutine ddi_gdlbnext(counter)
      integer(int64), intent(out) :: counter
      call ddi_gdlbnext_impl(counter)
   end subroutine ddi_gdlbnext

   ! ========================================================================
   ! Groups and Scope
   ! ========================================================================

   subroutine ddi_group_create(ngroups, world_id, group_id, master_id)
      integer(int32), intent(in) :: ngroups
      integer(int32), intent(out) :: world_id, group_id, master_id
      call ddi_group_create_impl(ngroups, world_id, group_id, master_id)
   end subroutine ddi_group_create

   subroutine ddi_scope(comm_id)
      integer(int32), intent(in) :: comm_id
      call ddi_scope_impl(comm_id)
   end subroutine ddi_scope

   subroutine ddi_ascope(comm_id)
      integer(int32), intent(in) :: comm_id
      call ddi_ascope_impl(comm_id)
   end subroutine ddi_ascope

   ! ========================================================================
   ! Timer
   ! ========================================================================

   subroutine ddi_timer_reset()
      call ddi_timer_reset_impl()
   end subroutine ddi_timer_reset

   subroutine ddi_timer_output()
      call ddi_timer_output_impl()
   end subroutine ddi_timer_output

   ! ========================================================================
   ! Communicator-specific operations
   ! ========================================================================

   subroutine ddi_comm_create(parent_id, color, new_comm_id)
      integer(int32), intent(in) :: parent_id, color
      integer(int32), intent(out) :: new_comm_id
      call ddi_comm_create_impl(parent_id, color, new_comm_id)
   end subroutine ddi_comm_create

   subroutine ddi_comm_destroy(comm_id)
      integer(int32), intent(in) :: comm_id
      call ddi_comm_destroy_impl(comm_id)
   end subroutine ddi_comm_destroy

   subroutine ddi_nproc_comm(comm_id, np, me)
      integer(int32), intent(in) :: comm_id
      integer(int32), intent(out) :: np, me
      call ddi_nproc_comm_impl(comm_id, np, me)
   end subroutine ddi_nproc_comm

   subroutine ddi_sync_comm(comm_id)
      integer(int32), intent(in) :: comm_id
      call ddi_sync_comm_impl(comm_id)
   end subroutine ddi_sync_comm

   subroutine ddi_gsumf_comm(comm_id, buffer, n)
      integer(int32), intent(in) :: comm_id, n
      real(dp), intent(inout) :: buffer(*)
      call ddi_gsumf_comm_impl(comm_id, buffer, n)
   end subroutine ddi_gsumf_comm

   subroutine ddi_gsumi_comm(comm_id, buffer, n)
      integer(int32), intent(in) :: comm_id, n
      integer(int64), intent(inout) :: buffer(*)  ! int64 for GAMESS -i8 compatibility
      call ddi_gsumi_comm_impl(comm_id, buffer, n)
   end subroutine ddi_gsumi_comm

   subroutine ddi_bcast_comm(comm_id, dtype, buffer, n, root)
      integer(int32), intent(in) :: comm_id, n, root
      character(len=1), intent(in) :: dtype
      real(dp), intent(inout), target :: buffer(*)
      call ddi_bcast_comm_impl(comm_id, dtype, buffer, n, root)
   end subroutine ddi_bcast_comm

   ! ========================================================================
   ! Array-wide operations
   ! ========================================================================

   subroutine ddi_arr_fill(handle, value)
      integer(int32), intent(in) :: handle
      real(dp), intent(in) :: value
      call ddi_arr_fill_impl(handle, value)
   end subroutine ddi_arr_fill

   subroutine ddi_arr_scale(handle, scale)
      integer(int32), intent(in) :: handle
      real(dp), intent(in) :: scale
      call ddi_arr_scale_impl(handle, scale)
   end subroutine ddi_arr_scale

   subroutine ddi_arr_min(handle, minval)
      integer(int32), intent(in) :: handle
      real(dp), intent(out) :: minval
      call ddi_arr_min_impl(handle, minval)
   end subroutine ddi_arr_min

   subroutine ddi_arr_max(handle, maxval)
      integer(int32), intent(in) :: handle
      real(dp), intent(out) :: maxval
      call ddi_arr_max_impl(handle, maxval)
   end subroutine ddi_arr_max

   subroutine ddi_arr_dot(handle_a, handle_b, dotprod)
      integer(int32), intent(in) :: handle_a, handle_b
      real(dp), intent(out) :: dotprod
      call ddi_arr_dot_impl(handle_a, handle_b, dotprod)
   end subroutine ddi_arr_dot

   subroutine ddi_arr_add(handle_a, handle_b, handle_c)
      integer(int32), intent(in) :: handle_a, handle_b, handle_c
      call ddi_arr_add_impl(handle_a, handle_b, handle_c)
   end subroutine ddi_arr_add

   subroutine ddi_arr_acc(handle_a, handle_b, alpha)
      integer(int32), intent(in) :: handle_a, handle_b
      real(dp), intent(in) :: alpha
      call ddi_arr_acc_impl(handle_a, handle_b, alpha)
   end subroutine ddi_arr_acc

   ! ========================================================================
   ! Combined operations
   ! ========================================================================

   subroutine ddi_getacc(handle, ilo, ihi, jlo, jhi, buffer)
      integer(int32), intent(in) :: handle, ilo, ihi, jlo, jhi
      real(dp), intent(inout) :: buffer(*)
      call ddi_getacc_impl(handle, ilo, ihi, jlo, jhi, buffer)
   end subroutine ddi_getacc

   ! ========================================================================
   ! Custom group creation
   ! ========================================================================

   subroutine ddi_group_create_custom(ngroups, group_sizes, world_id, group_id, master_id)
      integer(int32), intent(in) :: ngroups
      integer(int32), intent(in) :: group_sizes(*)
      integer(int32), intent(out) :: world_id, group_id, master_id
      call ddi_group_create_custom_impl(ngroups, group_sizes, world_id, group_id, master_id)
   end subroutine ddi_group_create_custom

end module ddi_compat
