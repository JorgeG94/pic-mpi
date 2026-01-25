!> External DDI API - linker-visible symbols for GAMESS
!!
!! This file provides external subroutines with bind(c, name="ddi_xxx_")
!! to match the Fortran name mangling expected by GAMESS. This allows
!! GAMESS code to call DDI functions via external linkage without
!! requiring 'use' statements.
!!
!! Note: GAMESS compiles with -i8 (8-byte integers), so all integer
!! parameters use int64 and are converted to int32 for internal use.
!!
!! Architecture:
!!   GAMESS Code (external calls, int64 integers)
!!        |
!!        v
!!   ddi_api.f90 (this file - bind(c) wrappers, int64 -> int32 conversion)
!!        |
!!        v
!!   ddi_impl.f90 (module implementation, int32 internally)
!!        |
!!        v
!!   darrays/groups (MPI-3 RMA operations)
!!

! ============================================================================
! Initialization / Finalization
! ============================================================================

subroutine ddi_init() bind(c, name="ddi_init_")
   use ddi_impl, only: ddi_init_impl
   call ddi_init_impl()
end subroutine ddi_init

subroutine ddi_finalize() bind(c, name="ddi_finalize_")
   use ddi_impl, only: ddi_finalize_impl
   call ddi_finalize_impl()
end subroutine ddi_finalize

subroutine ddi_memory(memddi) bind(c, name="ddi_memory_")
   use iso_fortran_env, only: int64
   use ddi_impl, only: ddi_memory_impl
   integer(int64), intent(in) :: memddi
   call ddi_memory_impl(memddi)
end subroutine ddi_memory

subroutine ddi_output(level) bind(c, name="ddi_output_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_output_impl
   integer(int64), intent(in) :: level
   call ddi_output_impl(int(level, int32))
end subroutine ddi_output

! ============================================================================
! Array Creation / Destruction
! ============================================================================

subroutine ddi_create(idim, jdim, handle) bind(c, name="ddi_create_")
   use iso_fortran_env, only: int64, int32, real64
   use ddi_impl, only: ddi_create_dp_impl
   integer(int64), intent(in) :: idim, jdim
   integer(int64), intent(out) :: handle
   integer(int32) :: handle32
   call ddi_create_dp_impl(int(idim, int32), int(jdim, int32), handle32, 0.0_real64)
   handle = int(handle32, int64)
end subroutine ddi_create

subroutine ddi_destroy(handle) bind(c, name="ddi_destroy_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_destroy_impl
   integer(int64), intent(in) :: handle
   call ddi_destroy_impl(int(handle, int32))
end subroutine ddi_destroy

subroutine ddi_distrib(handle, rank, ilo, ihi, jlo, jhi) bind(c, name="ddi_distrib_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_distrib_impl
   integer(int64), intent(in) :: handle, rank
   integer(int64), intent(out) :: ilo, ihi, jlo, jhi
   integer(int32) :: ilo32, ihi32, jlo32, jhi32
   call ddi_distrib_impl(int(handle, int32), int(rank, int32), ilo32, ihi32, jlo32, jhi32)
   ilo = int(ilo32, int64)
   ihi = int(ihi32, int64)
   jlo = int(jlo32, int64)
   jhi = int(jhi32, int64)
end subroutine ddi_distrib

subroutine ddi_zero(handle) bind(c, name="ddi_zero_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_zero_impl
   integer(int64), intent(in) :: handle
   call ddi_zero_impl(int(handle, int32))
end subroutine ddi_zero

! ============================================================================
! Data Transfer Operations
! ============================================================================

subroutine ddi_get(handle, ilo, ihi, jlo, jhi, buffer) bind(c, name="ddi_get_")
   use iso_fortran_env, only: int64, int32, real64
   use ddi_impl, only: ddi_get_dp_impl
   integer(int64), intent(in) :: handle, ilo, ihi, jlo, jhi
   real(real64), intent(out) :: buffer(*)
   call ddi_get_dp_impl(int(handle, int32), int(ilo, int32), int(ihi, int32), &
                        int(jlo, int32), int(jhi, int32), buffer)
end subroutine ddi_get

subroutine ddi_put(handle, ilo, ihi, jlo, jhi, buffer) bind(c, name="ddi_put_")
   use iso_fortran_env, only: int64, int32, real64
   use ddi_impl, only: ddi_put_dp_impl
   integer(int64), intent(in) :: handle, ilo, ihi, jlo, jhi
   real(real64), intent(in) :: buffer(*)
   call ddi_put_dp_impl(int(handle, int32), int(ilo, int32), int(ihi, int32), &
                        int(jlo, int32), int(jhi, int32), buffer)
end subroutine ddi_put

subroutine ddi_acc(handle, ilo, ihi, jlo, jhi, buffer) bind(c, name="ddi_acc_")
   use iso_fortran_env, only: int64, int32, real64
   use ddi_impl, only: ddi_acc_dp_impl
   integer(int64), intent(in) :: handle, ilo, ihi, jlo, jhi
   real(real64), intent(in) :: buffer(*)
   call ddi_acc_dp_impl(int(handle, int32), int(ilo, int32), int(ihi, int32), &
                        int(jlo, int32), int(jhi, int32), buffer)
end subroutine ddi_acc

! ============================================================================
! Global Operations
! ============================================================================

subroutine ddi_gsumf(tag, buffer, n) bind(c, name="ddi_gsumf_")
   use iso_fortran_env, only: int64, int32, real64
   use ddi_impl, only: ddi_gsumf_impl
   integer(int64), intent(in) :: tag, n
   real(real64), intent(inout) :: buffer(*)
   call ddi_gsumf_impl(int(tag, int32), buffer, int(n, int32))
end subroutine ddi_gsumf

subroutine ddi_gsumi(tag, buffer, n) bind(c, name="ddi_gsumi_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_gsumi_impl
   integer(int64), intent(in) :: tag, n
   integer(int64), intent(inout) :: buffer(*)
   ! Pass int64 buffer directly - impl uses MPI_INTEGER8
   call ddi_gsumi_impl(int(tag, int32), buffer, int(n, int32))
end subroutine ddi_gsumi

subroutine ddi_bcast(tag, dtype, buffer, n, root) bind(c, name="ddi_bcast_")
   use iso_fortran_env, only: int64, int32, real64
   use ddi_impl, only: ddi_bcast_impl
   integer(int64), intent(in) :: tag, n, root
   character(len=1), intent(in) :: dtype
   real(real64), intent(inout) :: buffer(*)
   call ddi_bcast_impl(int(tag, int32), dtype, buffer, int(n, int32), int(root, int32))
end subroutine ddi_bcast

! ============================================================================
! Point-to-Point Communication
! ============================================================================

subroutine ddi_send(buffer, size, to) bind(c, name="ddi_send_")
   use iso_fortran_env, only: int64, int32, real64
   use ddi_impl, only: ddi_send_impl
   real(real64), intent(in) :: buffer(*)
   integer(int64), intent(in) :: size, to
   call ddi_send_impl(buffer, int(size, int32), int(to, int32))
end subroutine ddi_send

subroutine ddi_recv(buffer, size, from) bind(c, name="ddi_recv_")
   use iso_fortran_env, only: int64, int32, real64
   use ddi_impl, only: ddi_recv_impl
   real(real64), intent(out) :: buffer(*)
   integer(int64), intent(in) :: size, from
   call ddi_recv_impl(buffer, int(size, int32), int(from, int32))
end subroutine ddi_recv

subroutine ddi_recvany(buffer, size, from) bind(c, name="ddi_recvany_")
   use iso_fortran_env, only: int64, int32, real64
   use ddi_impl, only: ddi_recvany_impl
   real(real64), intent(out) :: buffer(*)
   integer(int64), intent(in) :: size
   integer(int64), intent(out) :: from
   integer(int32) :: from32
   call ddi_recvany_impl(buffer, int(size, int32), from32)
   from = int(from32, int64)
end subroutine ddi_recvany

subroutine ddi_isend(buffer, size, to, req) bind(c, name="ddi_isend_")
   use iso_fortran_env, only: int64, int32, real64
   use ddi_impl, only: ddi_isend_impl
   real(real64), intent(in) :: buffer(*)
   integer(int64), intent(in) :: size, to
   integer(int64), intent(out) :: req
   integer(int32) :: req32
   call ddi_isend_impl(buffer, int(size, int32), int(to, int32), req32)
   req = int(req32, int64)
end subroutine ddi_isend

subroutine ddi_irecv(buffer, size, from, req) bind(c, name="ddi_irecv_")
   use iso_fortran_env, only: int64, int32, real64
   use ddi_impl, only: ddi_irecv_impl
   real(real64), intent(out) :: buffer(*)
   integer(int64), intent(in) :: size, from
   integer(int64), intent(out) :: req
   integer(int32) :: req32
   call ddi_irecv_impl(buffer, int(size, int32), int(from, int32), req32)
   req = int(req32, int64)
end subroutine ddi_irecv

subroutine ddi_wait(req) bind(c, name="ddi_wait_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_wait_impl
   integer(int64), intent(in) :: req
   call ddi_wait_impl(int(req, int32))
end subroutine ddi_wait

! ============================================================================
! Synchronization and Queries
! ============================================================================

subroutine ddi_sync(tag) bind(c, name="ddi_sync_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_sync_impl
   integer(int64), intent(in) :: tag
   call ddi_sync_impl(int(tag, int32))
end subroutine ddi_sync

subroutine ddi_nproc(np, me) bind(c, name="ddi_nproc_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_nproc_impl
   integer(int64), intent(out) :: np, me
   integer(int32) :: np32, me32
   call ddi_nproc_impl(np32, me32)
   np = int(np32, int64)
   me = int(me32, int64)
end subroutine ddi_nproc

subroutine ddi_nnode(nnodes, mynode) bind(c, name="ddi_nnode_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_nnode_impl
   integer(int64), intent(out) :: nnodes, mynode
   integer(int32) :: nnodes32, mynode32
   call ddi_nnode_impl(nnodes32, mynode32)
   nnodes = int(nnodes32, int64)
   mynode = int(mynode32, int64)
end subroutine ddi_nnode

subroutine ddi_ngroup(ngroups, mygroup) bind(c, name="ddi_ngroup_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_ngroup_impl
   integer(int64), intent(out) :: ngroups, mygroup
   integer(int32) :: ngroups32, mygroup32
   call ddi_ngroup_impl(ngroups32, mygroup32)
   ngroups = int(ngroups32, int64)
   mygroup = int(mygroup32, int64)
end subroutine ddi_ngroup

! ============================================================================
! Dynamic Load Balancing
! ============================================================================

subroutine ddi_dlbreset() bind(c, name="ddi_dlbreset_")
   use ddi_impl, only: ddi_dlbreset_impl
   call ddi_dlbreset_impl()
end subroutine ddi_dlbreset

subroutine ddi_dlbnext(counter) bind(c, name="ddi_dlbnext_")
   use iso_fortran_env, only: int64
   use ddi_impl, only: ddi_dlbnext_impl
   integer(int64), intent(out) :: counter
   call ddi_dlbnext_impl(counter)
end subroutine ddi_dlbnext

subroutine ddi_gdlbreset() bind(c, name="ddi_gdlbreset_")
   use ddi_impl, only: ddi_gdlbreset_impl
   call ddi_gdlbreset_impl()
end subroutine ddi_gdlbreset

subroutine ddi_gdlbnext(counter) bind(c, name="ddi_gdlbnext_")
   use iso_fortran_env, only: int64
   use ddi_impl, only: ddi_gdlbnext_impl
   integer(int64), intent(out) :: counter
   call ddi_gdlbnext_impl(counter)
end subroutine ddi_gdlbnext

! ============================================================================
! Groups and Scope
! ============================================================================

subroutine ddi_group_create(ngroups, world_id, group_id, master_id) bind(c, name="ddi_group_create_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_group_create_impl
   integer(int64), intent(in) :: ngroups
   integer(int64), intent(out) :: world_id, group_id, master_id
   integer(int32) :: world_id32, group_id32, master_id32
   call ddi_group_create_impl(int(ngroups, int32), world_id32, group_id32, master_id32)
   world_id = int(world_id32, int64)
   group_id = int(group_id32, int64)
   master_id = int(master_id32, int64)
end subroutine ddi_group_create

subroutine ddi_scope(comm_id) bind(c, name="ddi_scope_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_scope_impl
   integer(int64), intent(in) :: comm_id
   call ddi_scope_impl(int(comm_id, int32))
end subroutine ddi_scope

subroutine ddi_ascope(comm_id) bind(c, name="ddi_ascope_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_ascope_impl
   integer(int64), intent(in) :: comm_id
   call ddi_ascope_impl(int(comm_id, int32))
end subroutine ddi_ascope

! ============================================================================
! Timer Functions
! ============================================================================

subroutine ddi_timer_reset() bind(c, name="ddi_timer_reset_")
   use ddi_impl, only: ddi_timer_reset_impl
   call ddi_timer_reset_impl()
end subroutine ddi_timer_reset

subroutine ddi_timer_output() bind(c, name="ddi_timer_output_")
   use ddi_impl, only: ddi_timer_output_impl
   call ddi_timer_output_impl()
end subroutine ddi_timer_output

! ============================================================================
! Communicator-specific Operations
! ============================================================================

subroutine ddi_comm_create(parent_id, color, new_comm_id) bind(c, name="ddi_comm_create_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_comm_create_impl
   integer(int64), intent(in) :: parent_id, color
   integer(int64), intent(out) :: new_comm_id
   integer(int32) :: new_comm_id32
   call ddi_comm_create_impl(int(parent_id, int32), int(color, int32), new_comm_id32)
   new_comm_id = int(new_comm_id32, int64)
end subroutine ddi_comm_create

subroutine ddi_comm_destroy(comm_id) bind(c, name="ddi_comm_destroy_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_comm_destroy_impl
   integer(int64), intent(in) :: comm_id
   call ddi_comm_destroy_impl(int(comm_id, int32))
end subroutine ddi_comm_destroy

subroutine ddi_nproc_comm(comm_id, np, me) bind(c, name="ddi_nproc_comm_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_nproc_comm_impl
   integer(int64), intent(in) :: comm_id
   integer(int64), intent(out) :: np, me
   integer(int32) :: np32, me32
   call ddi_nproc_comm_impl(int(comm_id, int32), np32, me32)
   np = int(np32, int64)
   me = int(me32, int64)
end subroutine ddi_nproc_comm

subroutine ddi_sync_comm(comm_id) bind(c, name="ddi_sync_comm_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_sync_comm_impl
   integer(int64), intent(in) :: comm_id
   call ddi_sync_comm_impl(int(comm_id, int32))
end subroutine ddi_sync_comm

subroutine ddi_gsumf_comm(comm_id, buffer, n) bind(c, name="ddi_gsumf_comm_")
   use iso_fortran_env, only: int64, int32, real64
   use ddi_impl, only: ddi_gsumf_comm_impl
   integer(int64), intent(in) :: comm_id, n
   real(real64), intent(inout) :: buffer(*)
   call ddi_gsumf_comm_impl(int(comm_id, int32), buffer, int(n, int32))
end subroutine ddi_gsumf_comm

subroutine ddi_gsumi_comm(comm_id, buffer, n) bind(c, name="ddi_gsumi_comm_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_gsumi_comm_impl
   integer(int64), intent(in) :: comm_id, n
   integer(int64), intent(inout) :: buffer(*)
   ! Pass int64 buffer directly - impl now uses MPI_INTEGER8
   call ddi_gsumi_comm_impl(int(comm_id, int32), buffer, int(n, int32))
end subroutine ddi_gsumi_comm

subroutine ddi_bcast_comm(comm_id, dtype, buffer, n, root) bind(c, name="ddi_bcast_comm_")
   use iso_fortran_env, only: int64, int32, real64
   use ddi_impl, only: ddi_bcast_comm_impl
   integer(int64), intent(in) :: comm_id, n, root
   character(len=1), intent(in) :: dtype
   real(real64), intent(inout) :: buffer(*)
   call ddi_bcast_comm_impl(int(comm_id, int32), dtype, buffer, int(n, int32), int(root, int32))
end subroutine ddi_bcast_comm

! ============================================================================
! Array-wide Operations
! ============================================================================

subroutine ddi_arr_fill(handle, value) bind(c, name="ddi_arr_fill_")
   use iso_fortran_env, only: int64, int32, real64
   use ddi_impl, only: ddi_arr_fill_impl
   integer(int64), intent(in) :: handle
   real(real64), intent(in) :: value
   call ddi_arr_fill_impl(int(handle, int32), value)
end subroutine ddi_arr_fill

subroutine ddi_arr_scale(handle, scale) bind(c, name="ddi_arr_scale_")
   use iso_fortran_env, only: int64, int32, real64
   use ddi_impl, only: ddi_arr_scale_impl
   integer(int64), intent(in) :: handle
   real(real64), intent(in) :: scale
   call ddi_arr_scale_impl(int(handle, int32), scale)
end subroutine ddi_arr_scale

subroutine ddi_arr_min(handle, minval) bind(c, name="ddi_arr_min_")
   use iso_fortran_env, only: int64, int32, real64
   use ddi_impl, only: ddi_arr_min_impl
   integer(int64), intent(in) :: handle
   real(real64), intent(out) :: minval
   call ddi_arr_min_impl(int(handle, int32), minval)
end subroutine ddi_arr_min

subroutine ddi_arr_max(handle, maxval) bind(c, name="ddi_arr_max_")
   use iso_fortran_env, only: int64, int32, real64
   use ddi_impl, only: ddi_arr_max_impl
   integer(int64), intent(in) :: handle
   real(real64), intent(out) :: maxval
   call ddi_arr_max_impl(int(handle, int32), maxval)
end subroutine ddi_arr_max

subroutine ddi_arr_dot(handle_a, handle_b, dotprod) bind(c, name="ddi_arr_dot_")
   use iso_fortran_env, only: int64, int32, real64
   use ddi_impl, only: ddi_arr_dot_impl
   integer(int64), intent(in) :: handle_a, handle_b
   real(real64), intent(out) :: dotprod
   call ddi_arr_dot_impl(int(handle_a, int32), int(handle_b, int32), dotprod)
end subroutine ddi_arr_dot

subroutine ddi_arr_add(handle_a, handle_b, handle_c) bind(c, name="ddi_arr_add_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_arr_add_impl
   integer(int64), intent(in) :: handle_a, handle_b, handle_c
   call ddi_arr_add_impl(int(handle_a, int32), int(handle_b, int32), int(handle_c, int32))
end subroutine ddi_arr_add

subroutine ddi_arr_acc(handle_a, handle_b, alpha) bind(c, name="ddi_arr_acc_")
   use iso_fortran_env, only: int64, int32, real64
   use ddi_impl, only: ddi_arr_acc_impl
   integer(int64), intent(in) :: handle_a, handle_b
   real(real64), intent(in) :: alpha
   call ddi_arr_acc_impl(int(handle_a, int32), int(handle_b, int32), alpha)
end subroutine ddi_arr_acc

subroutine ddi_getacc(handle, ilo, ihi, jlo, jhi, buffer) bind(c, name="ddi_getacc_")
   use iso_fortran_env, only: int64, int32, real64
   use ddi_impl, only: ddi_getacc_impl
   integer(int64), intent(in) :: handle, ilo, ihi, jlo, jhi
   real(real64), intent(inout) :: buffer(*)
   call ddi_getacc_impl(int(handle, int32), int(ilo, int32), int(ihi, int32), &
                        int(jlo, int32), int(jhi, int32), buffer)
end subroutine ddi_getacc

subroutine ddi_group_create_custom(ngroups, group_sizes, world_id, group_id, master_id) &
   bind(c, name="ddi_group_create_custom_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_group_create_custom_impl
   integer(int64), intent(in) :: ngroups
   integer(int64), intent(in) :: group_sizes(*)
   integer(int64), intent(out) :: world_id, group_id, master_id
   integer(int32), allocatable :: group_sizes32(:)
   integer(int32) :: world_id32, group_id32, master_id32
   integer(int64) :: i
   allocate (group_sizes32(ngroups))
   do i = 1, ngroups
      group_sizes32(i) = int(group_sizes(i), int32)
   end do
   call ddi_group_create_custom_impl(int(ngroups, int32), group_sizes32, world_id32, group_id32, master_id32)
   world_id = int(world_id32, int64)
   group_id = int(group_id32, int64)
   master_id = int(master_id32, int64)
   deallocate (group_sizes32)
end subroutine ddi_group_create_custom

! ============================================================================
! SMP (Shared Memory) Operations
! ============================================================================

subroutine ddi_smp_create(size, handle) bind(c, name="ddi_smp_create_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_smp_create_impl
   integer(int64), intent(in) :: size
   integer(int64), intent(out) :: handle
   integer(int32) :: handle32
   call ddi_smp_create_impl(size, handle32)
   handle = int(handle32, int64)
end subroutine ddi_smp_create

subroutine ddi_smp_destroy(handle) bind(c, name="ddi_smp_destroy_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_smp_destroy_impl
   integer(int64), intent(in) :: handle
   call ddi_smp_destroy_impl(int(handle, int32))
end subroutine ddi_smp_destroy

subroutine ddi_smp_nproc(np, me) bind(c, name="ddi_smp_nproc_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_smp_nproc_impl
   integer(int64), intent(out) :: np, me
   integer(int32) :: np32, me32
   call ddi_smp_nproc_impl(np32, me32)
   np = int(np32, int64)
   me = int(me32, int64)
end subroutine ddi_smp_nproc

subroutine ddi_smp_offset(handle, addr_array, offset) bind(c, name="ddi_smp_offset_")
   use iso_fortran_env, only: int64, int32
   use iso_c_binding, only: c_loc, c_ptr
   use ddi_impl, only: ddi_smp_offset_impl
   integer(int64), intent(in) :: handle
   character(len=1), intent(in), target :: addr_array(*)
   integer(int64), intent(out) :: offset
   type(c_ptr) :: addr_ptr
   integer(int64) :: addr_base

   ! Get address of ADDR array
   addr_ptr = c_loc(addr_array(1))
   addr_base = transfer(addr_ptr, addr_base)

   call ddi_smp_offset_impl(int(handle, int32), addr_base, offset)
end subroutine ddi_smp_offset

subroutine ddi_smp_sync() bind(c, name="ddi_smp_sync_")
   use ddi_impl, only: ddi_smp_sync_impl
   call ddi_smp_sync_impl()
end subroutine ddi_smp_sync

subroutine ddi_smp_bcast(buffer, n, root) bind(c, name="ddi_smp_bcast_")
   use iso_fortran_env, only: int64, int32, real64
   use ddi_impl, only: ddi_smp_bcast_impl
   real(real64), intent(inout) :: buffer(*)
   integer(int64), intent(in) :: n, root
   call ddi_smp_bcast_impl(buffer, int(n, int32), int(root, int32))
end subroutine ddi_smp_bcast

subroutine ddi_smp_gsumf(buffer, n) bind(c, name="ddi_smp_gsumf_")
   use iso_fortran_env, only: int64, int32, real64
   use ddi_impl, only: ddi_smp_gsumf_impl
   real(real64), intent(inout) :: buffer(*)
   integer(int64), intent(in) :: n
   call ddi_smp_gsumf_impl(buffer, int(n, int32))
end subroutine ddi_smp_gsumf

! ============================================================================
! Masters Operations
! ============================================================================

subroutine ddi_masters_gsumf(buffer, n) bind(c, name="ddi_masters_gsumf_")
   use iso_fortran_env, only: int64, int32, real64
   use ddi_impl, only: ddi_masters_gsumf_impl
   real(real64), intent(inout) :: buffer(*)
   integer(int64), intent(in) :: n
   call ddi_masters_gsumf_impl(buffer, int(n, int32))
end subroutine ddi_masters_gsumf

subroutine ddi_masters_bcast(buffer, n, root) bind(c, name="ddi_masters_bcast_")
   use iso_fortran_env, only: int64, int32, real64
   use ddi_impl, only: ddi_masters_bcast_impl
   real(real64), intent(inout) :: buffer(*)
   integer(int64), intent(in) :: n, root
   call ddi_masters_bcast_impl(buffer, int(n, int32), int(root, int32))
end subroutine ddi_masters_bcast

! ============================================================================
! Communicator-specific Get/Put
! ============================================================================

subroutine ddi_get_comm(handle, ilo, ihi, jlo, jhi, buffer, comm_id) bind(c, name="ddi_get_comm_")
   use iso_fortran_env, only: int64, int32, real64
   use ddi_impl, only: ddi_get_comm_impl
   integer(int64), intent(in) :: handle, ilo, ihi, jlo, jhi, comm_id
   real(real64), intent(out) :: buffer(*)
   call ddi_get_comm_impl(int(handle, int32), int(ilo, int32), int(ihi, int32), &
                          int(jlo, int32), int(jhi, int32), buffer, int(comm_id, int32))
end subroutine ddi_get_comm

subroutine ddi_put_comm(handle, ilo, ihi, jlo, jhi, buffer, comm_id) bind(c, name="ddi_put_comm_")
   use iso_fortran_env, only: int64, int32, real64
   use ddi_impl, only: ddi_put_comm_impl
   integer(int64), intent(in) :: handle, ilo, ihi, jlo, jhi, comm_id
   real(real64), intent(in) :: buffer(*)
   call ddi_put_comm_impl(int(handle, int32), int(ilo, int32), int(ihi, int32), &
                          int(jlo, int32), int(jhi, int32), buffer, int(comm_id, int32))
end subroutine ddi_put_comm

subroutine ddi_get_dsid(dsid) bind(c, name="ddi_get_dsid_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_get_dsid_impl
   integer(int64), intent(out) :: dsid
   integer(int32) :: dsid32
   call ddi_get_dsid_impl(dsid32)
   dsid = int(dsid32, int64)
end subroutine ddi_get_dsid

subroutine ddi_scatter_acc(handle, buff, nelem, rows, cols) bind(c, name="ddi_scatter_acc_")
   use iso_fortran_env, only: int64, int32, real64
   use ddi_impl, only: ddi_scatter_acc_impl
   integer(int64), intent(in) :: handle, nelem
   real(real64), intent(in) :: buff(*)
   integer(int64), intent(in) :: rows(*), cols(*)
   integer(int32), allocatable :: rows32(:), cols32(:)
   integer(int64) :: i
   allocate (rows32(nelem), cols32(nelem))
   do i = 1, nelem
      rows32(i) = int(rows(i), int32)
      cols32(i) = int(cols(i), int32)
   end do
   call ddi_scatter_acc_impl(int(handle, int32), buff, int(nelem, int32), rows32, cols32)
   deallocate (rows32, cols32)
end subroutine ddi_scatter_acc

! ============================================================================
! Process-based Dynamic Load Balancing
! ============================================================================

subroutine ddi_procdlb_create(handle) bind(c, name="ddi_procdlb_create_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_procdlb_create_impl
   integer(int64), intent(out) :: handle
   integer(int32) :: handle32
   call ddi_procdlb_create_impl(handle32)
   handle = int(handle32, int64)
end subroutine ddi_procdlb_create

subroutine ddi_procdlb_destroy(handle) bind(c, name="ddi_procdlb_destroy_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_procdlb_destroy_impl
   integer(int64), intent(in) :: handle
   call ddi_procdlb_destroy_impl(int(handle, int32))
end subroutine ddi_procdlb_destroy

subroutine ddi_procdlb_reset(handle) bind(c, name="ddi_procdlb_reset_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_procdlb_reset_impl
   integer(int64), intent(in) :: handle
   call ddi_procdlb_reset_impl(int(handle, int32))
end subroutine ddi_procdlb_reset

subroutine ddi_procdlb_next(handle, iproc, counter) bind(c, name="ddi_procdlb_next_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_procdlb_next_impl
   integer(int64), intent(in) :: handle, iproc
   integer(int64), intent(out) :: counter
   call ddi_procdlb_next_impl(int(handle, int32), int(iproc, int32), counter)
end subroutine ddi_procdlb_next

! ============================================================================
! Miscellaneous
! ============================================================================

subroutine ddi_ndistrib(handle, node, ilo, ihi, jlo, jhi) bind(c, name="ddi_ndistrib_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_ndistrib_impl
   integer(int64), intent(in) :: handle, node
   integer(int64), intent(out) :: ilo, ihi, jlo, jhi
   integer(int32) :: ilo32, ihi32, jlo32, jhi32
   call ddi_ndistrib_impl(int(handle, int32), int(node, int32), ilo32, ihi32, jlo32, jhi32)
   ilo = int(ilo32, int64)
   ihi = int(ihi32, int64)
   jlo = int(jlo32, int64)
   jhi = int(jhi32, int64)
end subroutine ddi_ndistrib

subroutine ddi_level(level) bind(c, name="ddi_level_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_level_impl
   integer(int64), intent(out) :: level
   integer(int32) :: level32
   call ddi_level_impl(level32)
   level = int(level32, int64)
end subroutine ddi_level

! ============================================================================
! No-op Stubs (rarely used internal/debug functions)
! ============================================================================

subroutine ddi_pend() bind(c, name="ddi_pend_")
   use ddi_impl, only: ddi_pend_impl
   call ddi_pend_impl()
end subroutine ddi_pend

subroutine ddi_pbeg(section) bind(c, name="ddi_pbeg_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_pbeg_impl
   integer(int64), intent(in) :: section
   call ddi_pbeg_impl(int(section, int32))
end subroutine ddi_pbeg

subroutine ddi_addr_test(addr, valid) bind(c, name="ddi_addr_test_")
   use iso_fortran_env, only: int64, int32
   use ddi_impl, only: ddi_addr_test_impl
   integer(int64), intent(in) :: addr
   integer(int64), intent(out) :: valid
   integer(int32) :: valid32
   call ddi_addr_test_impl(addr, valid32)
   valid = int(valid32, int64)
end subroutine ddi_addr_test

subroutine ddi_gsum(tag, dtype, buffer, n) bind(c, name="ddi_gsum_")
   use iso_fortran_env, only: int64, int32, real64
   use ddi_impl, only: ddi_gsum_impl
   integer(int64), intent(in) :: tag, n
   character(len=1), intent(in) :: dtype
   real(real64), intent(inout) :: buffer(*)
   call ddi_gsum_impl(int(tag, int32), dtype, buffer, int(n, int32))
end subroutine ddi_gsum

subroutine ddi_gdlbreset_inner() bind(c, name="ddi_gdlbreset_inner_")
   use ddi_impl, only: ddi_gdlbreset_inner_impl
   call ddi_gdlbreset_inner_impl()
end subroutine ddi_gdlbreset_inner

subroutine ddi_gdlbnext_inner(counter) bind(c, name="ddi_gdlbnext_inner_")
   use iso_fortran_env, only: int64
   use ddi_impl, only: ddi_gdlbnext_inner_impl
   integer(int64), intent(out) :: counter
   call ddi_gdlbnext_inner_impl(counter)
end subroutine ddi_gdlbnext_inner
