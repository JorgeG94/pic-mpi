module pic_legacy_mpi_algorithms
   use mpi, only: MPI_STATUS_SIZE, MPI_ANY_SOURCE, MPI_ANY_TAG, &
                  MPI_SOURCE, MPI_TAG
   use mpi_comm_simple_legacy, only: comm_t, iprobe, recv, send
   use pic_types, only: dp, int32
   use pic_timer, only: timer_type
   implicit none
   private

   public :: process_fragment, calculate_exact_flops
   public :: global_coordinator, node_coordinator, node_worker
   public :: send_fragment_to_node, send_fragment_to_worker

contains

   subroutine process_fragment(fragment_idx, fragment_indices, fragment_size, matrix_size)
      integer, intent(in) :: fragment_idx, fragment_size, matrix_size
      integer, intent(in) :: fragment_indices(fragment_size)
      real(dp), allocatable :: A(:, :), B(:, :), C(:, :)
      integer :: i, j, k
      type(timer_type) :: gemm_timer
      real(dp) :: elapsed_time
      real(dp) :: sum_val
      integer :: dims
      integer :: error
      integer :: istat
      real(dp), parameter :: alpha = 1.0_dp
      real(dp), parameter :: beta = 0.0_dp

      dims = fragment_size*matrix_size

      ! Allocate and initialize fragment matrix
      allocate (A(dims, dims), B(dims, dims), C(dims, dims))

      do concurrent(j=1:dims, i=1:dims)
         A(i, j) = real(fragment_size*fragment_idx, dp)
         B(i, j) = real(fragment_size*fragment_idx, dp)
         C(i, j) = 0.0_dp
      end do

      call gemm_timer%start()
      call gemm_timer%stop()
      elapsed_time = gemm_timer%get_elapsed_time()

      !print *, "Gemm for fragment", fragment_indices, " was ", elapsed_time, " seconds with size ", dims

      deallocate (A, B, C)
   end subroutine process_fragment

   subroutine calculate_exact_flops(polymers, fragment_count, max_level, matrix_size, total_flops)
      integer, intent(in) :: polymers(:, :), fragment_count, max_level, matrix_size
      real(dp), intent(out) :: total_flops

      integer :: i, fragment_size
      integer, allocatable :: n_mers(:)
      real(dp), allocatable :: mer_flops(:)
      real(dp) :: mer_size

      ! Allocate counters for each n-mer level (1 to max_level)
      allocate (n_mers(max_level))
      allocate (mer_flops(max_level))

      n_mers = 0
      mer_flops = 0.0_dp

      ! Count fragments by size
      do i = 1, fragment_count
         fragment_size = count(polymers(i, :) > 0)
         if (fragment_size >= 1 .and. fragment_size <= max_level) then
            n_mers(fragment_size) = n_mers(fragment_size) + 1
         end if
      end do

      ! Calculate FLOPs for each n-mer level
      do i = 1, max_level
         mer_size = real(i*matrix_size, dp)
         mer_flops(i) = real(n_mers(i), dp)*2.0_dp*mer_size**3
      end do

      ! Total FLOPs
      total_flops = sum(mer_flops)

      ! Print breakdown
      print *, "Fragment breakdown:"
      do i = 1, max_level
         if (n_mers(i) > 0) then
            print '(a,i0,a,i0,a,f12.3,a)', "  ", i, "-mers:  ", n_mers(i), &
               " (", mer_flops(i)/1.0e9_dp, " GFLOP)"
         end if
      end do
      print '(a,i0,a,f12.3,a)', "  Total:    ", fragment_count, &
         " (", total_flops/1.0e9_dp, " GFLOP)"

      deallocate (n_mers, mer_flops)

   end subroutine calculate_exact_flops

   subroutine global_coordinator(world_comm, node_comm, total_fragments, polymers, max_level, node_leader_ranks, num_nodes)
      type(comm_t), intent(in) :: world_comm, node_comm
      integer, intent(in) :: total_fragments, max_level, num_nodes
      integer, intent(in) :: polymers(:, :), node_leader_ranks(:)

      integer :: current_fragment, finished_nodes
      integer :: request_source, dummy_msg
      integer :: status(MPI_STATUS_SIZE), local_status(MPI_STATUS_SIZE)
      logical :: handling_local_workers
      logical :: has_pending

      ! For local workers
      integer :: local_finished_workers, fragment_size, local_dummy
      integer, allocatable :: fragment_indices(:)

      current_fragment = total_fragments
      finished_nodes = 0
      local_finished_workers = 0
      handling_local_workers = (node_comm%size() > 1)

      print *, "Global coordinator starting with", total_fragments, "fragments for", num_nodes, "nodes"

      do while (finished_nodes < num_nodes)

         ! Remote node coordinator requests
         call iprobe(world_comm, MPI_ANY_SOURCE, 300, has_pending, status)
         if (has_pending) then
            call recv(world_comm, dummy_msg, status(MPI_SOURCE), 300)
            request_source = status(MPI_SOURCE)

            if (current_fragment >= 1) then
               call send_fragment_to_node(world_comm, current_fragment, polymers, max_level, request_source)
               current_fragment = current_fragment - 1
            else
               call send(world_comm, -1, request_source, 302)
               finished_nodes = finished_nodes + 1
            end if
         end if

         ! Local workers (shared memory)
         if (handling_local_workers .and. local_finished_workers < node_comm%size() - 1) then
            call iprobe(node_comm, MPI_ANY_SOURCE, 200, has_pending, local_status)
            if (has_pending) then
               call recv(node_comm, local_dummy, local_status(MPI_SOURCE), 200)

               if (current_fragment >= 1) then
                  call send_fragment_to_worker(node_comm, current_fragment, polymers, max_level, local_status(MPI_SOURCE))
                  current_fragment = current_fragment - 1
               else
                  call send(node_comm, -1, local_status(MPI_SOURCE), 202)
                  local_finished_workers = local_finished_workers + 1
               end if
            end if
         end if

         ! Finalize local worker completion
         if (handling_local_workers .and. local_finished_workers >= node_comm%size() - 1) then
            handling_local_workers = .false.
            if (num_nodes == 1) then
               finished_nodes = finished_nodes + 1
               print *, "Manually incremented finished_nodes for self"
            else
               finished_nodes = finished_nodes + 1
               print *, "Global coordinator finished local workers"
            end if
         end if

         !call sleep(0)
      end do

      print *, "Global coordinator finished all fragments"
   end subroutine global_coordinator

   subroutine send_fragment_to_node(world_comm, fragment_idx, polymers, max_level, dest_rank)
      type(comm_t), intent(in) :: world_comm
      integer, intent(in) :: fragment_idx, max_level, dest_rank
      integer, intent(in) :: polymers(:, :)
      integer :: fragment_size
      integer, allocatable :: fragment_indices(:)

      fragment_size = count(polymers(fragment_idx, :) > 0)
      allocate (fragment_indices(fragment_size))
      fragment_indices = polymers(fragment_idx, 1:fragment_size)

      call send(world_comm, fragment_idx, dest_rank, 301)
      call send(world_comm, fragment_size, dest_rank, 301)
      call send(world_comm, fragment_indices, dest_rank, 301)

      deallocate (fragment_indices)
   end subroutine send_fragment_to_node

   subroutine send_fragment_to_worker(node_comm, fragment_idx, polymers, max_level, dest_rank)
      type(comm_t), intent(in) :: node_comm
      integer, intent(in) :: fragment_idx, max_level, dest_rank
      integer, intent(in) :: polymers(:, :)
      integer :: fragment_size
      integer, allocatable :: fragment_indices(:)

      fragment_size = count(polymers(fragment_idx, :) > 0)
      allocate (fragment_indices(fragment_size))
      fragment_indices = polymers(fragment_idx, 1:fragment_size)

      call send(node_comm, fragment_idx, dest_rank, 201)
      call send(node_comm, fragment_size, dest_rank, 201)
      call send(node_comm, fragment_indices, dest_rank, 201)

      deallocate (fragment_indices)
   end subroutine send_fragment_to_worker

   subroutine node_coordinator(world_comm, node_comm, max_level)
      class(comm_t), intent(in) :: world_comm, node_comm
      integer(int32), intent(in) :: max_level

      integer(int32) :: fragment_idx, fragment_size, ierr, dummy_msg
      integer(int32) :: finished_workers
      integer(int32), allocatable :: fragment_indices(:)
      integer :: status(MPI_STATUS_SIZE), global_status(MPI_STATUS_SIZE)
      logical :: local_message_pending, more_fragments
      integer(int32) :: local_dummy

      finished_workers = 0
      more_fragments = .true.
      dummy_msg = 0

      do while (finished_workers < node_comm%size() - 1)
         call iprobe(node_comm, MPI_ANY_SOURCE, 200, local_message_pending, status)

         if (local_message_pending) then
            call recv(node_comm, local_dummy, MPI_ANY_SOURCE, 200)

            if (more_fragments) then
               call send(world_comm, dummy_msg, 0, 300)
               call recv(world_comm, fragment_idx, 0, MPI_ANY_TAG, global_status)

               if (global_status(MPI_TAG) == 301) then
                  call recv(world_comm, fragment_size, 0, 301, global_status)
                  allocate (fragment_indices(fragment_size))
                  call recv(world_comm, fragment_indices, 0, 301, global_status)

                  call send(node_comm, fragment_idx, status(MPI_SOURCE), 201)
                  call send(node_comm, fragment_size, status(MPI_SOURCE), 201)
                  call send(node_comm, fragment_indices, status(MPI_SOURCE), 201)

                  deallocate (fragment_indices)
               else
                  call send(node_comm, -1, status(MPI_SOURCE), 202)
                  finished_workers = finished_workers + 1
                  more_fragments = .false.
               end if
            else
               call send(node_comm, -1, status(MPI_SOURCE), 202)
               finished_workers = finished_workers + 1
            end if
         end if

         !call sleep(0)
      end do
   end subroutine node_coordinator

   subroutine node_worker(world_comm, node_comm, matrix_size, max_level)

      class(comm_t), intent(in) :: world_comm, node_comm
      integer, intent(in) :: matrix_size, max_level

      integer(int32) :: fragment_idx, fragment_size, ierr, dummy_msg
      integer(int32), allocatable :: fragment_indices(:)
      integer :: status(MPI_STATUS_SIZE)

      dummy_msg = 0

      do
         call send(node_comm, dummy_msg, 0, 200)
         call recv(node_comm, fragment_idx, 0, MPI_ANY_TAG, status)

         select case (status(MPI_TAG))
         case (201)
            call recv(node_comm, fragment_size, 0, 201, status)
            allocate (fragment_indices(fragment_size))
            call recv(node_comm, fragment_indices, 0, 201, status)

            call process_fragment(fragment_idx, fragment_indices, fragment_size, matrix_size)

            deallocate (fragment_indices)
         case (202)
            exit
         end select
      end do
   end subroutine node_worker

end module pic_legacy_mpi_algorithms
