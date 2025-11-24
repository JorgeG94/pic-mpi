program hierarchical_mpi_test_legacy
   use mpi, only: MPI_MAX_PROCESSOR_NAME, MPI_Init, MPI_Finalize
   use mpi_comm_simple_legacy, only: comm_t, comm_world, abort_comm, allgather, &
                                     get_processor_name, iprobe, recv, send, &
                                     MPI_Status, MPI_ANY_SOURCE, MPI_ANY_TAG
   use pic_timer, only: timer_type
   use pic_types, only: dp, default_int
   implicit none

   ! MPI wrappers
   type(comm_t) :: world_comm, node_comm
   character(len=MPI_MAX_PROCESSOR_NAME) :: hostname
   integer :: hostname_len

   ! Timing
   type(timer_type) :: timer
   real(dp) :: elapsed_time

   ! Node leader detection
   integer :: global_node_rank
   integer, allocatable :: all_node_leader_ranks(:), node_leader_ranks(:)
   integer :: i, j, num_nodes

   ! Test parameters
   integer(default_int) :: n_tasks, matrix_size
   integer :: narg
   character(len=32) :: arg

   !==============================
   ! MPI Initialization
   !==============================
   call MPI_Init(i)
   world_comm = comm_world()
   node_comm = world_comm%split()   ! shared memory communicator

   ! Get command line arguments
   narg = command_argument_count()
   n_tasks = 100      ! number of tasks to distribute
   matrix_size = 512  ! size for computation

   if (narg >= 1) then
      call get_command_argument(1, arg)
      read (arg, *) n_tasks
   end if
   if (narg >= 2) then
      call get_command_argument(2, arg)
      read (arg, *) matrix_size
   end if

   if (world_comm%size() < 2) then
      if (world_comm%leader()) print *, "This program requires at least 2 processes (LEGACY VERSION)."
      call abort_comm(world_comm, 1)
   end if

   if (world_comm%leader()) then
      call timer%start()
      print *, "========================================"
      print *, "Hierarchical MPI Test (LEGACY VERSION)"
      print *, "========================================"
      print *, "Total MPI ranks: ", world_comm%size()
      print *, "Tasks to distribute: ", n_tasks
      print *, "Matrix size: ", matrix_size
   end if

   !==============================
   ! Determine node leaders
   !==============================
   global_node_rank = -1
   if (node_comm%leader()) global_node_rank = world_comm%rank()

   allocate (all_node_leader_ranks(world_comm%size()))
   call allgather(world_comm, global_node_rank, all_node_leader_ranks)

   num_nodes = count(all_node_leader_ranks /= -1)

   ! Collect actual node leader ranks
   allocate (node_leader_ranks(num_nodes))
   i = 0
   do j = 1, world_comm%size()
      if (all_node_leader_ranks(j) /= -1) then
         i = i + 1
         node_leader_ranks(i) = all_node_leader_ranks(j)
      end if
   end do

   if (world_comm%leader()) then
      print *, "Number of nodes detected: ", num_nodes
      print *, "Node leader ranks: ", node_leader_ranks
   end if

   ! Multi-node check
   if (num_nodes > 1) then
      if (world_comm%size() < 3) then
         if (world_comm%leader()) print *, "Multi-node requires at least 3 processes."
         call abort_comm(world_comm, 1)
      end if
   end if

   !==============================
   ! Get processor name and setup devices
   !==============================
   call get_processor_name(hostname, hostname_len)

   if (node_comm%leader()) then
      print *, "Node leader on: ", trim(hostname), &
         " | node rank: ", node_comm%rank(), &
         " | world rank: ", world_comm%rank()
   else
      print *, "Node worker on: ", trim(hostname), &
         " | node rank: ", node_comm%rank(), &
         " | world rank: ", world_comm%rank()
      ! Set device for workers (commented out to avoid external dependencies)
      !call omp_set_default_device(node_comm%rank() - 1)
      !$acc set device_num(node_comm%rank() - 1)
   end if

   call world_comm%barrier()

   !==============================
   ! Role assignment and task distribution
   !==============================
   if (world_comm%leader() .and. node_comm%leader()) then
      print *, "========================================"
      print *, "Starting global coordinator"
      print *, "========================================"
      ! Global coordinator (rank 0 on node 0)
      call test_global_coordinator(world_comm, node_comm, n_tasks, node_leader_ranks, num_nodes, matrix_size)
   else if (node_comm%leader()) then
      ! Node coordinator (rank 0 on other nodes)
      call test_node_coordinator(world_comm, node_comm, matrix_size)
   else
      ! Worker process
      call test_node_worker(world_comm, node_comm, matrix_size)
   end if

   !==============================
   ! Final timing
   !==============================
   call world_comm%barrier()
   if (world_comm%leader()) then
      call timer%stop()
      elapsed_time = timer%get_elapsed_time()
      print *, "========================================"
      print *, "Total elapsed time: ", elapsed_time, " seconds"
      print *, "Tasks per second: ", real(n_tasks, dp)/elapsed_time
      print *, "========================================"
   end if

   !==============================
   ! Cleanup
   !==============================
   deallocate (all_node_leader_ranks, node_leader_ranks)
   call node_comm%finalize()
   call world_comm%finalize()
   call MPI_Finalize(i)

contains

   subroutine test_global_coordinator(world_comm, node_comm, total_tasks, node_leaders, num_nodes, matrix_size)
      type(comm_t), intent(in) :: world_comm, node_comm
      integer, intent(in) :: total_tasks, num_nodes, matrix_size
      integer, intent(in) :: node_leaders(:)

      integer :: current_task, finished_nodes
      integer :: request_source, dummy_msg
      type(MPI_Status) :: status, local_status
      logical :: handling_local_workers
      logical :: has_pending
      integer :: local_finished_workers, local_dummy

      current_task = total_tasks
      finished_nodes = 0
      local_finished_workers = 0
      handling_local_workers = (node_comm%size() > 1)

      print *, "Global coordinator managing ", total_tasks, " tasks for ", num_nodes, " nodes"

      do while (finished_nodes < num_nodes)
         ! Check for remote node coordinator requests
         call iprobe(world_comm, MPI_ANY_SOURCE, 300, has_pending, status)
         if (has_pending) then
            call recv(world_comm, dummy_msg, status%MPI_SOURCE, 300)
            request_source = status%MPI_SOURCE

            if (current_task >= 1) then
               ! Send task to remote node coordinator
               call send(world_comm, current_task, request_source, 301)
               current_task = current_task - 1
            else
               ! No more tasks, send termination signal
               call send(world_comm, -1, request_source, 302)
               finished_nodes = finished_nodes + 1
            end if
         end if

         ! Handle local workers (shared memory)
         if (handling_local_workers .and. local_finished_workers < node_comm%size() - 1) then
            call iprobe(node_comm, MPI_ANY_SOURCE, 200, has_pending, local_status)
            if (has_pending) then
               call recv(node_comm, local_dummy, local_status%MPI_SOURCE, 200)

               if (current_task >= 1) then
                  ! Send task to local worker
                  call send(node_comm, current_task, local_status%MPI_SOURCE, 201)
                  call send(node_comm, matrix_size, local_status%MPI_SOURCE, 201)
                  current_task = current_task - 1
               else
                  ! No more tasks
                  call send(node_comm, -1, local_status%MPI_SOURCE, 202)
                  local_finished_workers = local_finished_workers + 1
               end if
            end if
         end if

         ! Check if all local workers finished
         if (handling_local_workers .and. local_finished_workers >= node_comm%size() - 1) then
            handling_local_workers = .false.
            finished_nodes = finished_nodes + 1
            print *, "Global coordinator: all local workers finished"
         end if
      end do

      print *, "Global coordinator: all tasks distributed"
   end subroutine test_global_coordinator

   subroutine test_node_coordinator(world_comm, node_comm, matrix_size)
      type(comm_t), intent(in) :: world_comm, node_comm
      integer, intent(in) :: matrix_size

      integer :: task_id, dummy_msg
      integer :: finished_workers
      integer :: local_dummy
      type(MPI_Status) :: status, global_status
      logical :: local_message_pending, more_tasks

      finished_workers = 0
      more_tasks = .true.
      dummy_msg = 0

      do while (finished_workers < node_comm%size() - 1)
         ! Check for local worker requests
         call iprobe(node_comm, MPI_ANY_SOURCE, 200, local_message_pending, status)

         if (local_message_pending) then
            call recv(node_comm, local_dummy, MPI_ANY_SOURCE, 200)

            if (more_tasks) then
               ! Request task from global coordinator
               call send(world_comm, dummy_msg, 0, 300)
               call recv(world_comm, task_id, 0, MPI_ANY_TAG, global_status)

               if (global_status%MPI_TAG == 301) then
                  ! Forward task to local worker
                  call send(node_comm, task_id, status%MPI_SOURCE, 201)
                  call send(node_comm, matrix_size, status%MPI_SOURCE, 201)
               else
                  ! No more tasks
                  call send(node_comm, -1, status%MPI_SOURCE, 202)
                  finished_workers = finished_workers + 1
                  more_tasks = .false.
               end if
            else
               ! Already out of tasks
               call send(node_comm, -1, status%MPI_SOURCE, 202)
               finished_workers = finished_workers + 1
            end if
         end if
      end do
   end subroutine test_node_coordinator

   subroutine test_node_worker(world_comm, node_comm, matrix_size)
      type(comm_t), intent(in) :: world_comm, node_comm
      integer, intent(in) :: matrix_size

      integer :: task_id, task_size, dummy_msg
      type(MPI_Status) :: status
      real(dp), allocatable :: A(:, :), B(:, :), C(:, :)
      integer :: i, j, k, dims
      real(dp), parameter :: alpha = 17.0_dp

      dummy_msg = 0

      do
         ! Request work from node coordinator
         call send(node_comm, dummy_msg, 0, 200)
         call recv(node_comm, task_id, 0, MPI_ANY_TAG, status)

         select case (status%MPI_TAG)
         case (201)
            ! Received task
            call recv(node_comm, task_size, 0, 201, status)

            ! Simulate work with matrix multiplication
            dims = task_size
            allocate (A(dims, dims), B(dims, dims), C(dims, dims))

            ! Initialize matrices
            do concurrent(j=1:dims, i=1:dims)
               A(i, j) = real(task_id, dp)
               B(i, j) = real(task_id, dp)
               C(i, j) = 0.0_dp
            end do

            ! Simple matrix multiplication (naive implementation - jaxpy)
            do concurrent(j=1:dims, i=1:dims)
               C(i, j) = alpha*A(i, j) + B(i, j)
            end do

            deallocate (A, B, C)

         case (202)
            ! Termination signal
            exit
         end select
      end do
   end subroutine test_node_worker

end program hierarchical_mpi_test_legacy
