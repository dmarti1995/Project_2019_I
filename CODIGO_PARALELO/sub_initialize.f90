
  ! N       : number of particles 
  ! T       : temperatute
  ! v       : velocities of particles --> (N,3) array
  ! numproc : number of processors
  ! taskid  : integer to identify each processor
  ! ierr    : error identification parameter

  subroutine initialize_v(v,T,N,numproc,taskid,ierr)

  use mpi
  
  implicit none

  integer, parameter                    :: dp = 8
  real(dp), dimension(N,3), intent(out) :: v
  real(dp), intent(in)                  :: T
  integer, intent(in)                   :: N, numproc, taskid
  integer, intent(inout)                :: ierr
  real(dp), dimension(:,:), allocatable :: local
  integer, dimension(2)                 :: starts, size_total, size_sub
  integer, dimension(:), allocatable    :: seed
  integer                               :: blocktype, resizedtype, realsize, long, i, partition, size_seed
  real(dp)                              :: kinetic_partial, kinetic, rand(3), suma_local(3), suma(3)
  integer(kind=MPI_Address_kind)        :: start, extent

  ! For random number generator...
  call random_seed(size=size_seed)
  allocate(seed(size_seed))

  ! Each processor works with 'partition' number of particles.
  ! In general (N / numproc /= integer_number), we will need 
  ! to check and take into account the rest of particles...
  ! That will be done later
  partition = N / numproc
  allocate(local(partition,3))  
  
  starts = (/0,0/)
  size_total = (/N,3/)
  size_sub = (/partition,3/)

  ! ____________________________________________________________

  ! Create subarray structure
  call MPI_Type_create_subarray(2, size_total, size_sub, starts, &
                               MPI_ORDER_FORTRAN, MPI_REAL8,     &
                               blocktype, ierr)

  ! Realsize
  call MPI_Type_size(MPI_REAL8, realsize, ierr)

  ! Memory of each subarray
  extent = realsize*partition

  start = 0

  ! Create resized type array
  call MPI_Type_create_resized(blocktype, start, extent, resizedtype, ierr)
  call MPI_Type_commit(resizedtype, ierr)

  ! _____________________________________________________________

  ! Each processor creates random numbers for their
  ! local array : independent random numbers for each processor

  seed = 377*(1+taskid)
  call random_seed(put=seed)
  call random_number(local)
  local = local - 0.5_dp


  suma_local(1) = sum(local(:,1))
  suma_local(2) = sum(local(:,2))
  suma_local(3) = sum(local(:,3))


  ! All processor must know the sum of all the velocities
  ! for x,y and z directions to set the total momentum 
  ! to zero

  call MPI_allreduce(suma_local,suma,size(suma),MPI_REAL8, &
                     MPI_SUM,MPI_COMM_WORLD,ierr) 


  ! If the number of particles cannot be divided
  ! by the number of processor, we must complete 
  ! the v array 

  long = partition * numproc

  if (long/=N ) then
     do i=long+1,N
         ! same seed for all processors
         ! if not, we would have a different v
         ! array in each processor
         seed = 134
         call random_seed(put=seed)
         call random_number(rand)
         v(i,:) = rand - 0.5_dp
     enddo

     suma = suma + (/ sum(v(long+1:,1)), sum(v(long+1:,2)), sum(v(long+1:,3)) /)

  endif


  ! The processors set the momentum of their local vectors to zero 
  
  local(:,1) = local(:,1) - suma(1) / N
  local(:,2) = local(:,2) - suma(2) / N
  local(:,3) = local(:,3) - suma(3) / N

  ! Once again, if the number of particles cannot be divided
  ! by the number of processor, we must complete 
  ! the v array --> set momentum to zero

   if (long/=N ) then
        do i = long+1,N
           v(i,1) = v(i,1) - suma(1) / N 
           v(i,2) = v(i,2) - suma(2) / N 
           v(i,3) = v(i,3) - suma(3) / N
        enddo
   endif
  
  ! Calculate the kinetic energy: each processor
  ! makes a partial sum and then we collect everything
  ! with MPI_allreduce ---> mpi_sum

  kinetic_partial = sum(local**2.0_dp)

  call MPI_allreduce(kinetic_partial,kinetic,1,MPI_REAL8, &
                     MPI_SUM,MPI_COMM_WORLD,ierr)

  ! Check if there is more to sum ---> long/=N ?
  if (long/=N ) then
     kinetic = kinetic + sum(v(long+1:,:)**2.0_dp)
  endif

  
  kinetic = kinetic * 0.5_dp
 
  ! Rescale velocities with temperature
  local = local * sqrt(1.5_dp * N * T / kinetic)

  ! Before finishing the subroutine, all the processors
  ! must know the v array
  call MPI_allgather(local, size(local), MPI_REAL8,  &
               v, 1, resizedtype, MPI_COMM_WORLD, ierr) 

  ! Check long/=N ...
  if (long/=N ) then
         v(long+1:,:) = v(long+1:,:) * sqrt(1.5_dp * N * T / kinetic)
  endif


  deallocate(local)

  end subroutine initialize_v


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! N   : number of particles 
  ! a   : lattice spacing
  ! T   : temperatute
  ! pos : position of particles   --> (N,3) array
  ! v   : velocities of particles --> (N,3) array
  ! rho : density of the system

  subroutine initialize_r(pos,rho,N,L)

  implicit none
  integer, parameter    :: dp = 8
  integer, intent(in)   :: N
  real(dp), intent(in)  :: rho
  real(dp), intent(out) :: pos(N,3)
  real(dp), intent(out) :: L
  real(dp)              :: a, kinetic
  integer               :: i, j, k, cont, M
  integer, dimension(3) :: R

  ! Initialization of the structure

  M = nint((N/4.0_dp)**(1.0_dp/3.0_dp))
  L = (N/rho)**(1.0_dp/3.0_dp)
  a = L/M 

  cont = 0

  do i = 0, M-1
      do j = 0, M-1
          do k = 0, M-1
              R = (/i,j,k/)
              pos(4*cont+1,:) = a*R
              pos(4*cont+2,:) = a*R + pos(1,:) + a*0.5_dp*(/1.0_dp,1.0_dp,0.0_dp/)
              pos(4*cont+3,:) = a*R + pos(1,:) + a*0.5_dp*(/0.0_dp,1.0_dp,1.0_dp/)
              pos(4*cont+4,:) = a*R + pos(1,:) + a*0.5_dp*(/1.0_dp,0.0_dp,1.0_dp/)
              cont = cont + 1
          enddo
      enddo
  enddo

end subroutine initialize_r
