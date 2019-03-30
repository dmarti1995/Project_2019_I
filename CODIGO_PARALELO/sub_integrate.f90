! INTEGRACION: VELOCITY VERLET
! ----------------------------

! INPUT:
!       + N               :: Numero de particulas
!       + dt              :: paso de tiempo
!       + Rcut            :: cut-off radio de las interacciones
!       + L               :: lado de la caja
!       + numproc         :: numero de procesadores
!       + taskid          :: id del procesador que esta ejectudando
!       + max_length      :: cantidad maxima de elementos con los que trabaja un procesador
!       + displs(numproc) :: indices desde los cuales comienza a trabajar cada procesador
!       + counts(numproc) :: cantidad de elementos con los que trabaja cada procesador
!       + index_local(numproc,2) :: tabla en la que se indica cuales son los indices entre los que
!                                   trabaja cada procesador (respecto a la matriz global)
! INPUT y OUTPUT:
!       + r(N,3) :: posiciones de las particulas
!       + v(N,3) :: velocidades en el instante t
!       + F(N,3) :: fuerzas sobre cada particula (Entra la fuerza a tiempo t y sale la fuerza a t+dt)
!       + t      :: tiempo
!       + ierror :: parametro para controlar errores de MPI
! OUTPUT:
!       + press  :: presion
!       + Epot   :: energia potencial
!       + Ekin   :: energia cinetica

subroutine v_verlet_step (N, r, v, t, dt, Rcut, L, F, press, Epot, Ekin, &
                          numproc, taskid, max_length, displs, counts, index_local, ierror)
use mpi
implicit none

integer, parameter              :: dp = 8
real(dp), intent(inout)         :: r(N,3), v(N,3), F(N,3), t
real(dp), intent(in)            :: dt, Rcut, L
real(dp), intent(out)           :: press, Epot, Ekin
integer, intent(in)             :: N
integer, intent(in)             :: numproc, taskid, max_length
integer, intent(in)             :: displs(numproc), counts(numproc), index_local(numproc,2)
integer, intent(inout)          :: ierror
real(dp), dimension(max_length) :: local_r, local_v, local_f
integer                         :: my_N_elem
integer                         :: i, j, k

my_N_elem = counts(taskid+1)

! CREAR SUBVECTOR LOCAL
! ACTUALIZACION 1: POSICIONES Y VELOCIDAD PROVISIONAL
i = index_local(taskid+1, 1)
j = index_local(taskid+1, 2)
do k = 1, my_N_elem
  local_r(k) = r(i,j) + v(i,j)*dt + 0.5_dp*F(i,j)*dt*dt
  local_v(k) = v(i,j) + 0.5_dp*F(i,j)*dt

  i = i + 1
  if (i > N) then
    i = 1
    j = j + 1
  endif
enddo


! COMUNICACION ENTRE PROCESADORES: POSICIONES
call MPI_ALLGATHERV (local_r, my_N_elem, MPI_REAL8,  &
                     r, counts, displs, MPI_REAL8,   &
                     MPI_COMM_WORLD, ierror)


! CALCULO FUERZAS
call force (N, L, Rcut, r, v, F, press, Epot, Ekin)


! CREAR SUBVECTOR LOCAL FUERZAS
! ACTUALIZACION 2: VELOCIDADES
i = index_local(taskid+1, 1)
j = index_local(taskid+1, 2)
do k = 1, my_N_elem
  local_f(k) = F(i,j)
  local_v(k) = local_v(k) + 0.5_dp*F(i,j)*dt
  i = i + 1
  if (i > N) then
    i = 1
    j = j + 1
  endif
enddo


! COMUNICACION ENTRE PROCESADORES: VELOCIDADES
call MPI_ALLGATHERV (local_v, my_N_elem, MPI_REAL8,  &
                     v, counts, displs, MPI_REAL8,   &
                     MPI_COMM_WORLD, ierror)


t = t + dt

end subroutine v_verlet_step



! TERMOSTATO DE ANDERSEN
! ----------------------

subroutine thermostat (N, v, nu, temperature,numproc,taskid,max_length, &
displs, counts, index_local, ierror)
use mpi
implicit none
integer, parameter  :: dp = 8
real(dp), intent(inout) :: v(N,3)
real(dp), intent(in)    :: nu, temperature
integer, intent(in) :: N,numproc,taskid
integer, intent(in)         :: max_length
integer, intent(in)         :: displs(numproc), counts(numproc), index_local(numproc,2)
real(dp), dimension(max_length) :: local_r, local_v, local_f

integer ::             request(numproc),ierror,stat(MPI_STATUS_SIZE),partner
real(dp)                :: u, z1, z2
integer             :: i, j,my_N_elem,k



my_N_elem = counts(taskid+1)

i = index_local(taskid+1, 1)
j = index_local(taskid+1, 2)


!do i = 1, N
do k=1,my_N_elem         !each
    local_v(k) = v(i,j)
    call random_number(u)
    if (u < nu) then
        call normal_distr(0.0_dp, sqrt(temperature), z1, z2)
        local_v(k) = z1                  !changing velocities to a gaussian distribution, randomly
        local_v(k) = z2
        call normal_distr(0.0_dp, sqrt(temperature), z1, z2)
        local_v(k) = z1
    endif

    i = i + 1
    if (i > N) then
        i = 1
        j = j + 1
    endif


enddo


call MPI_ALLGATHERV (local_v, my_N_elem, MPI_REAL8,  &
v, counts, displs, MPI_REAL8,   &
MPI_COMM_WORLD, ierror)



end subroutine thermostat




! GENERADOR DE DISTRIBUCION GAUSSIANA
! -----------------------------------

subroutine normal_distr(mu, sigma, z1, z2)
implicit none 
integer, parameter    :: dp = 8
real(dp), intent(in)  :: mu, sigma
real(dp), intent(out) :: z1, z2
real(dp), parameter   :: pi = acos(-1.0)
real(dp)              :: x, y

call random_number(x)
call random_number(y)
z1 = sqrt(-2.0_dp * log(x)) * cos(2.0_dp * pi * y)
z2 = sqrt(-2.0_dp * log(x)) * sin(2.0_dp * pi * y)
z1 = z1 * sigma + mu
z2 = z2 * sigma + mu

end subroutine normal_distr

