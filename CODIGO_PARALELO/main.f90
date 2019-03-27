program DinMo

use mpi

implicit none


real, allocatable :: pos(:,:),vel(:,:),f_par(:,:)
real, allocatable :: grad(:),posini(:,:),gmean(:),gmean_total(:)
real              :: rc,tbath,pext,press,ppot,rho,length,time,dt,nu,tcalc,eps,mass,sigma,length2
real              :: utime,utemp,upress,udens,epot,ekin,conteg,contep,meansq,interv,rho2,dgr,tmelt, p_mean
integer           :: npar, dim,ii,timesteps,outg,oute,equi,nbox
integer              :: taskid, numproc, ierror, partition1,partition2
integer, allocatable :: pairindex(:,:)
integer, allocatable :: table_index1(:,:),table_index2(:,:)
integer           :: i,j
integer, allocatable :: seed(:)
integer :: nseed
integer, allocatable :: displs(:), counts(:), index_local(:,:)
integer              :: ini, ind_i, ind_j, partition_jon, max_length

! -------------------------
! INTIALIZE MPI ENVIRONMENT
! -------------------------
call MPI_INIT(ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD,taskid,ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)


open(1,file='input.txt',status='old')   !Input file

read(1,*) equi          ! equilibration time before starting measures
read(1,*) timesteps     ! integer timesteps
read(1,*) dim           ! spatial dimension
read(1,*) npar          ! number of particles
read(1,*) sigma         ! energy epsilon
read(1,*) eps           ! energy epsilon
read(1,*) mass          ! atom mass g/mol
read(1,*) dt            ! simulation timestep
read(1,*) nu            ! Andersen thermostat coupling parameter
read(1,*) rho           ! density
read(1,*) rc            ! cutoff radious, in reduce units, where 1d0 is the particle radius
read(1,*) tbath         ! thermal bath temperature, In kelvin
read(1,*) pext          ! system pressure, in atm
read(1,*) oute,outg     ! number of timesteps to measure g(r) and MSD
read(1,*) nbox          ! number of positions to calculate radial distribution function



call random_seed(size = nseed) !initialisation of a random seed for each processor
allocate(seed(nseed))
call random_seed(get=seed)


! Alocation of main variables:

    ! pos(npar,dim) --------> atoms position
    ! vel(npar,dim) --------> velocities at time t
    ! fpar(npar,dim) -------> calculation of forces acting on each particle
    ! posini(npar,dim) -----> fixed position to calculate mean square displacement
    ! pairindex((npar*(npar-1))/2,2) ----> vector of pairs of particles
allocate(pos(npar,dim),vel(npar,dim),f_par(npar,dim),posini(npar,dim))   !allocation of variables
allocate(grad(nbox),gmean(nbox),gmean_total(nbox))
allocate(pairindex((npar*(npar-1))/2,2))

!Define vector of pairs of particles
ii = 1
do i = 1,npar-1
    do j = i+1,npar
        pairindex(ii,:) = (/i,j/)
        ii = ii + 1
    enddo
enddo

! Files where data will be saved:
! Units output: time --> ps, energy --> J/mol, temperature --> Kelvin, pressure --> Atmospheres
if (taskid == 0) then
    open(10,file='data_EK_EP_T_P.dat',status='unknown')
    open(20,file='mean_square_disp.dat',status='unknown')
    open(30,file='rad_dist_func.dat',status='unknown')
endif


! Calculation of the magnitudes in reduced units
call reduced_units(eps,mass,sigma,utime,utemp,upress,udens)

rho   = rho   / udens
pext  = pext  / upress
tbath = tbath / utemp

tmelt = 100.0 * tbath

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! DESCOMENTAR ESTO PARA HACER EL CALCULO DEL ARGON QUE HACIAMOS EN MOMO
eps   = 998.0
sigma = 3.4
rho   = 0.8
mass  = 40.0
Tbath = 10.0
Pext  = 2.5
tmelt = 100.0 * tbath
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

!DIVISIÓN DE LA INFORMACIÓN POR WORKERS

! Tabla con los indices de los atomos que corresponden a cada worker
allocate(table_index1(0:numproc-1,2))
allocate(table_index2(0:numproc-1,2))
partition1 = npar/numproc
partition2 = (npar*(npar-1)/2)/numproc
do ii = 0, numproc-2
  table_index1(ii,1) = ii*partition1+1
  table_index2(ii,1) = ii*partition2+1
  table_index1(ii,2) = (ii+1)*partition1
  table_index2(ii,2) = (ii+1)*partition2
enddo
ii = numproc-1
table_index1(ii,1) = ii*partition1+1
table_index2(ii,1) = ii*partition2+1
table_index1(ii,2) = npar
table_index2(ii,2) = npar*(npar-1)/2


! ------------------------------------------------------------------------

! DIVISION TRABAJO: SUBRUTINA INTEGRACION
partition_jon  = 3*npar/numproc                    
max_length = partition_jon + mod(3*npar,numproc)   

allocate(displs(numproc), counts(numproc), index_local(numproc,2))

do i = 1, numproc
  displs(i) = (i-1)*partition_jon
  
  ini = displs(i) + 1

  ind_i = mod(ini, npar)
  ind_j = ini/npar + 1
  if (ind_i == 0) then
    ind_i = npar
    ind_j = ind_j - 1
  endif
  index_local(i,1) = ind_i
  index_local(i,2) = ind_j

enddo
counts          = partition_jon
counts(numproc) = max_length

! ------------------------------------------------------------------------



! Initialisation of the system velocity and position

call initialize_r(pos,rho,Npar,length)
call initialize_v(vel,Tmelt,Npar,numproc,taskid,ierror)

! Cutoff radio as a funtion of the length of the box
rc = 0.48*length


! Equilibration and melting of the system: We leave a certain number of timesteps to destroy initial order

time=0.0

do ii=1,equi

    call force(npar,length,rc,pos,vel,f_par,press,epot,ekin)
    call v_verlet_step (npar, pos, vel, time, dt, Rc, Length, F_par, press, Epot, Ekin, &
                        numproc, taskid, max_length, displs, counts, index_local, ierror,&
                        tbath,nu)
!  Everyone worked on their subrutine but we mixed with the v_verlet_step in order to avoid comunications
!    call PBC(npar,length,pos,numproc,taskid,counts,displs,&
!         max_length,ierror)
!    call thermostat (npar, vel, nu, tbath,numproc,taskid,max_length, &
!         displs, counts, index_local, ierror)

enddo

! Reinitialize velocities according to the bath temperature

call initialize_v(vel,Tmelt,Npar,numproc,taskid,ierror)


! Equilibrate with real bath temperature

time = 0.0

do ii = 1, equi
    call force(npar,length,rc,pos,vel,f_par,press,epot,ekin)
    call v_verlet_step (npar, pos, vel, time, dt, Rc, Length, F_par, press, Epot, Ekin, &
                        numproc, taskid, max_length, displs, counts, index_local, ierror, &
                        tbath,nu)
!    call PBC(npar,length,pos,numproc,taskid,counts,displs,&
!         max_length,ierror)
!    call thermostat (npar, vel, nu, tbath,numproc,taskid,max_length, &
!    displs, counts, index_local, ierror)
enddo


posini = pos
contep = 0.0
conteg = 0.0
time   = 0.0
gmean  = 0.0
p_mean = 0.0

! Loop of times where velocity-Verlet algorithm is called
! conteg : counts the times we call g(r) radial distribution function
! gmean  : accumulated distribution function

do ii = 1, timesteps
    
    call force(npar,length,rc,pos,vel,f_par,press,epot,ekin)
    call v_verlet_step (npar, pos, vel, time, dt, Rc, Length, F_par, press, Epot, Ekin, &
                        numproc, taskid, max_length, displs, counts, index_local, ierror,&
                        tbath,nu)
!    call PBC(npar,length,pos,numproc,taskid,counts,displs,&
!         max_length,ierror)
!    call thermostat (npar, vel, nu, tbath,numproc,taskid,max_length, &
!    displs, counts, index_local, ierror)
    
    if (mod(ii,oute).eq.0) then
        write(10,'(5f20.8,i12)')  time*utime, ekin*eps, epot*eps, tcalc*utemp, press*upress, ii
    endif

    tcalc = 2.0 * Ekin / (3.0*npar)
    p_mean = p_mean + press        

    if (mod(ii,outg).eq.0) then
        
        conteg = conteg + 1.0                  
        
        call msdisplacement(numproc,taskid,table_index1,npar,dim,posini,pos,length,meansq)
        call gr(numproc,taskid,table_index2,pairindex,npar,dim,rho,length,pos,1.0*length,nbox,grad,dgr)
        if (taskid == 0) then
            write(20,*) time * utime, meansq * sigma**2.0
            print*,ii
        endif

        gmean = gmean + grad      
    
    endif

enddo

call MPI_REDUCE(gmean,gmean_total,size(gmean),MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierror)


!End of simulation: Now we write the radial distribution function
if (taskid == 0) then
do ii = 1, nbox

    write(30,*) sigma * dgr * real(ii), gmean_total(ii) / conteg

enddo
endif

call MPI_FINALIZE(ierror)

end


!-----------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------
!------------------------------------------SUBROUTINE REDUCED UNITS-----------------------------------------------
!-----------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------
!----------------------------Calculates all required unit magnitudes of the simulation----------------------------

subroutine reduced_units(uener,mass,sigma,utime,utemp,upress,udens)

implicit none
real :: mass,sigma,utime,utemp,uener,upress,udens

upress = (uener/(sigma**3.0))*16.3886  !pressure in atm (J/mol)/A^3=10^30/6.022E23 Pa * 1atm/101325Pa
utemp  = uener / 8.31                               ! temperature in Kelvin units (uener-->J/mol, 8.31---> J/(K·mol)
udens  = (1.0 / 0.6022) * (sigma**3.0) / (mass)     ! 0.6022: factor (10^-8)^3*(6.022^10^(-23))
utime  = sigma * sqrt((mass*1.0E-3)/uener) * 100.0  ! unit time in picoseconds

return

end
