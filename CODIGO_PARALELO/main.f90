program DinMo

implicit none

include 'mpif.h'

real, allocatable :: pos(:,:),vel(:,:),f_par(:,:)
real, allocatable :: grad(:),posini(:,:),gmean(:)
real              :: rc,tbath,pext,press,ppot,rho,length,time,dt,nu,tcalc,eps,mass,sigma,length2
real              :: utime,utemp,upress,udens,epot,ekin,conteg,contep,meansq,interv,rho2,dgr,tmelt, p_mean
integer           :: npar, dim,ii,timesteps,outg,oute,equi,nbox
integer              :: taskid, numproc, ierror, partition
integer, allocatable :: table_index(:,:)


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


! Alocation of main variables:

    ! pos(npar,dim) --------> atoms position
    ! vel(npar,dim) --------> velocities at time t
    ! fpar(npar,dim) -------> calculation of forces acting on each particle
    ! posini(npar,dim) -----> fixed position to calculate mean square displacement

allocate(pos(npar,dim),vel(npar,dim),f_par(npar,dim),posini(npar,dim))   !allocation of variables
allocate(grad(nbox),gmean(nbox))


! Files where data will be saved:
! Units output: time --> ps, energy --> J/mol, temperature --> Kelvin, pressure --> Atmospheres

open(10,file='data_EK_EP_T_P.dat',status='unknown')
open(20,file='mean_square_disp.dat',status='unknown')
open(30,file='rad_dist_func.dat',status='unknown')


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
!eps   = 998.0
!sigma = 3.4
!rho   = 0.8
!mass  = 40.0
!Tbath = 10.0
!Pext  = 2.5
!tmelt = 100.0 * tbath
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------


! -------------------------
! INTIALIZE MPI ENVIRONMENT
! -------------------------
call MPI_INIT(ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD,taskid,ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)

! Tabla con los indices de los atomos que corresponden a cada worker
allocate(table_index(0:numproc-1,2))
partition = npar/numproc
do ii = 0, numproc-2
  table_index(ii,1) = ii*partition+1
  table_index(ii,2) = (ii+1)*partition
enddo
ii = numproc-1
table_index(ii,1) = ii*partition+1
table_index(ii,2) = npar


! Initialisation of the system velocity and position

call initialize_r(pos,rho,Npar,length)
call initialize_v(vel,Tmelt,Npar)

! Cutoff radio as a funtion of the length of the box
rc = 0.48*length


! Equilibration and melting of the system: We leave a certain number of timesteps to destroy initial order

time=0.0

do ii=1,equi

    call force(npar,length,rc,pos,vel,f_par,press,epot,ekin)
    call v_verlet_step (npar, pos, vel, time, dt, rc, length, f_par,press,epot,ekin)
    call PBC(npar, length, pos)
    call thermostat (npar, vel, nu, Tmelt)

enddo

! Reinitialize velocities according to the bath temperature

call initialize_v(vel,tbath,Npar)

! Equilibrate with real bath temperature

time = 0.0

do ii = 1, equi
    call force(npar,length,rc,pos,vel,f_par,press,epot,ekin)
    call v_verlet_step (npar, pos, vel, time, dt, rc, length, f_par,press,epot,ekin)
    call PBC(npar, length, pos)
    call thermostat (npar, vel, nu, Tbath)
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
    call v_verlet_step (npar, pos, vel, time, dt, rc, length, f_par,press,epot,ekin)
    call PBC(npar, length, pos)
    call thermostat (npar, vel, nu, Tbath)
    
    if (mod(ii,oute).eq.0) then
        write(10,'(5f20.8,i12)')  time*utime, ekin*eps, epot*eps, tcalc*utemp, press*upress, ii
    endif

    tcalc = 2.0 * Ekin / (3.0*npar)
    p_mean = p_mean + press        

    if (mod(ii,outg).eq.0) then
        
        conteg = conteg + 1.0                  
        
        call msdisplacement(npar,dim,posini,pos,length,meansq)
        call gr(npar,dim,rho,length,pos,1.0*length,nbox,grad,dgr)
        write(20,*) time * utime, meansq * sigma**2.0
      
        gmean(:) = gmean(:) + grad(:)       
    
    endif

enddo

call MPI_FINALIZE(ierror)

! End of simulation: Now we write the radial distribution function

do ii = 1, nbox

    write(30,*) sigma * dgr * real(ii), gmean(ii) / conteg

enddo

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
