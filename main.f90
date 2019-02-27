program DinMo
implicit none
real, allocatable :: pos(:,:),vel(:,:),f_par(:,:),velinf(:,:),velsup(:,:)
real, allocatable :: grad(:),posini(:,:),gmean(:)
real :: rc,tbath,pext,press,ppot,rho,length,time,dt,taup,taut,tcalc,eps,mass,sigma,length2
real :: uener,utime,utemp,upress,udens,epot,ekin,conteg,contep,meansq,interv,rho2
integer :: npar, dim,ii,timesteps,outg,oute,equi,nbox


open(1,file='input.dat',status='old')   !Input file

read(1,*) equi          !equilibration time before starting measures
read(1,*) timesteps     !integer timesteps
read(1,*) dim        !spatial dimension
read(1,*) npar       !number of particles
read(1,*) sigma         !energy epsilon
read(1,*) eps         !energy epsilon
read(1,*) mass        !atom mass g/mol
read(1,*) dt         !simulation timestep
read(1,*) taut        !Berensden coupling parameter
read(1,*) taup        !Berensden coupling parameter
read(1,*) rho        !density
read(1,*) rc         !cutoff radious, in reduce units, where 1d0 is the particle radius
read(1,*) tbath      !thermal bath temperature, In kelvin
read(1,*) pext       !system pressure, in atm
read(1,*) oute,outg   !number of timesteps to measure g(r) and MSD
read(1,*) nbox        !Number of positions to calculate radial distribution function

!Alocation of main variables:
    !pos(npar,dim) ------> atoms position
    !vel(npar,dim) --------> velocities at time t
    !velinf(npar,dim)--------->velocities at time t-deltaT/2
    !velsup(npar,dim) ------->velocities at time t+deltaT/2
    !fpar(npar,dim) ------> calculation of forces acting on each particle
    !posini(npar,dim) ------> fixed position to calculate mean square displacement
allocate(pos(npar,dim),vel(npar,dim),velinf(npar,dim),velsup(npar,dim),f_par(npar,dim),posini(npar,dim))   !allocation of variables
allocate(grad(nbox),gmean(nbox))

!Files where data will be saved:
!units output: time ---> ps, energy --> KJ/mol, temperature --> kelvin, pressure ---> Atmospheres
open(10,file='data_EK_EP_T_P.dat',status='unknown')
open(20,file='mean_square_disp.dat',status='unknown')


!Calculation of the magnitudes in reduced units---
call reduced_units(eps,mass,sigma,utime,utemp,upress,udens)
rho=rho/udens
pext=pext/upress
tbath=tbath/utemp

!Initialisation of the system velocity and position
call initialize(pos,velinf,rho,tbath,npar,length)
velsup=velinf

time=0d0
!Equilibration of the system: We leave a certain number of timesteps to destroy initial order
do ii=1,equi
    call force(npar,length,rc,pos,f_par,ppot,epot)
    call leap_frog(npar,dim,pos,vel,velinf,velsup,time,f_par,ppot,dt,taut,taup,rc,length,tbath,pext,press,tcalc,ekin)
enddo

posini=pos
contep=0d0
conteg=0d0
time=0d0
!Bucle of times where leap-frog algorithm is called
do ii=1,timesteps
    call force(npar,length,rc,pos,f_par,ppot,epot)
    call leap_frog(npar,dim,pos,vel,velinf,velsup,time,f_par,ppot,dt,taut,taup,rc,length,tbath,pext,press,tcalc,ekin)
    if (mod(ii,oute).eq.0) then
        write(10,*) time*utime,ekin*eps*1d-3,epot*eps*1d-3,tcalc*utemp,press*upress
       ! print*, time*utime,ekin*eps*1d-3,epot*eps*1d-3,tcalc*utemp,press*upress,nbox
    endif
    if (mod(ii,outg).eq.0) then
    conteg=conteg+1d0
   ! call MSDISPLACEMENT(npar,dim,posini,pos,length,meansq)
  !  call gr(npar,dim,rho,length,pos,nbox,grad,interv)
  !  write(20,*) time*utime,meansq*sigma*sigma
   ! print*, grad(100),nbox,length,press*upres
 print*, time*utime,ekin*eps*1d-3,epot*eps*1d-3,tcalc*utemp,press*upress,length
    endif
enddo

!-------End of simulation: Now we write the radial distribution function
do ii=1,timesteps


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

upress=(uener/(sigma**3.0))*(1d0/1.013d5)     !pressure units in atmosphere
utemp=uener/8.31d0        !temperature in Kelvin units
udens=(1d0/0.6022)*(sigma**3)/(mass)   !0.6022: factor (10^-8)^3*(6.022^10^(-23))
utime=sigma*sqrt((mass*1d-3)/uener)*100d0 !unit time in picoseconds
return
end
