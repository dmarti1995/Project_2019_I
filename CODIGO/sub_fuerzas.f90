!***************************************************************
!In this subroutine we compute Forces,Pressure,Potential and kinetic energy.
!***************************************************************
!Input variables:
!N: number of particles
!L: box length
!cutoff: maximum distance for force computation
!x: particles positions
!v: particles velocity
!***************************************************************
!Output variables:
!F: partciles forces
!press: pressure
!epot: potential energy
!ekin: kinetic energy
!***************************************************************
!***************************************************************

subroutine force(dim,N,L,cutoff,x, v, F,press,epot,ekin,celllist,ipc,ncells)
implicit none
integer :: dim
integer :: ipc(N,dim+1),ncells(3)
integer :: ipar,j,N,icell,jcell,kcell,icelltmp,jcelltmp,kcelltmp,kpar,jpar
real :: L,dx,dy,dz,rij2,rij8,rij14,cutoff,forcedist,press,potdist,rij12,rij6,epot,ekin
integer :: celllist(0:10,ncells(1),ncells(2),ncells(3))
real, dimension(N,3) :: x,F,v
real :: shift,shift14,shift8,shift6,shift12,shiftenergy
F=0.0
press=0d0
epot=0d0

!do i=1,N main do 1
!  do j=1,N  main do 2
do ipar=1,N          !MAIN DO!

    do icelltmp=ipc(ipar,1)-1,ipc(ipar,1)+1  !i loop  we iterate for all the neighbour cells to that
    do jcelltmp=ipc(ipar,2)-1,ipc(ipar,2)+1  !j loop  of i particle (ipar)
    do kcelltmp=ipc(ipar,3)-1,ipc(ipar,3)+1  !k loop

        icell=icelltmp-ncells(1)*floor(real(icelltmp-1)/ncells(1))  !apply boundary conditions
        jcell=jcelltmp-ncells(2)*floor(real(jcelltmp-1)/ncells(2))  !to all the cells
        kcell=kcelltmp-ncells(3)*floor(real(kcelltmp-1)/ncells(3))


        do kpar=1,celllist(0,icell,jcell,kcell)   !for each cell, we iterate for the particles
                                                  !celllist(0,icell,jcell,kcell  contains
                                                  !number of particles in cell i,j,k
            jpar=celllist(kpar,icell,jcell,kcell)   !celllist(kpar,icell,jcell,kcell) contains
            if (ipar.eq.jpar) cycle                 !the kpar index of particle i,j,k
            dx=x(ipar,1)-x(jpar,1)                  !kpar wont be greater than 9, so that's why
                                                    !it goes form 1 to 9.
            dx=dx-nint(dx/L)*L !PBC for distances
            dy=x(ipar,2)-x(jpar,2)
            dy=dy-nint(dy/L)*L
            dz=x(ipar,3)-x(jpar,3)
            dz=dz-nint(dz/L)*L
            rij2=dx**2.0D0+dy**2.0D0+dz**2.0D0
            if (rij2<cutoff**2.0) then !Computation only for shorter than cutoff
                rij8=rij2**4.0
                rij14=rij2**7.0
                rij6=rij2**3.0
                rij12=rij2**6
                shift14=1d0/(cutoff**14.0)
                shift8=1d0/(cutoff**8.0)
                shift=(48.0D0*shift14 - 24.0D0*shift8)
                forcedist=(48.0D0/rij14 - 24.0D0/rij8)-shift  !shift force --> force has to be 0 at rc
                shift6=1/(cutoff**6)
                shift12=1/(cutoff**12)
                shiftenergy=2d0*((1d0*shift12)-(1d0*shift6))
              !  print*,shift
                potdist=2d0*((1d0/rij12)-(1d0/rij6))-shiftenergy !shift energy --> energy has to be 0 at rc
                F(ipar,1)=F(ipar,1)+forcedist*dx     !force matrix
                F(ipar,2)=F(ipar,2)+forcedist*dy
                F(ipar,3)=F(ipar,3)+forcedist*dz
                epot=epot+potdist              !potential energy calculation
            endif
        enddo
   enddo                                !i loop
   enddo
   enddo


  ! Potential pressure
  do jpar = 1, 3
    press = press +x(ipar,jpar)*F(ipar,jpar)
  enddo

enddo    !MAIN DO!

! kinetic energy
ekin = 0.5*sum(v**2)

! total pressure
press = press + 2.0*ekin
press = (1.0/3.0)*(1.0/L**3)*press


end subroutine
