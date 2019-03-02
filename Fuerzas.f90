subroutine force(N,L,cutoff,x, v, F,press,epot,ekin)
implicit none
integer :: i,j,N
real :: L,dx,dy,dz,rij2,rij8,rij14,cutoff,forcedist,press,potdist,rij12,rij6,epot,ekin
real, dimension(N,3) :: x,F,v
F=0.0
press=0d0
epot=0d0

do i=1,N
  do j=1,N
    if (i/=j) then
      dx=x(i,1)-x(j,1)
      dx=dx-nint(dx/L)*L
      dy=x(i,2)-x(j,2)
      dy=dy-nint(dy/L)*L
      dz=x(i,3)-x(j,3)
      dz=dz-nint(dz/L)*L
      rij2=dx**2.0D0+dy**2.0D0+dz**2.0D0
      if (rij2<cutoff**2.0) then
        rij8=rij2**4.0
        rij14=rij2**7.0
        rij6=rij2**3.0
        rij12=rij2**6
        forcedist=(48.0D0/rij14 - 24.0D0/rij8)
        potdist=2d0*((1d0/rij12)-(1d0/rij6))
        F(i,1)=F(i,1)+forcedist*dx
        F(i,2)=F(i,2)+forcedist*dy
        F(i,3)=F(i,3)+forcedist*dz
        epot=epot+potdist              !potential energy calculation
      endif
    endif
  enddo
  ! Potential pressure
  do j = 1, 3
    press = press +x(i,j)*F(i,j)
  enddo
enddo

! kinetic energy
ekin = 0.5*sum(v**2)

! total pressure
press = press + 2.0*ekin
press = (1.0/3.0)*(1.0/L**3)*press


end subroutine
