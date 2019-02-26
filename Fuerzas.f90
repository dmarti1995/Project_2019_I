subroutine force(N,L,cutoff,x,F,press,epot)
implicit none
integer :: i,j,N
real :: L,dx,dy,dz,rij2,rij8,rij14,cutoff,forcedist,press,potdist,rij12,rij6,epot
real, dimension(N,3) :: x,F
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
        potdist=4d0*((1d0/rij12)-(1d0/rij6))
        F(i,1)=F(i,1)+forcedist*dx
        F(i,2)=F(i,2)+forcedist*dy
        F(i,3)=F(i,3)+forcedist*dz
        press=press+rij2*forcedist     !pressure calculation
        epot=epot+potdist              !potential energy calculation
      endif
    endif
  enddo
enddo
end subroutine
