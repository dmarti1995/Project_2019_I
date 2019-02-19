subroutine force()
implicit none
integer :: i,j,N
real :: L,dx,dy,dz,rij
real :: 
do i=1,N
  do j=1,N
    if (i/=j) then
      dx=x(i,1)-x(j,1)
      dx=dx-nint(dx/L)*L
      dy=x(i,2)-x(j,2)
      dy=dy-nint(dy/L)*L
      dz=x(i,3)-x(j,3)
      dz=dz-nint(dz/L)*L
      rij=sqrt(dx**2.0D0+dy**2.0D0+dz**2.0D0)
      if (rij<cutoff) then
        F(i,1)=F(i,1)+(48.0D0/rij**14.0D0 - 24.0D0/rij**8.0D0)*dx 
        F(i,2)=F(i,2)+(48.0D0/rij**14.0D0 - 24.0D0/rij**8.0D0)*dy 
        F(i,3)=F(i,3)+(48.0D0/rij**14.0D0 - 24.0D0/rij**8.0D0)*dz 
      endif
    endif
  enddo
enddo

