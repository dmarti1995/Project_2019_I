subroutine force(N,L,cutoff,x,F)
implicit none
integer :: i,j,N
real :: L,dx,dy,dz,rij2,rij8,rij14,cutoff
real, dimension(N,3) :: x,F
F=0.0
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
        F(i,1)=F(i,1)+(48.0D0/rij14 - 24.0D0/rij8)*dx 
        F(i,2)=F(i,2)+(48.0D0/rij14 - 24.0D0/rij8)*dy 
        F(i,3)=F(i,3)+(48.0D0/rij14 - 24.0D0/rij8)*dz 
      endif
    endif
  enddo
enddo
end subroutine
