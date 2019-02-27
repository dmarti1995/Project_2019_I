!------------------------------------------------------------------------
!Rcurrent(in): position of the particles at the current time step
!Rfixed(in): fixed configuration of the particles from which the meansquaredisp is computed
!L(in): system length
!meansquaredisp(out): Value of the meansquaredisp of this time step
subroutine MSDISPLACEMENT(npar,dim,Rcurrent,Rfixed,L,meansquaredisp)
implicit none
integer,intent(in) :: npar,dim
real,intent(in):: Rcurrent(npar,dim),Rfixed(npar,dim)
real,intent(in):: L
real,intent(out):: meansquaredisp
real:: R(3)
real:: suma
integer:: i


suma = 0.0
do i = 1,size(Rcurrent,1)
	R = Rcurrent(i,:)-Rfixed(i,:)
	call PBC(1,L,R) 
	suma = suma + (sum(R**2))
enddo

meansquaredisp = suma/dble(size(Rcurrent,1))
end subroutine MSDISPLACEMENT
!------------------------------------------------------------------------
