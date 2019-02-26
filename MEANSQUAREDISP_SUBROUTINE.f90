!------------------------------------------------------------------------
subroutine MSDISPLACEMENT(Rcurrent,Rfixed,L,meansquaredisp)
real(8),intent(in):: Rcurrent(:,:),Rfixed(:,:)
real(8),intent(in):: L
real(8),intent(out):: meansquaredisp
real(8):: R(3)
real(8):: suma
integer:: i


suma = 0.0D0
do i = 1,size(Rcurrent,1)
	R = Rcurrent(i,:)-Rfixed(i,:)
	call PBC(1,L,R) 
	suma = suma + (sum(R**2))
enddo

meansquaredisp = suma/dble(size(Rcurrent,1))
end subroutine MSDISPLACEMENT
!------------------------------------------------------------------------