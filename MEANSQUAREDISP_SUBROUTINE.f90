!------------------------------------------------------------------------
subroutine MSDISPLACEMENT(N,RN,Rfixed,meansquaredisp,L)
integer:: N
real(8):: RN(3,N),Rfixed(3,N)
real(8):: rj_imatge(3)
real(8):: suma,meansquaredisp
real(8):: L
integer:: i
real(8):: D2

suma = 0.0D0
do i = 1,N
	call PBCD(RN(:,i),Rfixed(:,i),L,D2,rj_imatge)
	suma = suma + D2
enddo

meansquaredisp = suma/dble(N)
end subroutine MSDISPLACEMENT
!------------------------------------------------------------------------