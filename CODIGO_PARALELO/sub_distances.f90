subroutine dist(N,dim,L,COORD,R2ij,MATRR2ij)
!N(in); numero de part
!dim(in): dimension sistema
!L(in): system length
integer,intent(in):: N,dim
real,intent(in):: L
real,intent(in):: COORD(N,dim)
real,intent(out):: R2ij(N*(N-1)/2),MATRR2ij(i,j)
real:: Rij(dim)
integer:: i,j,cont


R2ij = 0.0
MATRR2ij = 0.0

cont = 1

do i = 1,N-1
	do j = i+1,N
		Rij = (COORD(j,:)-COORD(i,:)) 
		Rij = Rij - nint(Rij/L)*L
		R2ij(cont) = (sum(Rij**2))
		MATRR2ij(i,j) = R2ij(cont)
		MATRR2ij(j,i) = R2ij(cont)
		cont = cont + 1
	enddo
enddo

return
end subroutine dist