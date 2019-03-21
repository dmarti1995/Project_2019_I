subroutine dist(N,dim,L,COORD,R2ij)
!N(in); numero de part
!dim(in): dimension sistema
!L(in): system length
!COORD(in): vector of the positions of the system
!R2ij(out): Array of distances between pairs of particles i,j
integer,intent(in):: N,dim
real,intent(in):: L
real,intent(in):: COORD(N,dim)
real,intent(out):: R2ij(N,N)
real:: Rij(dim)
integer:: i,j


R2ij = 0.0

do i = 1,N-1
	do j = i+1,N
		Rij = (COORD(j,:)-COORD(i,:)) 
		Rij = Rij -nint(Rij/L)*L
		R2ij(i,j) = (sum(Rij**2))
		R2ij(j,i) = R2ij(i,j)
	enddo
enddo

return
end subroutine dist