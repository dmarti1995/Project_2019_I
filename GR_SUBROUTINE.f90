!------------------------------------------------------------------------
subroutine gr(N,dim,density,L,COORD,rmax,ncajas,g,dr)
!N(in); numero de part
!dim(in): dimension sistema
!density(in): system density
!L(in): system length
!COORD(in): vector of the positions of the system
!rmax(in): radio máximo hasta el que se calcula g(r)
!ncajas(in): number of r points in which the g(r) is computed
!g(out): vector g(ncajas) of this time step
!dr(out): valor del intervalo dr
implicit none
real,parameter:: pi = acos(-1.0)
real,parameter:: con = (4.0/3.0)*pi
integer,intent(in):: N,dim
real,intent(in):: density,L
real,intent(in):: COORD(N,dim)
real,intent(in):: rmax 
integer,intent(in):: ncajas
real:: contint
integer:: Histogram(ncajas),caja
real,intent(out):: g(ncajas)
real:: dr
integer:: Nx,Ny,Nz
integer:: i,j
real:: r(3),rL(3)
real:: R2,R2_caja


g = 0.0d0
dr = rmax/dble(ncajas)


Histogram = 0

!All possible particle pairs
do i = 1,N-1
	do j = i+1,N
		r = (COORD(j,:)-COORD(i,:))
!We study our box and the nearest neighbour boxes
		do Nx = -1,1
			do Ny = -1,1
				do Nz = -1,1
!we make an histogram of all the images of the particle pairs.
					rL = r + L*(/Nx,Ny,Nz/)
					R2 = sqrt(sum(rL**2))
					caja = int(r2/dr) + 1
					if (caja<0) then
						cycle
					elseif(caja>ncajas) then
						cycle
					else
					Histogram(caja) = Histogram(caja) + 1
					endif
				enddo
			enddo
		enddo
	enddo
enddo

!From the Histogram we compute the g(r) of this time step 
contint = 0.0d0
do j = 1,ncajas
	if (Histogram(j)/=0) then
		g(j) =  dble(Histogram(j))/&
		(2*N*(con*((contint+dr)**3-(contint)**3))*density)
	endif
	contint = contint + dr
enddo

return
end subroutine gr
!------------------------------------------------------------------------
