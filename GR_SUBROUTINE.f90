!------------------------------------------------------------------------
subroutine gr(N,dim,density,L,COORD,ncajas,g,interv)
!N(in); numero de part
!dim(in): dimension sistema
!density(in): system density
!L(in): system length
!COORD(in): vector of the positions of the system
!ncajas(in): number of r points in which the g(r) is computed
!g(out): vector g(ncajas) of this time step
!interv(out): valor del intervalo dr
implicit none
real,parameter:: pi = acos(-1.0)
real,parameter:: con = (4.0/3.0)*pi
integer,intent(in):: N,dim
real(8),intent(in):: density,L
real(8),intent(in):: COORD(N,dim)
integer,intent(in):: ncajas
real(8):: xmin, xmax,contint
integer:: Histogram(ncajas),caja
real(8),intent(out):: g(ncajas),interv
integer:: Nx,Ny,Nz
integer:: i,j
real(8):: r(3),rL(3)
real(8):: R2,R2_caja



xmin = 0.0d0
xmax = 2.0D0*L
interv = (xmax-xmin)/dble(ncajas)


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
					caja = int(r2/interv) + 1
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
contint = xmin
do j = 1,ncajas
	if (Histogram(j)/=0) then
		g(j) =  dble(Histogram(j))/&
		((con*((contint+interv)**3-(contint)**3))*density)
	endif
	contint = contint + interv
enddo

return
end subroutine gr
!------------------------------------------------------------------------