!------------------------------------------------------------------------
!density(in): system density
!L(in): system length
!R(in): vector of the positions of the system
!ncajas(in): number of r points in which the g(r) is computed
!g(out): vector g(ncajas) of this time step
subroutine gr(density,L,R,ncajas,g,interv)
real,parameter:: pi = acos(-1.0)
real,parameter:: con = (4.0/3.0)*pi
real,intent(in):: density,L
real,intent(in):: R(:,:)
integer,intent(in):: ncajas
integer:: N
real:: xmin, xmax,contint
integer:: Histogram(ncajas),caja
real,intent(out):: g(ncajas)
real,intent(out):: interv
integer:: Nx,Ny,Nz
integer:: i,j
real:: rad(3),rL(3)
real:: R2

N = size(R,1)

xmin = 0.0
xmax = 2.0*L
interv = (xmax-xmin)/dble(ncajas)


Histogram = 0

!All possible particle pairs
do i = 1,N-1
	do j = i+1,N
!Distance between them
		rad = (R(j,:)-R(i,:))
!We study our box and the nearest neighbour boxes
		do Nx = -1,1
			do Ny = -1,1
				do Nz = -1,1
!we make an histogram of all the images of the particle pairs.
					rL = rad + L*(/Nx,Ny,Nz/)
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