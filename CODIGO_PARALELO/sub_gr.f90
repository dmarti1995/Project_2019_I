!------------------------------------------------------------------------
subroutine gr(numproc,taskid,table_index2,pairindex,N,dim,density,L,COORD,rmax,ncajas,g,dr)
!N(in); numero de part
!dim(in): dimension sistema
!density(in): system density
!L(in): system length
!COORD(in): vector of the positions of the system
!rmax(in): radio máximo hasta el que se calcula g(r)
!ncajas(in): number of r points in which the g(r) is computed
!g(out): vector g(ncajas) of this time step
!dr(out): valor del intervalo dr
  USE MPI
implicit none
integer, parameter   :: dp = 8
real(dp),parameter   :: pi = acos(-1.0_dp)
real(dp),parameter   :: con = (4.0_dp/3.0_dp)*pi
integer,intent(in)   :: numproc,taskid,table_index2(0:numproc-1,2),pairindex((N*(N-1))/2,2)
integer,intent(in)   :: N,dim
real(dp),intent(in)  :: density,L
real(dp),intent(in)  :: COORD(N,dim)
real(dp),intent(in)  :: rmax 
integer,intent(in)   :: ncajas
real(dp)             :: contint
integer              :: Histogram(ncajas),caja
real(dp),intent(out) :: g(ncajas)
real(dp)             :: dr
integer              :: Nx,Ny,Nz
integer              :: i,j
real(dp)             :: r(3),rL(3)
real(dp)             :: R2,R2_caja
integer              :: imin,imax,cont

imin = table_index2(taskid,1)
imax = table_index2(taskid,2)

g = 0.0_dp
dr = rmax/dble(ncajas)


Histogram = 0

!Study the selected particles

do cont= imin,imax
i = pairindex(cont,1)
j = pairindex(cont,2)
!print*,i,j,cont
		r = (COORD(j,:)-COORD(i,:))
!We study our box and the nearest neighbour boxes
	do Nx = -1,1
		do Ny = -1,1
			do Nz = -1,1
!we make an histogram of all the images of the particle pairs.
				rL = r + L*(/Nx,Ny,Nz/)
				R2 = sqrt(sum(rL**2.0_dp))
				caja = int(r2/dr) + 1
				if (caja<0) then
					cycle
				elseif(caja>ncajas) then
					cycle
				else
				Histogram(caja) = Histogram(caja) + 2
				endif
			enddo
		enddo
	enddo
enddo

!From the Histogram we compute the g(r) of this time step 
contint = 0.0_dp
do j = 1,ncajas
	if (Histogram(j)/=0) then
		g(j) =  dble(Histogram(j))/&
		(N*(con*((contint+dr)**3.0_dp-(contint)**3.0_dp))*density)
	endif
	contint = contint + dr
enddo

return
end subroutine gr
!------------------------------------------------------------------------


!------------------------------------------------------------------------
subroutine gr2(numproc,taskid,table_index2,pairindex,N,dim,L,COORD,rmax,ncajas,Histogram,dr)
!N(in); numero de part
!dim(in): dimension sistema
!density(in): system density
!L(in): system length
!COORD(in): vector of the positions of the system
!rmax(in): radio máximo hasta el que se calcula g(r)
!ncajas(in): number of r points in which the g(r) is computed
!g(out): vector g(ncajas) of this time step
!dr(out): valor del intervalo dr
  USE MPI
implicit none
integer, parameter  :: dp = 8
integer,intent(in)  :: numproc,taskid,table_index2(0:numproc-1,2),pairindex((N*(N-1))/2,2)
integer,intent(in)  :: N,dim
real(dp),intent(in) :: L
real(dp),intent(in) :: COORD(N,dim)
real(dp),intent(in) :: rmax 
integer,intent(in)  :: ncajas
integer,intent(out) :: Histogram(ncajas)
integer             :: caja
real(dp)            :: dr
integer             :: Nx,Ny,Nz
integer             :: i,j
real(dp)            :: r(3),rL(3)
real(dp)            :: R2,R2_caja
integer             :: imin,imax,cont

imin = table_index2(taskid,1)
imax = table_index2(taskid,2)


dr = rmax/dble(ncajas)


Histogram = 0

!Study the selected particles

do cont= imin,imax
i = pairindex(cont,1)
j = pairindex(cont,2)
!print*,i,j,cont
		r = (COORD(j,:)-COORD(i,:))
!We study our box and the nearest neighbour boxes
	do Nx = -1,1
		do Ny = -1,1
			do Nz = -1,1
!we make an histogram of all the images of the particle pairs.
				rL = r + L*(/Nx,Ny,Nz/)
				R2 = sqrt(sum(rL**2.0_dp))
				caja = int(r2/dr) + 1
				if (caja<0) then
					cycle
				elseif(caja>ncajas) then
					cycle
				else
				Histogram(caja) = Histogram(caja) + 2
				endif
			enddo
		enddo
	enddo
enddo

return
end subroutine gr2
!------------------------------------------------------------------------
