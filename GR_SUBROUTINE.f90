!------------------------------------------------------------------------
subroutine gr(density,L,R,ncajas,g)
real(8),parameter:: pi = acos(-1.0d0)
real(8),parameter:: con = (4.0d0/3.0d0)*pi
real(8),intent(in):: density,L
real(8),intent(in):: R(:,:)
integer,intent(in):: ncajas
integer:: N
real(8):: xmin, xmax,interv,contint
integer:: Histogram(ncajas),caja
real(8),intent(out):: g(ncajas)
integer:: Nx,Ny,Nz
integer:: i,j
real(8):: rad(3),rL(3)
real(8):: R2

N = size(R,1)

xmin = 0.0d0
xmax = 2.0D0*L
interv = (xmax-xmin)/dble(ncajas)

g = 0.0d0
Histogram = 0

do i = 1,N-1
	do j = i+1,N
		rad = (R(j,:)-R(i,:))
		do Nx = -1,1
			do Ny = -1,1
				do Nz = -1,1
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


contint = xmin
do j = 1,ncajas
	if (Histogram(j)/=0) then
		g(j) = g(j) + dble(Histogram(j))/&
		((con*((contint+interv)**3-(contint)**3))*density)
	endif
	contint = contint + interv
enddo

return
end subroutine gr
!------------------------------------------------------------------------