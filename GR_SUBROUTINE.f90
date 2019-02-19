!------------------------------------------------------------------------
subroutine gr(N,L,COORD,g,density)
integer:: N
real(8):: density
real(8):: L
real(8):: COORD(3,N)
real(8):: xmin, xmax,interv,contint
integer:: Histogram(ncajas),caja
real(8):: g(ncajas)
integer:: Nx,Ny,Nz
integer:: i,j
real(8):: r(3),rL(3)
real(8):: R2,R2_caja
integer:: cont


xmin = 0.0d0
xmax = 2.0D0*L
interv = (xmax-xmin)/dble(ncajas)

g = 0.0d0
Histogram = 0
cont = 0

do i = 1,N
	do j = i+1,N
		r = (COORD(:,j)-COORD(:,i))
		do Nx = -1,1
			do Ny = -1,1
				do Nz = -1,1
					rL = r + L*(/Nx,Ny,Nz/)
					R2 = sqrt(sum(rL**2))
					caja = int(r2/interv) + 1
					if (caja<0) then
						cycle
					elseif(caja>ncajas) then
						cycle
					else
					Histogram(caja) = Histogram(caja) + 1
					cont = cont + 1
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