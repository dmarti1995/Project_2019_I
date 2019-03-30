
!------------------------------------------------------------------------
!Rcurrent(in): position of the particles at the current time step
!Rfixed(in): fixed configuration of the particles from which the meansquaredisp is computed
!L(in): system length
!meansquaredisp(out): Value of the meansquaredisp of this time step
subroutine MSDISPLACEMENT(numproc,taskid,table_index1,npar,dim,Rcurrent,Rfixed,L,meansquaredisp)
  USE MPI
implicit none
integer, parameter   :: dp = 8
integer,intent(in)   :: numproc,taskid,table_index1(0:numproc-1,2)
integer,intent(in)   :: npar,dim
real(dp),intent(in)  :: Rcurrent(npar,dim),Rfixed(npar,dim)
real(dp),intent(in)  :: L
real(dp),intent(out) :: meansquaredisp
real(dp)             :: R(3)
real(dp)             :: suma_local,suma_total
integer              :: i
integer              :: imin,imax
integer              :: ierror

imin = table_index1(taskid,1)
imax = table_index1(taskid,2)

suma_local = 0.0_dp
do i = imin,imax
	R = Rcurrent(i,:)-Rfixed(i,:)
	R = R - L*nint(R/L)
	suma_local = suma_local + (sum(R**2.0_dp))
enddo


call MPI_REDUCE(suma_local,suma_total,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierror)

if (taskid==0) then
	meansquaredisp = suma_total/dble(npar)
else
	meansquaredisp = 0.0_dp
endif
return
end subroutine MSDISPLACEMENT
!------------------------------------------------------------------------
