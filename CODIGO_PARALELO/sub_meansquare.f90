
!------------------------------------------------------------------------
!Rcurrent(in): position of the particles at the current time step
!Rfixed(in): fixed configuration of the particles from which the meansquaredisp is computed
!L(in): system length
!meansquaredisp(out): Value of the meansquaredisp of this time step
subroutine MSDISPLACEMENT(numproc,taskid,table_index1,npar,dim,Rcurrent,Rfixed,L,meansquaredisp)
  USE MPI
implicit none
integer,intent(in):: numproc,taskid,table_index1(0:numproc-1,2)
integer,intent(in) :: npar,dim
real,intent(in):: Rcurrent(npar,dim),Rfixed(npar,dim)
real,intent(in):: L
real,intent(out):: meansquaredisp
real:: R(3)
real:: suma_local,suma_total
integer:: i
integer:: imin,imax
integer:: ierror

imin = table_index1(taskid,1)
imax = table_index1(taskid,2)

suma_local = 0.0
do i = imin,imax
	R = Rcurrent(i,:)-Rfixed(i,:)
	R = R - L*nint(R/L)
	suma_local = suma_local + (sum(R**2))
enddo


call MPI_REDUCE(suma_local,suma_total,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierror)

if (taskid==0) then
	meansquaredisp = suma_total/dble(npar)
else
	meansquaredisp = 0
endif
return
end subroutine MSDISPLACEMENT
!------------------------------------------------------------------------