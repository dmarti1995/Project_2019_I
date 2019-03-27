subroutine ekinpress(numproc,taskid,table_index1,N,dim,L,x,v,F,reskin,respress)
use mpi
implicit none
integer,  intent(in) :: numproc,dim,N,taskid,table_index1(0:numproc-1,2)
integer :: iii,k,Nini,Nfin,ierror
real,intent(in) :: L
real,intent(in), dimension(N,dim) :: x,v,F
real,intent(out) :: reskin,respress
real :: press,ekin


Nini=table_index1(taskid,1)
Nfin=table_index1(taskid,2)
ekin=0.0
press=0.0
do iii=Nini,Nfin
  ! kinetic energy local
   ekin = ekin + 0.5*(sum(v(iii,:)**2))
  ! Potential pressure local
   do k=1,3
    press = press +x(iii,k)*F(iii,k)
   enddo
enddo

! total pressure local
press = press + 2.0*ekin
press = (1.0/3.0)*(1.0/L**3)*press

! kinetic energy global
call MPI_REDUCE(ekin, reskin, 1, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD,ierror)

! total pressure global
call MPI_REDUCE(press, respress, 1, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD,ierror)
end subroutine ekinpress



