subroutine ekinpress(numproc,taskid,table_index1,N,dim,L,x,v,F,reskin,respress)
use mpi
implicit none
integer, parameter                    :: dp = 8
integer,  intent(in)                  :: numproc,dim,N,taskid,table_index1(0:numproc-1,2)
integer                               :: iii,k,Nini,Nfin,ierror
real(dp),intent(in)                   :: L
real(dp),intent(in), dimension(N,dim) :: x,v,F
real(dp),intent(out)                  :: reskin,respress
real(dp)                              :: press,ekin


Nini=table_index1(taskid,1)
Nfin=table_index1(taskid,2)
ekin=0.0_dp
press=0.0_dp
do iii=Nini,Nfin
  ! kinetic energy local
   ekin = ekin + 0.5_dp*(sum(v(iii,:)**2.0_dp))
  ! Potential pressure local
   do k=1,3
    press = press +x(iii,k)*F(iii,k)
   enddo
enddo

! total pressure local
press = press + 2.0_dp*ekin
press = (1.0_dp/3.0_dp)*(1.0_dp/L**3.0_dp)*press

! kinetic energy global
call MPI_REDUCE(ekin, reskin, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD,ierror)

! total pressure global
call MPI_REDUCE(press, respress, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD,ierror)
end subroutine ekinpress



