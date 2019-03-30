SUBROUTINE PBC(nat,L,X,numproc,taskid,counts,displs,maxlength,&
           ierror) !Periodic boundary conditions
use mpi
implicit none
integer, parameter :: dp = 8
integer     :: nat,ii,jj,numproc,taskid,kk,counts(numproc)
integer     :: displs(numproc),maxlength,ierror,my_N_elem
real(dp)    :: L, X(nat,3),local_X(maxlength)
my_N_elem=counts(taskid+1)
kk=1+(displs(taskid+1)+1)/nat
jj=displs(taskid+1)+1-(kk-1)*nat
do ii=1,my_N_elem
      local_X(ii)=X(jj,kk)-floor(X(jj,kk)/L)*L
      if (jj==nat) then
         jj=0
         kk=kk+1
      endif
      jj=jj+1
enddo
call MPI_ALLGATHERV (local_X,my_N_elem , MPI_REAL8,  &
                     X, counts, displs, MPI_REAL8,   &
                     MPI_COMM_WORLD, ierror)
END SUBROUTINE
