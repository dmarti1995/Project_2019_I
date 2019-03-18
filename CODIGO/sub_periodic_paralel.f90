SUBROUTINE PBC(nat,L,X,numproc,taskid) !Periodic boundary conditions
implicit none
integer :: nat,ii,jj,numproc,taskid,partition
real    :: L, X(nat,3)
partition=nat/numproc
if (taskid.lt.(numproc-1)) then
   do ii=taskid*partition+1,partition*(taskid+1)
      do jj=1,3
         X(ii,jj)=X(ii,jj)-floor(X(ii,jj)/L)*L
      enddo
   enddo
else
   do ii=taskid*partition+1,nat
      do jj=1,3
         X(ii,jj)=X(ii,jj)-floor(X(ii,jj)/L)*L
      enddo
   enddo
endif
END SUBROUTINE
