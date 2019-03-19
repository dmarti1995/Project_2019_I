SUBROUTINE PBC(nat,L,X,table,nproc,taskid) !Periodic boundary conditions
implicit none
integer :: nat,ii,jj,numproc,taskid,
integer :: table(0:nproc-1,2),imin,imax
real    :: L, X(nat,3)
imin=table(taskid,1)
imax=table(taskid,2)
do ii=imin,imax
   do jj=1,3
      X(ii,jj)=X(ii,jj)-floor(X(ii,jj)/L)*L
   enddo
enddo
END SUBROUTINE
