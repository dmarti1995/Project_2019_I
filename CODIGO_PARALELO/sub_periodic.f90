SUBROUTINE PBC(nat,L,X) !Periodic boundary conditions
implicit none
integer :: nat,ii,jj
real    :: L, X(nat,3)
do ii=1,nat
   do jj=1,3
      X(ii,jj)=X(ii,jj)-floor(X(ii,jj)/L)*L
   enddo
enddo
END SUBROUTINE
