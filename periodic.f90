SUBROUTINE PBC(nat,L,X) !Periodic boundary conditions
implicit none
integer:: nat,ii,jj
real(8):: L, X(nat,3)
do ii=1,nat
   do jj=1,3
      X(ii,jj)=X(ii,jj)-nint(X(ii,jj)/L)*L
   enddo
enddo
END SUBROUTINE
