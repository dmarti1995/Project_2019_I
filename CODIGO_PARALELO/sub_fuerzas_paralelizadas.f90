subroutine force(numproc,taskid,table_index2,N,Pair,dim,L,cutoff,x,v,Ftot,resepot)
use mpi
implicit none
integer,intent(in):: numproc,taskid,table_index2(0:numproc-1,2),Pair((N*(N-1))/2,2)
integer,intent(in):: N,dim
integer :: iii,i,j,k,Nini,Nfin, ierror
real,intent(in) :: L,cutoff
real :: dx,dy,dz,rij2,rij8,rij14,forcedist,potdist,rij12,rij6,epot
real,intent(in), dimension(N,dim) :: x,v
real, dimension(N,3) :: F
real,intent(out), dimension(N,dim) :: Ftot
real,intent(out) :: resepot

Nini=table_index2(taskid,1)
Nfin=table_index2(taskid,2)

F=0.0D0
epot=0.0D0
do iii=Nini,Nfin !Pair particles for this worker
   i=Pair(iii,1) !Particle i on the interaction
   j=Pair(iii,2) !Particle j on the interaction
   dx=x(i,1)-x(j,1)
   dx=dx-nint(dx/L)*L !PBC for distances
   dy=x(i,2)-x(j,2)
   dy=dy-nint(dy/L)*L
   dz=x(i,3)-x(j,3)
   dz=dz-nint(dz/L)*L
   rij2=dx**2.0D0+dy**2.0D0+dz**2.0D0
   if (rij2<cutoff**2.0) then !Computation only for shorter than cutoff 
     rij8=rij2**4.0
     rij14=rij2**7.0
     rij6=rij2**3.0
     rij12=rij2**6
     forcedist=(48.0D0/rij14 - 24.0D0/rij8)
     potdist=2d0*((1d0/rij12)-(1d0/rij6))
     F(i,1)=F(i,1)+forcedist*dx     !force matrix
     F(i,2)=F(i,2)+forcedist*dy
     F(i,3)=F(i,3)+forcedist*dz
     epot=epot+2*potdist              !potential energy calculation
     F(j,1)=F(j,1)-forcedist*dx     !force matrix
     F(j,2)=F(j,2)-forcedist*dy
     F(j,3)=F(j,3)-forcedist*dz
   endif
enddo

! Suma de todas las fuerzas
call MPI_ALLREDUCE(F,Ftot,size(Ftot),MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierror)

! potential energy global
call MPI_REDUCE(epot, resepot, 1, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD,ierror)

!if (taskid==0) write(6,*) "Resepot", resepot
end subroutine force


