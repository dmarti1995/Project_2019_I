
 !---------------make_celllist-----------------------------------
 !Subroutine that creates the cell lists
 subroutine make_celllist(npar,dim,ncells,soc,ipc,pos,celllist)
 implicit none
 integer :: npar,dim,ncells(3)
 real :: soc,pos(npar,dim)
 integer :: ipc(npar,dim+1)
 integer, intent(out) :: celllist(0:10,ncells(1),ncells(2),ncells(3))
 integer :: ipar, iinc, ix, iy, iz,ip
 celllist(:,:,:,:) = 0

 do ipar = 1, npar

    ix = ceiling(pos(ipar,1)/soc)    !calculation of the cell at which pertains each particle
    iy = ceiling(pos(ipar,2)/soc)
    iz = ceiling(pos(ipar,3)/soc)

    ip = celllist(0,ix,iy,iz) + 1

    celllist(0,ix,iy,iz) = ip               !counting the number of particles at cell ix,iy,iz
    celllist(ip, ix, iy,iz) = ipar          !at ip, we assign the index of the particle at cell ix,iy,iz

    ipc(ipar,1:4) = (/ix,iy,iz,ip/)        !we define ipc variable, which will allow us to actualise the cells
 enddo 

 end subroutine make_celllist

!subroutine update_cell_list: in case the cells configuration is modified, this subroutine
! allows to actualise the cell at wich partile ipar pertains
subroutine update_cell_list(dim,npar,ipc,ipar,ipcnew,celllist,ncells)
implicit none
integer :: npar,dim,ncells(3)
integer, intent(in) :: ipar, ipcnew(3)
integer, intent(inout) :: celllist(0:10,ncells(1),ncells(2),ncells(3))
integer :: jj, ipcold(4)
integer :: ipc(npar,dim+1)

jj=ipc(ipar,4)
ipcold(:)=ipc(ipar,1:4)    !old cell coordinates of particle ipar

celllist(jj,ipcold(1), ipcold(2),ipcold(3)) = &                                            
    celllist(celllist(0, ipcold(1), ipcold(2),ipcold(3)), ipcold(1), ipcold(2),ipcold(3) )
ipc((celllist(jj, ipcold(1), ipcold(2),ipcold(3))),4) = jj

celllist(0, ipcold(1), ipcold(2),ipcold(3)) = celllist(0, ipcold(1), ipcold(2),ipcold(3))-1  !we eliminate the particle that has moved from old cell
celllist(celllist(0, ipcold(1), ipcold(2),ipcold(3))+1,ipcold(1), ipcold(2),ipcold(3))=0

celllist(0, ipcnew(1), ipcnew(2),ipcnew(3)) = celllist(0, ipcnew(1), ipcnew(2),ipcnew(3))+1  !we actualise the particle that has moved to new cell
celllist(celllist(0, ipcnew(1), ipcnew(2),ipcnew(3)), ipcnew(1), ipcnew(2),ipcnew(3))=ipar

ipc(ipar,1:3)=ipcnew(:)
ipc(ipar,4)=celllist(0, ipcnew(1), ipcnew(2),ipcnew(3))

end subroutine update_cell_list


!subroutine to call in main program. Actualises all particles cells

subroutine actualise_cells(npar,dim,pos,celllist,ncells,soc,ipc)
integer :: npar,dim,ncells(3)
real :: soc,pos(npar,dim)
integer :: ipc(npar,dim+1),ipc_new(3),jj
integer, intent(out) :: celllist(0:10,ncells(1),ncells(2),ncells(3))

do jj=1,npar
    ipc_new(:) = ceiling(pos(jj,:)/soc)
    if (any(ipc_new .ne. ipc(jj,1:3))) then  ! update cell list
        call update_cell_list(dim,npar,ipc,jj,ipc_new,celllist,ncells)  !if a particle gets out of the cell we reintroduce it
    endif
enddo

return
end

