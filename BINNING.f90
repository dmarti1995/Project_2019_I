module BINNING 
implicit none
contains
!SUBRUTINA QUE CALCULA EL VALOR ESPERAT I LA VARIANÇA D'UN VECTOR.
Subroutine Val_esp(vec,Val,Var)
real(8):: vec(:)
real(8):: Val,Var,suma
integer:: dim,i

dim = size(vec)


Val = sum(vec)/dble(dim)

suma = 0.0d0
do i = 1,dim
	suma = suma + (vec(i) - Val)**2

enddo

Var = suma/(dble(dim)*(dble(dim) - 1))

end Subroutine Val_esp

!SUBRUTINA QUE REALITZA EL BINNING DE L'ENERGIA I DE LA MAGNETITZACIÓ.
!AQUEST PROGRAMA REP numBINmi NOMBRE MÉS PETIT DE CAIXES A LES QUALS ES
!FARÀ EL BINNING, ARCHin NOM DE L'ARXIU DE DADES D'ON S'OBTÉ EL VALOR
!DE L'ENERGIA I DE LA MAGNETITZACIÓ PER CADA STEP DE MONTECARLO,ARCHout
!NOM DE L'ARXIU ON ES GUARDARÀ EL VALOR DEL VALOR ESPERAT DE E I M, I LES
!SEVES RESPECTIVES VARIANCES EN FUNCIÓ DE LA MIDA m DE LA CAIXA DEL BINNING.
!NMC NOMBRE DE PASSOS DE MONTECARLO QUE CONTÉ L'ARXIU (NOMBRE DE FILES)
!Nini NOMBE DE PASSOS DE MONTECARLO INICIALS QUE NO ES TENEN EN COMPTE PER
!FER EL BINNING
!L MIDA DEL SISTEMA

Subroutine BINNING_SUBR(numBINmin,ARCHin,ARCHout,NMC,Nini,L)
integer,intent(in):: numBINmin
character(len = *),intent(in):: ARCHin,ARCHout
integer:: NMC,Nini,L
real(8):: dum1,dum2,dum3,dum4,dum5,dum6,dum7
real(8),allocatable:: BINNEDP(:),AUXP(:)
real(8),allocatable:: BINNEDEp(:),AUXEp(:)
real(8),allocatable:: BINNEDEc(:),AUXEc(:)
real(8),allocatable:: BINNEDEm(:),AUXEm(:)
real(8):: sumaP,sumaEc,sumaEp,sumaEm
integer:: n,m,i,j,cont,numBIN,r
real(8):: ValP, VarP
real(8):: ValEc, VarEc
real(8):: ValEp, VarEp
real(8):: ValEm, VarEm


open(20,FILE=ARCHin)
open(432,FILE=ARCHout)
n = NMC - Nini
allocate(AUXP(n))
allocate(AUXEc(n))
allocate(AUXEp(n))
allocate(AUXEm(n))
cont = 0
do i =1,NMC
	read(20,*) dum1,dum2,dum3,dum4,dum5,dum6,dum7
	print*,dum7
	if (i > Nini) then
		cont = cont + 1
!ES GUARDA EN UN VECTOR 
		AUXEc(cont) = dum2
		AUXEp(cont) = dum3
		AUXEm(cont) = dum4
		AUXP(cont) = dum7
	endif
enddo

! CALCUL DEL VALOR ESPERAT I DE LA VARIANÇA DE E I M PER CAIXES DE MIDA m=1
m = 1
call Val_esp(AUXP,ValP,VarP)
call Val_esp(AUXEc,ValEc,VarEc)
call Val_esp(AUXEp,ValEp,VarEp)
call Val_esp(AUXEm,ValEm,VarEm)
write(432,*) m,ValP,sqrt(VarP),ValEc,sqrt(VarEc),ValEp,sqrt(VarEp),ValEm,sqrt(VarEm)

!FENT CAIXES CADA COP EL DOBLE DE GRANS ES FA EL PROCES DE BINNING

m = 2
numBIN = n/m
do while (numBIN > numBINmin)

! CONTROL DE SI EL NOMBRE DE DADES ES MULTIPLE DE LA MIDA DE LA CAIXA I, PER TANT,
! CONTROL DE SI HI S'HA D'AFEGIR UNA CAIXA EXTRA AMB LES DADES QUE SOBREN QUE TÉ
! UNA MIDA DIFERENT DE LA RESTA.
	r = size(AUXP) - numBIN*2

	if (r == 0) then
!SI TOTES LES DADES ENTREN EN CAIXES DE LA MATEIXA MIDA:
		allocate(BINNEDP(numBIN))
		allocate(BINNEDEc(numBIN))
		allocate(BINNEDEp(numBIN))
		allocate(BINNEDEm(numBIN))
	elseif (r/= 0) then
!SI SOBREN DADES QUE NO CABEN EN UNA CAIXA DE MIDA n/m ES GENERA UNA CAIXA NOVA
! DE MIDA VARIABLE SEGONS EL NOMBRE DE DADES QUE SOBREN.
		allocate(BINNEDP(numBIN+1))
		allocate(BINNEDEc(numBIN+1))
		allocate(BINNEDEp(numBIN+1))
		allocate(BINNEDEm(numBIN+1))
	endif


	cont = 0
	BINNEDP = 0.0d0
	BINNEDEc = 0.0d0
	BINNEDEp = 0.0d0
	BINNEDEm = 0.0d0
!S'OMPLEN LES CAIXES AMB PARELLES DE VALORS PROMITJATS DE CADA DUES CAIXES DEL 
! PAS ANTERIOR DEL BINNIG.
	do i = 1,numBIN
		sumaP = 0.0d0
		sumaEc = 0.0d0
		sumaEp = 0.0d0
		sumaEm = 0.0d0
		do j = 1,2
			cont = cont+1
			sumap = sumap + AUXp(cont)
			sumaEc = sumaEc + AUXEc(cont)
			sumaEp = sumaEp + AUXEp(cont)
			sumaEm = sumaEm + AUXEm(cont)
		enddo
		BINNEDP(i) = dble(sumaP)/2.0d0
		BINNEDEc(i) = dble(sumaEc)/2.0d0
		BINNEDEp(i) = dble(sumaEp)/2.0d0
		BINNEDEm(i) = dble(sumaEm)/2.0d0
	enddo
!PEL CAS EN QUE SOBRAVEN DADES ES CALCULA LA ÚLTIMA CAIXA:
	do i =1,r
		cont = cont+1
		BINNEDP(numBIN + 1) = BINNEDP(numBIN + 1) + AUXP(cont)/r
		BINNEDEc(numBIN + 1) = BINNEDEc(numBIN + 1) + AUXEc(cont)/r
		BINNEDEp(numBIN + 1) = BINNEDEp(numBIN + 1) + AUXEp(cont)/r
		BINNEDEm(numBIN + 1) = BINNEDEm(numBIN + 1) + AUXEm(cont)/r
	enddo
!ES CALCULA EL VALOR ESPERAT I LA VARIANÇA DE E I M AMB AQUEST NOU RESAMPLEIG
	call Val_esp(BINNEDP,ValP,VarP)
	call Val_esp(BINNEDEc,ValEc,VarEc)
	call Val_esp(BINNEDEp,ValEp,VarEp)
	call Val_esp(BINNEDEm,ValEm,VarEm)
	write(432,*) m,ValP,sqrt(VarP),ValEc,sqrt(VarEc),ValEp,sqrt(VarEp),ValEm,sqrt(VarEm)

	deallocate(AUXP)
	deallocate(AUXEc)
	deallocate(AUXEp)
	deallocate(AUXEm)
	allocate(AUXP(size(BINNEDP)))
	allocate(AUXEc(size(BINNEDEc)))
	allocate(AUXEp(size(BINNEDEp)))
	allocate(AUXEm(size(BINNEDEm)))
	AUXP = BINNEDP
	AUXEc = BINNEDEc
	AUXEp = BINNEDEp
	AUXEm = BINNEDEm
	deallocate(BINNEDP)
	deallocate(BINNEDEc)
	deallocate(BINNEDEp)
	deallocate(BINNEDEm)

	m = m*2
	numBIN = size(AUXP)/2
enddo

deallocate(AUXP)
deallocate(AUXEc)
deallocate(AUXEp)
deallocate(AUXEm)


close(432)

end Subroutine BINNING_SUBR
end module BINNING

!_________________________________________________________________________
!_________________________________________________________________________
!_________________________________________________________________________
!				PROGRAMA PRINCIPAL D'ANALISIS DE DADES
!_________________________________________________________________________
!_________________________________________________________________________
!_________________________________________________________________________
PROGRAM BINNING_MAIN_PROG
use BINNING
implicit none
!PARÂMETRES DEL PROGRAMA:
integer,parameter:: Npassos = 1000000 , Nini = 100 , Npart = 108

call BINNING_SUBR(2,"VELOCITY_VERLET_ANDERSEN_0.8.txt",&
					"BINNING_0.8.txt",Npassos,Nini,Npart)


END PROGRAM BINNING_MAIN_PROG