! INPUT y OUTPUT:
!	+ r(N,3)    :: posiciones de las particulas
!	+ v(N,3)    :: velocidades en el instante (t)
!	+ vinf(N,3) :: velocidades en el instante (t-dt/2)
!	+ vsup(N,3) :: velocidades en el instante (t+dt/2)
!	+ t         :: tiempo
!	+ L         :: longitud de la caja
! INPUT:
!	+ F(N,3)    :: fuerzas que actuan sobre cada particula en el instante (t)
!	+ Ppot      :: contribucion a la presion debida a la interaccion entre particulas
!	+ dt        :: paso de tiempo
!	+ taup/taut :: tiempo caracteristico del termostato/barostato
!	+ Rcut      :: cut-off de la interaccion entre particulas
!	+ Text      :: temperatura del baño termico
!	+ Pext      :: presion externa
! OUTPUT:
!	+ Pcalc     :: presion total del sistema
!	+ Tcalc     :: temperatura del sistema
!	+ Ekin      :: energia cinetica del sistema

subroutine leap_frog (N, r, v, vinf, vsup, t, F, Ppot, dt, taut, taup, Rcut, L, Text, Pext, Pcalc, Tcalc, Ekin)
implicit none
integer, intent(in) :: N
real, intent(inout) :: r(N,3), vinf(N,3), vsup(N,3), v(N,3), t, L,ekin
real, intent(in)    :: F(N,3), Ppot, dt, Rcut, taut, taup, Text, Pext
real                :: Pcalc, Tcalc, lambda, mu
integer             :: i, j


! TERMOSTATO
Ekin   = 0.5 * sum(vinf**2)
Tcalc  = 2.0 * Ekin / (3.0*N)
lambda = sqrt(1.0 + dt/taut * (Text/Tcalc - 1.0))

! BAROSTATO
Pcalc = (1.0/3.0)*(1.0/L**3)*(Ppot + 2.0*Ekin)
mu    = (1.0 + dt/taup*(Pcalc-Pext))**(1.0/3.0)

! Leap frog + reescalados (velocidad)
vsup = lambda * (vinf + F*dt)
r    = r + vsup*dt

! Reescalar r y L
L = L * mu
r = r * mu

! v(t)
v = 0.5 * (vinf + vsup)

vinf=vsup

t = t + dt

end subroutine leap_frog
