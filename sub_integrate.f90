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
!	+ tau       :: tiempo caracteristico del termostato/barostato (tau = 20*dt)
!	+ Rcut      :: cut-off de la interaccion entre particulas
!	+ Text      :: temperatura del baño termico
!	+ Ppot      :: presion externa

subroutine leap_frog (r, v, vinf, vsup, t, F, Ppot, dt, tau, Rcut, L, Text, Pext)
implicit none
real, intent(inout) :: r(:,:), vinf(:,:), vsup(:,:), v(:,:), t, L
real, intent(in)    :: F(:,:), Ppot, dt, Rcut, tau, Text, Pext
real                :: Ekin, Pcalc, Tcalc, lambda, mu
integer             :: i, j, N

N = size(r,1)

! TERMOSTATO
Ekin   = 0.5 * sum(vinf**2)
Tcalc  = 2.0 * Ekin / (3.0*N)
lambda = sqrt(1.0 + dt/tau * (Text/Tcalc - 1.0))

! BAROSTATO
Pcalc = Ppot + 2.0*Ekin/(3.0 * L**3)
mu    = (1.0 + dt/tau*(Pcalc-Pext))**(1.0/3.0)

! Leap frog + reescalados
vinf = vsup
vsup = lambda * (vinf + F*dt)
r    = r + vsup*dt

! Reescalar r y L
L = L * mu
r = r * mu

! v(t)
v = 0.5 * (vinf + vsup)

! Condiciones periodicas
call PBC (N, L, r)

t = t + dt

end subroutine leap_frog
