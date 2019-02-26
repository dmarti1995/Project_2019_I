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

subroutine leap_frog (npar,dim,r, v, vinf, vsup, t, F, Ppot, dt, taut,taup, Rcut, L, Text, Pext,pcalc,tcalc,ekin)
implicit none
real, intent(inout) :: r(npar,dim), vinf(npar,dim), vsup(npar,dim), v(npar,dim), t, L,ekin
real, intent(in)    :: F(npar,dim), Ppot, dt, Rcut, taut,taup, Text, Pext
real                :: Pcalc, Tcalc, lambda, mu
integer             :: i, j, npar,dim

! TERMOSTATO
Ekin   = 0.5 * sum(vinf**2)
Tcalc  = 2.0 * Ekin / (3.0*npar)
!print*, tcalc,ekin,vinf(3,1)
lambda = sqrt(1.0 + dt/taut * (Text/Tcalc - 1.0))
! BAROSTATO
Pcalc = (1d0/3d0)*(1/L**3)*(Ppot + 2d0*Ekin)
mu    = (1.0 + dt/taup*(Pcalc-Pext))**(1.0/3.0)
!print*,lambda,text,tcalc,mu,pcalc
! Leap frog + reescalados (velocidad)
!vinf = vsup
vsup = lambda * (vinf + F*dt)
r    = r + vsup*dt
! Reescalar r y L
L = L * mu
r = r * mu

! v(t)
v = 0.5 * (vinf + vsup)

vinf=vsup

! Condiciones periodicas
call PBC (npar, L, r)

t = t + dt

end subroutine leap_frog
