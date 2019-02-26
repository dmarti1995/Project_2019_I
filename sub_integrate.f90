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
