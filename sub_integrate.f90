subroutine v_verlet_step (N, r, v, t, dt, Rcut, L, F, press, Epot, Ekin)
implicit none
real, intent(inout) :: r(N,3), v(N,3), F(N,3), t
real, intent(in)    :: dt, Rcut, L
real, intent(out)   :: press, Epot, Ekin
integer             :: N

r = r + v*dt + 0.5*F*dt**2
v = v + 0.5*F*dt

call force (N, L, Rcut, r, v, F, press, Epot, Ekin)
v = v + 0.5*F*dt

t = t + dt

end subroutine v_verlet_step


subroutine thermostat (N, v, nu, temperature)
implicit none
real, intent(inout) :: v(N,3)
real, intent(in)    :: nu, temperature
integer, intent(in) :: N
real                :: u, z1, z2
integer             :: i, j

do i = 1, N
  call random_number(u)
  if (u < nu) then
    call normal_distr(0.0, sqrt(temperature), z1, z2)
    v(i,1) = z1
    v(i,2) = z2
    call normal_distr(0.0, sqrt(temperature), z1, z2)
    v(i,3) = z1
  endif
enddo

end subroutine thermostat


subroutine normal_distr(mu, sigma, z1, z2)
implicit none 
real, intent(in)  :: mu, sigma
real, intent(out) :: z1, z2
real, parameter   :: pi = acos(-1.0)
real              :: x, y

call random_number(x)
call random_number(y)
z1 = sqrt(-2.0 * log(x)) * cos(2.0 * pi * y)
z2 = sqrt(-2.0 * log(x)) * sin(2.0 * pi * y)
z1 = z1 * sigma + mu
z2 = z2 * sigma + mu

end subroutine normal_distr

