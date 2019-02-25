subroutine v_verlet_step (r, v, t, dt, Rcut, L, F)
implicit none
real, intent(inout) :: r(:,:), v(:,:), t
real, intent(in)    :: dt, Rcut, L
real, intent(out)   :: F(:,:)
integer             :: N

N = size(r,1)
call force (N, L, Rcut, r, F)
r = r + v*dt + 0.5*F*dt**2
v = v + 0.5*F*dt

call PBC (N, L, r)

call force (N, L, Rcut, r, F)
v = v + 0.5*F*dt

end subroutine v_verlet_step


subroutine thermostat (v, nu, temperature)
implicit none
real, intent(inout) :: v(:,:)
real, intent(in)    :: nu, temperature
real                :: u, z1, z2
integer             :: i, j

do i = 1, size(v,1)
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

