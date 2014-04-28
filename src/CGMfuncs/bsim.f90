!>  @file  bsim.f90
!!  @brief Basic subprograms for a sharp interface method (BSIM) for Cartesian grid data structure
!< 

function sim_geth(phi)
  implicit none
  real          :: sim_geth
  real          :: phi
  sim_geth = 1.0
  if( phi < 0 ) then
    sim_geth = 0.0
  endif
  return
end function sim_geth

function sim_geth_s(phi, d)
  implicit none
  real          :: sim_geth_s
  real          :: phi
	real					:: d
	sim_geth_s = 1.0
	if( phi < -d ) then
		sim_geth_s = 0.0
	else if( phi > d ) then
		sim_geth_s = 1.0
	else
		sim_geth_s = 0.5*(phi + d)/d
	end if
  return
end function sim_geth_s

function sim_geth_2(phi, d)
  implicit none
  real          :: sim_geth_2
  real          :: phi
	real					:: d
	real					:: pi
	pi = atan(1.0)*4.0
	sim_geth_2 = 0.5 + 0.5*max(-1.0, min(1.0, phi/d + 1.0/pi*sin(pi*phi/d)))
  return
end function sim_geth_2

function sim_getd(phi0, phi1)
  implicit none
  real          :: sim_getd
  real          :: phi0, phi1
  sim_getd = abs(phi0)/(abs(phi0) + abs(phi1))
  if( phi0 >= 0 .and. phi1 >= 0 ) then
    sim_getd = 1.0
  endif
  if( phi0 < 0 .and. phi1 < 0 ) then
    sim_getd = 1.0
  endif
  return
end function sim_getd

function sim_getminmod(a, b)
  implicit none
  real          :: sim_getminmod
  real          :: a, b
!  sim_getminmod = 0.5*(sign(a, a) + sign(b, b))*min(abs(a), abs(b))
  sim_getminmod = 0.0
  if( abs(a) .lt. abs(b) ) then
    sim_getminmod = a
  else if( abs(a) .gt. abs(b) ) then
    sim_getminmod = b
  endif
  return
end function sim_getminmod

function sim_getupwind(u, n, p)
  implicit none
  real          :: sim_getupwind
  real          :: u
  real          :: n, p
  sim_getupwind = 0.0
  if( u .gt. 0 ) then
    sim_getupwind = n
  else if( u .lt. 0 ) then
    sim_getupwind = p
  endif
  return
end function sim_getupwind

function sim_getweno3(v1, v2, v3)
  implicit none
  real          :: sim_getweno3
  real          :: v1, v2, v3
  real          :: r
  real          :: w
  real          :: eps
  eps = 0.001
  r = (eps + (v2-v1)**2)/(eps + (v3-v2)**2)
  w = 1.0/(1.0 + 2*r*r)
  sim_getweno3 = 0.5*(v2 + v3) - 0.5*w*(v1 - 2.0*v2 + v3)
  return
end function sim_getweno3

function sim_getux(x, y, z, Omega, zmin, rmin)
	implicit none
	real					:: sim_getux
	real					:: x, y, z, Omega, zmin, rmin
	real					:: r2, rl, theta
	real					:: nx, ny, nz
	real					:: Ux, Uy, Uz
	r2 = x*x + y*y
	rl = sqrt(r2)
	theta = atan(x)
	nx = x/rl
	ny = y/rl
	nz = z/rl
	Ux =-rl*Omega*ny
	Uy =+rl*Omega*nx
	Uz = 0.0
	Ux =-Omega*y
	Uy =+Omega*x
	Uz = 0.0
	if( rl < rmin .and. z > zmin ) then
		Ux = 0.0
		Uy = 0.0
		Uz = 0.0
	end if
	sim_getux = Ux
	return
end function sim_getux

function sim_getuy(x, y, z, Omega, zmin, rmin)
	implicit none
	real					:: sim_getuy
	real					:: x, y, z, Omega, zmin, rmin
	real					:: r2, rl, theta
	real					:: nx, ny, nz
	real					:: Ux, Uy, Uz
	r2 = x*x + y*y
	rl = sqrt(r2)
	theta = atan(x)
	nx = x/rl
	ny = y/rl
	nz = z/rl
	Ux =-rl*Omega*ny
	Uy =+rl*Omega*nx
	Uz = 0.0
	Ux =-Omega*y
	Uy =+Omega*x
	Uz = 0.0
	if( rl < rmin .and. z > zmin ) then
		Ux = 0.0
		Uy = 0.0
		Uz = 0.0
	end if
	sim_getuy = Uy
	return
end function sim_getuy

function sim_getuz(x, y, z, Omega, zmin, rmin)
	implicit none
	real					:: sim_getuz
	real					:: x, y, z, Omega, zmin, rmin
	real					:: r2, rl, theta
	real					:: nx, ny, nz
	real					:: Ux, Uy, Uz
	r2 = x*x + y*y
	rl = sqrt(r2)
	theta = atan(x)
	nx = x/rl
	ny = y/rl
	nz = z/rl
	Ux =-rl*Omega*ny
	Uy =+rl*Omega*nx
	Uz = 0.0
	Ux =-Omega*y
	Uy =+Omega*x
	Uz = 0.0
	if( rl < rmin .and. z > zmin ) then
		Ux = 0.0
		Uy = 0.0
		Uz = 0.0
	end if
	sim_getuz = Uz
	return
end function sim_getuz

subroutine sim_calc_c( &
                fc &
              , f &
              , v &
              , phi &
              , psi &
							, org &
              , sz, g)
  implicit none
  integer                  :: i, j, k
  integer                  :: ix, jx, kx
  integer                  :: g
  integer, dimension(3)    :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: fc
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: f
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:6)  :: v
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: phi, psi
	real, dimension(3)			:: org
  real                    :: f00
  real                    :: f11, f13, f12, f14, f15, f16
  real                    :: f21, f23, f22, f24, f25, f26
  real                    :: f1, f3, f2, f4, f5, f6
  real                    :: vx1, vx3, vy2, vy4, vz5, vz6
  real                    :: q1, q3, q2, q4, q5, q6
  real                    :: psi0, psi1, psi2, psi3, psi4, psi5, psi6
  real                    :: Hpsi0, Hpsi1, Hpsi2, Hpsi3, Hpsi4, Hpsi5, Hpsi6
  real                    :: dpsi1, dpsi2, dpsi3, dpsi4, dpsi5, dpsi6
	real										:: x, y, z, xi, yi, zi
	real										:: Uw
  real                    :: sim_getweno3
  real                    :: sim_geth
  real                    :: sim_getd
	real										:: sim_getux
	real										:: sim_getuy
	real										:: sim_getuz
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
!$omp parallel private(i, j, k) &
!$omp          ,private(f00) &
!$omp          ,private(f11, f13, f12, f14, f15, f16) &
!$omp          ,private(f21, f23, f22, f24, f25, f26) &
!$omp          ,private(f1, f3, f2, f4, f5, f6) &
!$omp          ,private(vx1, vx3, vy2, vy4, vz5, vz6) &
!$omp          ,private(q1, q3, q2, q4, q5, q6) &
!$omp          ,private(psi0, psi1, psi2, psi3, psi4, psi5, psi6) &
!$omp          ,private(Hpsi0, Hpsi1, Hpsi2, Hpsi3, Hpsi4, Hpsi5, Hpsi6) &
!$omp          ,private(dpsi1, dpsi2, dpsi3, dpsi4, dpsi5, dpsi6) &
!$omp					 ,private(x, y, z, xi, yi, zi) &
!$omp          ,private(Uw) 
!$omp do schedule(dynamic, 1)
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
		x = org(1) + float(i) - 0.5
		y = org(2) + float(j) - 0.5
		z = org(3) + float(k) - 0.5

    vx1 = v(i, j, k, 1)
    vx3 = v(i, j, k, 3)
    vy2 = v(i, j, k, 2)
    vy4 = v(i, j, k, 4)
    vz5 = v(i, j, k, 5)
    vz6 = v(i, j, k, 6)

    psi0 = psi(i, j, k)
    psi1 = psi(i+1, j, k)
    psi3 = psi(i-1, j, k)
    psi2 = psi(i, j+1, k)
    psi4 = psi(i, j-1, k)
    psi5 = psi(i, j, k+1)
    psi6 = psi(i, j, k-1)

    Hpsi0 = sim_geth(psi0)
    Hpsi1 = sim_geth(psi1)
    Hpsi3 = sim_geth(psi3)
    Hpsi2 = sim_geth(psi2)
    Hpsi4 = sim_geth(psi4)
    Hpsi5 = sim_geth(psi5)
    Hpsi6 = sim_geth(psi6)

    dpsi1 = sim_getd(psi0, psi1)
    dpsi3 = sim_getd(psi0, psi3)
    dpsi2 = sim_getd(psi0, psi2)
    dpsi4 = sim_getd(psi0, psi4)
    dpsi5 = sim_getd(psi0, psi5)
    dpsi6 = sim_getd(psi0, psi6)

    f00 = f(i, j, k)
    f11 = f(i+1, j, k)
    f13 = f(i-1, j, k)
    f12 = f(i, j+1, k)
    f14 = f(i, j-1, k)
    f15 = f(i, j, k+1)
    f16 = f(i, j, k-1)
		if( psi(i, j, k) >= 0 ) then
		else
			if( psi(i+1, j, k) >= 0 ) then
				xi = x + dpsi1
				yi = y
				zi = z
				f11 = 0.0
			endif
			if( psi(i-1, j, k) >= 0 ) then
				xi = x - dpsi3
				yi = y
				zi = z
				f13 = 0.0
			endif
			if( psi(i, j+1, k) >= 0 ) then
				xi = x
				yi = y + dpsi2
				zi = z
				f12 = 0.0
			endif
			if( psi(i, j-1, k) >= 0 ) then
				xi = x
				yi = y - dpsi4
				zi = z
				f14 = 0.0
			endif
			if( psi(i, j, k+1) >= 0 ) then
				xi = x
				yi = y
				zi = z + dpsi5
				f15 = 0.0
			endif
			if( psi(i, j, k-1) >= 0 ) then
				xi = x
				yi = y
				zi = z - dpsi6
				f16 = 0.0
			endif
		endif

		f1 = 0.5*(f00 + f11)
		f3 = 0.5*(f00 + f13)
		f2 = 0.5*(f00 + f12)
		f4 = 0.5*(f00 + f14)
		f5 = 0.5*(f00 + f15)
		f6 = 0.5*(f00 + f16)

    q1 = f1*vx1
    q3 = f3*vx3
    q2 = f2*vy2
    q4 = f4*vy4
    q5 = f5*vz5
    q6 = f6*vz6

    if( psi1 >= 0 .and. psi3 >= 0 ) then
      q1 = 0.0
      q3 = 0.0
    endif
    if( psi2 >= 0 .and. psi4 >= 0 ) then
      q2 = 0.0
      q4 = 0.0
    endif
    if( psi5 >= 0 .and. psi6 >= 0 ) then
      q5 = 0.0
      q6 = 0.0
    endif

    fc(i, j, k) = q1 - q3 + q2 - q4 + q5 - q6

  end do
  end do
  end do
!$omp end do
!$omp end parallel
end subroutine sim_calc_c

subroutine sim_calc_c_u( &
                fc &
              , f &
              , v &
              , phi &
              , psi &
							, dir &
							, Omega &
							, zmin &
							, zmax &
							, rmin &
							, rmax &
							, org &
              , sz, g)
  implicit none
  integer                  :: i, j, k
  integer                  :: ix, jx, kx
  integer                  :: g
  integer, dimension(3)    :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: fc
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: f
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:6)  :: v
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: phi, psi
	integer									:: dir
	real										:: Omega, zmin, zmax, rmin, rmax
	real, dimension(3)			:: org
  real                    :: f00
  real                    :: f11, f13, f12, f14, f15, f16
  real                    :: f21, f23, f22, f24, f25, f26
  real                    :: f1, f3, f2, f4, f5, f6
  real                    :: vx1, vx3, vy2, vy4, vz5, vz6
  real                    :: q1, q3, q2, q4, q5, q6
  real                    :: psi0, psi1, psi2, psi3, psi4, psi5, psi6
  real                    :: Hpsi0, Hpsi1, Hpsi2, Hpsi3, Hpsi4, Hpsi5, Hpsi6
  real                    :: dpsi1, dpsi2, dpsi3, dpsi4, dpsi5, dpsi6
	real										:: x, y, z, xi, yi, zi
	real										:: Uw
  real                    :: sim_getweno3
  real                    :: sim_geth
  real                    :: sim_getd
	real										:: sim_getux
	real										:: sim_getuy
	real										:: sim_getuz
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
!$omp parallel private(i, j, k) &
!$omp          ,private(f00) &
!$omp          ,private(f11, f13, f12, f14, f15, f16) &
!$omp          ,private(f21, f23, f22, f24, f25, f26) &
!$omp          ,private(f1, f3, f2, f4, f5, f6) &
!$omp          ,private(vx1, vx3, vy2, vy4, vz5, vz6) &
!$omp          ,private(q1, q3, q2, q4, q5, q6) &
!$omp          ,private(psi0, psi1, psi2, psi3, psi4, psi5, psi6) &
!$omp          ,private(Hpsi0, Hpsi1, Hpsi2, Hpsi3, Hpsi4, Hpsi5, Hpsi6) &
!$omp          ,private(dpsi1, dpsi2, dpsi3, dpsi4, dpsi5, dpsi6) &
!$omp					 ,private(x, y, z, xi, yi, zi) &
!$omp          ,private(Uw) 
!$omp do schedule(dynamic, 1)
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
		x = org(1) + float(i) - 0.5
		y = org(2) + float(j) - 0.5
		z = org(3) + float(k) - 0.5

    vx1 = v(i, j, k, 1)
    vx3 = v(i, j, k, 3)
    vy2 = v(i, j, k, 2)
    vy4 = v(i, j, k, 4)
    vz5 = v(i, j, k, 5)
    vz6 = v(i, j, k, 6)

    psi0 = psi(i, j, k)
    psi1 = psi(i+1, j, k)
    psi3 = psi(i-1, j, k)
    psi2 = psi(i, j+1, k)
    psi4 = psi(i, j-1, k)
    psi5 = psi(i, j, k+1)
    psi6 = psi(i, j, k-1)

    Hpsi0 = sim_geth(psi0)
    Hpsi1 = sim_geth(psi1)
    Hpsi3 = sim_geth(psi3)
    Hpsi2 = sim_geth(psi2)
    Hpsi4 = sim_geth(psi4)
    Hpsi5 = sim_geth(psi5)
    Hpsi6 = sim_geth(psi6)

    dpsi1 = sim_getd(psi0, psi1)
    dpsi3 = sim_getd(psi0, psi3)
    dpsi2 = sim_getd(psi0, psi2)
    dpsi4 = sim_getd(psi0, psi4)
    dpsi5 = sim_getd(psi0, psi5)
    dpsi6 = sim_getd(psi0, psi6)

    f00 = f(i, j, k)
    f11 = f(i+1, j, k)
    f13 = f(i-1, j, k)
    f12 = f(i, j+1, k)
    f14 = f(i, j-1, k)
    f15 = f(i, j, k+1)
    f16 = f(i, j, k-1)
		if( psi(i, j, k) >= 0 ) then
		else
			if( psi(i+1, j, k) >= 0 ) then
				xi = x + dpsi1
				yi = y
				zi = z
				if( dir == 0 ) then
					Uw = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 1 ) then
					Uw = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 2 ) then
					Uw = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
				endif
				f11 = f00 + (Uw - f00)/dpsi1
			endif
			if( psi(i-1, j, k) >= 0 ) then
				xi = x - dpsi3
				yi = y
				zi = z
				if( dir == 0 ) then
					Uw = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 1 ) then
					Uw = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 2 ) then
					Uw = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
				endif
				f13 = f00 + (Uw - f00)/dpsi3
			endif
			if( psi(i, j+1, k) >= 0 ) then
				xi = x
				yi = y + dpsi2
				zi = z
				if( dir == 0 ) then
					Uw = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 1 ) then
					Uw = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 2 ) then
					Uw = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
				endif
				f12 = f00 + (Uw - f00)/dpsi2
			endif
			if( psi(i, j-1, k) >= 0 ) then
				xi = x
				yi = y - dpsi4
				zi = z
				if( dir == 0 ) then
					Uw = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 1 ) then
					Uw = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 2 ) then
					Uw = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
				endif
				f14 = f00 + (Uw - f00)/dpsi4
			endif
			if( psi(i, j, k+1) >= 0 ) then
				xi = x
				yi = y
				zi = z + dpsi5
				if( dir == 0 ) then
					Uw = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 1 ) then
					Uw = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 2 ) then
					Uw = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
				endif
				f15 = f00 + (Uw - f00)/dpsi5
			endif
			if( psi(i, j, k-1) >= 0 ) then
				xi = x
				yi = y
				zi = z - dpsi6
				if( dir == 0 ) then
					Uw = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 1 ) then
					Uw = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 2 ) then
					Uw = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
				endif
				f16 = f00 + (Uw - f00)/dpsi6
			endif
		endif

		f1 = 0.5*(f00 + f11)
		f3 = 0.5*(f00 + f13)
		f2 = 0.5*(f00 + f12)
		f4 = 0.5*(f00 + f14)
		f5 = 0.5*(f00 + f15)
		f6 = 0.5*(f00 + f16)

    q1 = f1*vx1
    q3 = f3*vx3
    q2 = f2*vy2
    q4 = f4*vy4
    q5 = f5*vz5
    q6 = f6*vz6

    if( psi1 >= 0 .and. psi3 >= 0 ) then
      q1 = 0.0
      q3 = 0.0
    endif
    if( psi2 >= 0 .and. psi4 >= 0 ) then
      q2 = 0.0
      q4 = 0.0
    endif
    if( psi5 >= 0 .and. psi6 >= 0 ) then
      q5 = 0.0
      q6 = 0.0
    endif

    fc(i, j, k) = q1 - q3 + q2 - q4 + q5 - q6

  end do
  end do
  end do
!$omp end do
!$omp end parallel
end subroutine sim_calc_c_u

subroutine sim_calc_c_u_u1( &
                fc &
              , f &
              , v &
              , phi &
              , psi &
							, dir &
							, Omega &
							, zmin &
							, zmax &
							, rmin &
							, rmax &
							, org &
              , sz, g)
  implicit none
  integer                  :: i, j, k
  integer                  :: ix, jx, kx
  integer                  :: g
  integer, dimension(3)    :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: fc
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: f
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:6)  :: v
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: phi, psi
	integer									:: dir
	real										:: Omega, zmin, zmax, rmin, rmax
	real, dimension(3)			:: org
  real                    :: f00
  real                    :: f11, f13, f12, f14, f15, f16
  real                    :: f21, f23, f22, f24, f25, f26
  real                    :: f1, f3, f2, f4, f5, f6
  real                    :: vx1, vx3, vy2, vy4, vz5, vz6
  real                    :: q1, q3, q2, q4, q5, q6
  real                    :: psi0, psi1, psi2, psi3, psi4, psi5, psi6
  real                    :: Hpsi0, Hpsi1, Hpsi2, Hpsi3, Hpsi4, Hpsi5, Hpsi6
  real                    :: dpsi1, dpsi2, dpsi3, dpsi4, dpsi5, dpsi6
	real										:: x, y, z, xi, yi, zi
	real										:: Uw
  real                    :: sim_getweno3
  real                    :: sim_geth
  real                    :: sim_getd
	real										:: sim_getux
	real										:: sim_getuy
	real										:: sim_getuz
	real										:: sim_getupwind
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
!$omp parallel private(i, j, k) &
!$omp          ,private(f00) &
!$omp          ,private(f11, f13, f12, f14, f15, f16) &
!$omp          ,private(f21, f23, f22, f24, f25, f26) &
!$omp          ,private(f1, f3, f2, f4, f5, f6) &
!$omp          ,private(vx1, vx3, vy2, vy4, vz5, vz6) &
!$omp          ,private(q1, q3, q2, q4, q5, q6) &
!$omp          ,private(psi0, psi1, psi2, psi3, psi4, psi5, psi6) &
!$omp          ,private(Hpsi0, Hpsi1, Hpsi2, Hpsi3, Hpsi4, Hpsi5, Hpsi6) &
!$omp          ,private(dpsi1, dpsi2, dpsi3, dpsi4, dpsi5, dpsi6) &
!$omp					 ,private(x, y, z, xi, yi, zi) &
!$omp          ,private(Uw) 
!$omp do schedule(dynamic, 1)
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
		x = org(1) + float(i) - 0.5
		y = org(2) + float(j) - 0.5
		z = org(3) + float(k) - 0.5

    vx1 = v(i, j, k, 1)
    vx3 = v(i, j, k, 3)
    vy2 = v(i, j, k, 2)
    vy4 = v(i, j, k, 4)
    vz5 = v(i, j, k, 5)
    vz6 = v(i, j, k, 6)

    psi0 = psi(i, j, k)
    psi1 = psi(i+1, j, k)
    psi3 = psi(i-1, j, k)
    psi2 = psi(i, j+1, k)
    psi4 = psi(i, j-1, k)
    psi5 = psi(i, j, k+1)
    psi6 = psi(i, j, k-1)

    Hpsi0 = sim_geth(psi0)
    Hpsi1 = sim_geth(psi1)
    Hpsi3 = sim_geth(psi3)
    Hpsi2 = sim_geth(psi2)
    Hpsi4 = sim_geth(psi4)
    Hpsi5 = sim_geth(psi5)
    Hpsi6 = sim_geth(psi6)

    dpsi1 = sim_getd(psi0, psi1)
    dpsi3 = sim_getd(psi0, psi3)
    dpsi2 = sim_getd(psi0, psi2)
    dpsi4 = sim_getd(psi0, psi4)
    dpsi5 = sim_getd(psi0, psi5)
    dpsi6 = sim_getd(psi0, psi6)

    f00 = f(i, j, k)
    f11 = f(i+1, j, k)
    f13 = f(i-1, j, k)
    f12 = f(i, j+1, k)
    f14 = f(i, j-1, k)
    f15 = f(i, j, k+1)
    f16 = f(i, j, k-1)
		if( psi(i, j, k) >= 0 ) then
		else
			if( psi(i+1, j, k) >= 0 ) then
				xi = x + dpsi1
				yi = y
				zi = z
				if( dir == 0 ) then
					Uw = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 1 ) then
					Uw = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 2 ) then
					Uw = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
				endif
				f11 = f00 + (Uw - f00)/dpsi1
				f11 = f00
				f11 = Uw
			endif
			if( psi(i-1, j, k) >= 0 ) then
				xi = x - dpsi3
				yi = y
				zi = z
				if( dir == 0 ) then
					Uw = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 1 ) then
					Uw = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 2 ) then
					Uw = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
				endif
				f13 = f00 + (Uw - f00)/dpsi3
				f13 = f00
				f13 = Uw
			endif
			if( psi(i, j+1, k) >= 0 ) then
				xi = x
				yi = y + dpsi2
				zi = z
				if( dir == 0 ) then
					Uw = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 1 ) then
					Uw = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 2 ) then
					Uw = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
				endif
				f12 = f00 + (Uw - f00)/dpsi2
				f12 = f00
				f12 = Uw
			endif
			if( psi(i, j-1, k) >= 0 ) then
				xi = x
				yi = y - dpsi4
				zi = z
				if( dir == 0 ) then
					Uw = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 1 ) then
					Uw = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 2 ) then
					Uw = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
				endif
				f14 = f00 + (Uw - f00)/dpsi4
				f14 = f00
				f14 = Uw
			endif
			if( psi(i, j, k+1) >= 0 ) then
				xi = x
				yi = y
				zi = z + dpsi5
				if( dir == 0 ) then
					Uw = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 1 ) then
					Uw = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 2 ) then
					Uw = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
				endif
				f15 = f00 + (Uw - f00)/dpsi5
				f15 = f00
				f15 = Uw
			endif
			if( psi(i, j, k-1) >= 0 ) then
				xi = x
				yi = y
				zi = z - dpsi6
				if( dir == 0 ) then
					Uw = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 1 ) then
					Uw = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 2 ) then
					Uw = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
				endif
				f16 = f00 + (Uw - f00)/dpsi6
				f16 = f00
				f16 = Uw
			endif
		endif

		f1 = sim_getupwind(vx1, f00, f11)
		f3 = sim_getupwind(vx3, f13, f00)
		f2 = sim_getupwind(vy2, f00, f12)
		f4 = sim_getupwind(vy4, f14, f00)
		f5 = sim_getupwind(vz5, f00, f15)
		f6 = sim_getupwind(vz6, f16, f00)

    q1 = f1*vx1
    q3 = f3*vx3
    q2 = f2*vy2
    q4 = f4*vy4
    q5 = f5*vz5
    q6 = f6*vz6

    if( psi1 >= 0 .and. psi3 >= 0 ) then
      q1 = 0.0
      q3 = 0.0
    endif
    if( psi2 >= 0 .and. psi4 >= 0 ) then
      q2 = 0.0
      q4 = 0.0
    endif
    if( psi5 >= 0 .and. psi6 >= 0 ) then
      q5 = 0.0
      q6 = 0.0
    endif

    fc(i, j, k) = q1 - q3 + q2 - q4 + q5 - q6

  end do
  end do
  end do
!$omp end do
!$omp end parallel
end subroutine sim_calc_c_u_u1

subroutine sim_calc_c_u_quick( &
                fc &
              , f &
              , v &
              , phi &
              , psi &
							, dir &
							, Omega &
							, zmin &
							, zmax &
							, rmin &
							, rmax &
							, org &
              , sz, g)
  implicit none
  integer                  :: i, j, k
  integer                  :: ix, jx, kx
  integer                  :: g
  integer, dimension(3)    :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: fc
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: f
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:6)  :: v
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: phi, psi
	integer									:: dir
	real										:: Omega, zmin, zmax, rmin, rmax
	real, dimension(3)			:: org
  real                    :: f00
  real                    :: f11, f13, f12, f14, f15, f16
  real                    :: f21, f23, f22, f24, f25, f26
  real                    :: f1, f3, f2, f4, f5, f6
  real                    :: vx1, vx3, vy2, vy4, vz5, vz6
  real                    :: q1, q3, q2, q4, q5, q6
  real                    :: psi0, psi1, psi2, psi3, psi4, psi5, psi6
  real                    :: Hpsi0, Hpsi1, Hpsi2, Hpsi3, Hpsi4, Hpsi5, Hpsi6
  real                    :: dpsi1, dpsi2, dpsi3, dpsi4, dpsi5, dpsi6
	real										:: x, y, z, xi, yi, zi
	real										:: Uw
  real                    :: sim_getweno3
  real                    :: sim_geth
  real                    :: sim_getd
	real										:: sim_getux
	real										:: sim_getuy
	real										:: sim_getuz
  real                    :: sim_getupwind
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
!$omp parallel private(i, j, k) &
!$omp          ,private(f00) &
!$omp          ,private(f11, f13, f12, f14, f15, f16) &
!$omp          ,private(f21, f23, f22, f24, f25, f26) &
!$omp          ,private(f1, f3, f2, f4, f5, f6) &
!$omp          ,private(vx1, vx3, vy2, vy4, vz5, vz6) &
!$omp          ,private(q1, q3, q2, q4, q5, q6) &
!$omp          ,private(psi0, psi1, psi2, psi3, psi4, psi5, psi6) &
!$omp          ,private(Hpsi0, Hpsi1, Hpsi2, Hpsi3, Hpsi4, Hpsi5, Hpsi6) &
!$omp          ,private(dpsi1, dpsi2, dpsi3, dpsi4, dpsi5, dpsi6) &
!$omp					 ,private(x, y, z, xi, yi, zi) &
!$omp          ,private(Uw) 
!$omp do schedule(dynamic, 1)
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
		x = org(1) + float(i) - 0.5
		y = org(2) + float(j) - 0.5
		z = org(3) + float(k) - 0.5

    vx1 = v(i, j, k, 1)
    vx3 = v(i, j, k, 3)
    vy2 = v(i, j, k, 2)
    vy4 = v(i, j, k, 4)
    vz5 = v(i, j, k, 5)
    vz6 = v(i, j, k, 6)

    psi0 = psi(i, j, k)
    psi1 = psi(i+1, j, k)
    psi3 = psi(i-1, j, k)
    psi2 = psi(i, j+1, k)
    psi4 = psi(i, j-1, k)
    psi5 = psi(i, j, k+1)
    psi6 = psi(i, j, k-1)

    Hpsi0 = sim_geth(psi0)
    Hpsi1 = sim_geth(psi1)
    Hpsi3 = sim_geth(psi3)
    Hpsi2 = sim_geth(psi2)
    Hpsi4 = sim_geth(psi4)
    Hpsi5 = sim_geth(psi5)
    Hpsi6 = sim_geth(psi6)

    dpsi1 = sim_getd(psi0, psi1)
    dpsi3 = sim_getd(psi0, psi3)
    dpsi2 = sim_getd(psi0, psi2)
    dpsi4 = sim_getd(psi0, psi4)
    dpsi5 = sim_getd(psi0, psi5)
    dpsi6 = sim_getd(psi0, psi6)

    f00 = f(i, j, k)
    f11 = f(i+1, j, k)
    f13 = f(i-1, j, k)
    f12 = f(i, j+1, k)
    f14 = f(i, j-1, k)
    f15 = f(i, j, k+1)
    f16 = f(i, j, k-1)
    f21 = f(i+2, j, k)
    f23 = f(i-2, j, k)
    f22 = f(i, j+2, k)
    f24 = f(i, j-2, k)
    f25 = f(i, j, k+2)
    f26 = f(i, j, k-2)
		if( psi(i+1, j, k) >= 0 ) then
			f11 = f00
			f21 = f00
		endif
		if( psi(i-1, j, k) >= 0 ) then
			f13 = f00
			f23 = f00
		endif
		if( psi(i, j+1, k) >= 0 ) then
			f12 = f00
			f22 = f00
		endif
		if( psi(i, j-1, k) >= 0 ) then
			f14 = f00
			f24 = f00
		endif
		if( psi(i, j, k+1) >= 0 ) then
			f15 = f00
			f25 = f00
		endif
		if( psi(i, j, k-1) >= 0 ) then
			f16 = f00
			f26 = f00
		endif

		if( psi(i+1, j, k) >= 0 ) then
			xi = x + dpsi1
			yi = y
			zi = z
			if( dir == 0 ) then
				Uw = sim_getux(xi, yi, zi, Omega, zmin, rmin)
			else if ( dir == 1 ) then
				Uw = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
			else if ( dir == 2 ) then
				Uw = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
			endif
			f11 = Uw
			f21 = Uw
		endif
		if( psi(i-1, j, k) >= 0 ) then
			xi = x - dpsi3
			yi = y
			zi = z
			if( dir == 0 ) then
				Uw = sim_getux(xi, yi, zi, Omega, zmin, rmin)
			else if ( dir == 1 ) then
				Uw = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
			else if ( dir == 2 ) then
				Uw = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
			endif
			f13 = Uw
			f23 = Uw
		endif
		if( psi(i, j+1, k) >= 0 ) then
			xi = x
			yi = y + dpsi2
			zi = z
			if( dir == 0 ) then
				Uw = sim_getux(xi, yi, zi, Omega, zmin, rmin)
			else if ( dir == 1 ) then
				Uw = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
			else if ( dir == 2 ) then
				Uw = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
			endif
			f12 = Uw
			f22 = Uw
		endif
		if( psi(i, j-1, k) >= 0 ) then
			xi = x
			yi = y - dpsi4
			zi = z
			if( dir == 0 ) then
				Uw = sim_getux(xi, yi, zi, Omega, zmin, rmin)
			else if ( dir == 1 ) then
				Uw = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
			else if ( dir == 2 ) then
				Uw = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
			endif
			f14 = Uw
			f24 = Uw
		endif
		if( psi(i, j, k+1) >= 0 ) then
			xi = x
			yi = y
			zi = z + dpsi5
			if( dir == 0 ) then
				Uw = sim_getux(xi, yi, zi, Omega, zmin, rmin)
			else if ( dir == 1 ) then
				Uw = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
			else if ( dir == 2 ) then
				Uw = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
			endif
			f15 = Uw
			f25 = Uw
		endif
		if( psi(i, j, k-1) >= 0 ) then
			xi = x
			yi = y
			zi = z - dpsi6
			if( dir == 0 ) then
				Uw = sim_getux(xi, yi, zi, Omega, zmin, rmin)
			else if ( dir == 1 ) then
				Uw = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
			else if ( dir == 2 ) then
				Uw = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
			endif
			f16 = Uw
			f26 = Uw
		endif
		if( psi(i+2, j, k) >= 0 ) then
			f21 = f00
			f21 = f11
		endif
		if( psi(i-2, j, k) >= 0 ) then
			f23 = f00
			f23 = f13
		endif
		if( psi(i, j+2, k) >= 0 ) then
			f22 = f00
			f22 = f12
		endif
		if( psi(i, j-2, k) >= 0 ) then
			f24 = f00
			f24 = f14
		endif
		if( psi(i, j, k+2) >= 0 ) then
			f25 = f00
			f25 = f15
		endif
		if( psi(i, j, k-2) >= 0 ) then
			f26 = f00
			f26 = f16
		endif

    f1 = sim_getupwind(vx1, 0.125*(6.0*f00 + 3.0*f11 - f13), 0.125*(3.0*f00 + 6.0*f11 - f21))
    f3 = sim_getupwind(vx3, 0.125*(3.0*f00 + 6.0*f13 - f23), 0.125*(6.0*f00 + 3.0*f13 - f11))
    f2 = sim_getupwind(vy2, 0.125*(6.0*f00 + 3.0*f12 - f14), 0.125*(3.0*f00 + 6.0*f12 - f22))
    f4 = sim_getupwind(vy4, 0.125*(3.0*f00 + 6.0*f14 - f24), 0.125*(6.0*f00 + 3.0*f14 - f12))
    f5 = sim_getupwind(vz5, 0.125*(6.0*f00 + 3.0*f15 - f16), 0.125*(3.0*f00 + 6.0*f15 - f25))
    f6 = sim_getupwind(vz6, 0.125*(3.0*f00 + 6.0*f16 - f26), 0.125*(6.0*f00 + 3.0*f16 - f15))

    q1 = f1*vx1
    q3 = f3*vx3
    q2 = f2*vy2
    q4 = f4*vy4
    q5 = f5*vz5
    q6 = f6*vz6

    if( psi1 >= 0 .and. psi3 >= 0 ) then
      q1 = 0.0
      q3 = 0.0
    endif
    if( psi2 >= 0 .and. psi4 >= 0 ) then
      q2 = 0.0
      q4 = 0.0
    endif
    if( psi5 >= 0 .and. psi6 >= 0 ) then
      q5 = 0.0
      q6 = 0.0
    endif

    fc(i, j, k) = q1 - q3 + q2 - q4 + q5 - q6

  end do
  end do
  end do
!$omp end do
!$omp end parallel
end subroutine sim_calc_c_u_quick

subroutine sim_calc_c_u_test( &
                fc &
              , f &
              , v &
              , phi &
              , psi &
							, dir &
							, Omega &
							, zmin &
							, zmax &
							, rmin &
							, rmax &
							, org &
              , sz, g)
  implicit none
  integer                  :: i, j, k
  integer                  :: ix, jx, kx
  integer                  :: g
  integer, dimension(3)    :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: fc
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: f
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:6)  :: v
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: phi, psi
	integer									:: dir
	real										:: Omega, zmin, zmax, rmin, rmax
	real, dimension(3)			:: org
  real                    :: f00
  real                    :: f11, f13, f12, f14, f15, f16
  real                    :: f21, f23, f22, f24, f25, f26
  real                    :: f1, f3, f2, f4, f5, f6
  real                    :: vx1, vx3, vy2, vy4, vz5, vz6
  real                    :: q1, q3, q2, q4, q5, q6
  real                    :: psi0, psi1, psi2, psi3, psi4, psi5, psi6
  real                    :: Hpsi0, Hpsi1, Hpsi2, Hpsi3, Hpsi4, Hpsi5, Hpsi6
  real                    :: dpsi1, dpsi2, dpsi3, dpsi4, dpsi5, dpsi6
	real										:: x, y, z, xi, yi, zi
	real										:: Uw
  real                    :: sim_getweno3
  real                    :: sim_geth
  real                    :: sim_getd
	real										:: sim_getux
	real										:: sim_getuy
	real										:: sim_getuz
  real                    :: sim_getupwind
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
!$omp parallel private(i, j, k) &
!$omp          ,private(f00) &
!$omp          ,private(f11, f13, f12, f14, f15, f16) &
!$omp          ,private(f21, f23, f22, f24, f25, f26) &
!$omp          ,private(f1, f3, f2, f4, f5, f6) &
!$omp          ,private(vx1, vx3, vy2, vy4, vz5, vz6) &
!$omp          ,private(q1, q3, q2, q4, q5, q6) &
!$omp          ,private(psi0, psi1, psi2, psi3, psi4, psi5, psi6) &
!$omp          ,private(Hpsi0, Hpsi1, Hpsi2, Hpsi3, Hpsi4, Hpsi5, Hpsi6) &
!$omp          ,private(dpsi1, dpsi2, dpsi3, dpsi4, dpsi5, dpsi6) &
!$omp					 ,private(x, y, z, xi, yi, zi) &
!$omp          ,private(Uw) 
!$omp do schedule(dynamic, 1)
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
		x = org(1) + float(i) - 0.5
		y = org(2) + float(j) - 0.5
		z = org(3) + float(k) - 0.5

    vx1 = v(i, j, k, 1)
    vx3 = v(i, j, k, 3)
    vy2 = v(i, j, k, 2)
    vy4 = v(i, j, k, 4)
    vz5 = v(i, j, k, 5)
    vz6 = v(i, j, k, 6)

    psi0 = psi(i, j, k)
    psi1 = psi(i+1, j, k)
    psi3 = psi(i-1, j, k)
    psi2 = psi(i, j+1, k)
    psi4 = psi(i, j-1, k)
    psi5 = psi(i, j, k+1)
    psi6 = psi(i, j, k-1)

    Hpsi0 = sim_geth(psi0)
    Hpsi1 = sim_geth(psi1)
    Hpsi3 = sim_geth(psi3)
    Hpsi2 = sim_geth(psi2)
    Hpsi4 = sim_geth(psi4)
    Hpsi5 = sim_geth(psi5)
    Hpsi6 = sim_geth(psi6)

    dpsi1 = sim_getd(psi0, psi1)
    dpsi3 = sim_getd(psi0, psi3)
    dpsi2 = sim_getd(psi0, psi2)
    dpsi4 = sim_getd(psi0, psi4)
    dpsi5 = sim_getd(psi0, psi5)
    dpsi6 = sim_getd(psi0, psi6)

    f00 = f(i, j, k)
    f11 = f(i+1, j, k)
    f13 = f(i-1, j, k)
    f12 = f(i, j+1, k)
    f14 = f(i, j-1, k)
    f15 = f(i, j, k+1)
    f16 = f(i, j, k-1)

		f1 = 0.5*(f00 + f11)*0.95
		f3 = 0.5*(f00 + f13)*0.95
		f2 = 0.5*(f00 + f12)*0.95
		f4 = 0.5*(f00 + f14)*0.95
		f5 = 0.5*(f00 + f15)*0.95
		f6 = 0.5*(f00 + f16)*0.95

		f1 = f1 + sim_getupwind(vx1, f00, f11)*0.05
		f3 = f3 + sim_getupwind(vx3, f13, f00)*0.05
		f2 = f2 + sim_getupwind(vy2, f00, f12)*0.05
		f4 = f4 + sim_getupwind(vy4, f14, f00)*0.05
		f5 = f5 + sim_getupwind(vz5, f00, f15)*0.05
		f6 = f6 + sim_getupwind(vz6, f16, f00)*0.05

		if( psi(i+1, j, k) >= 0 ) then
			xi = x + dpsi1
			yi = y
			zi = z
			if( dir == 0 ) then
				Uw = sim_getux(xi, yi, zi, Omega, zmin, rmin)
			else if ( dir == 1 ) then
				Uw = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
			else if ( dir == 2 ) then
				Uw = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
			endif
			f11 = Uw
			f1 = 0.5*(f00 + f11)
		endif
		if( psi(i-1, j, k) >= 0 ) then
			xi = x - dpsi3
			yi = y
			zi = z
			if( dir == 0 ) then
				Uw = sim_getux(xi, yi, zi, Omega, zmin, rmin)
			else if ( dir == 1 ) then
				Uw = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
			else if ( dir == 2 ) then
				Uw = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
			endif
			f13 = Uw
			f3 = 0.5*(f00 + f13)
		endif
		if( psi(i, j+1, k) >= 0 ) then
			xi = x
			yi = y + dpsi2
			zi = z
			if( dir == 0 ) then
				Uw = sim_getux(xi, yi, zi, Omega, zmin, rmin)
			else if ( dir == 1 ) then
				Uw = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
			else if ( dir == 2 ) then
				Uw = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
			endif
			f12 = Uw
			f2 = 0.5*(f00 + f12)
		endif
		if( psi(i, j-1, k) >= 0 ) then
			xi = x
			yi = y - dpsi4
			zi = z
			if( dir == 0 ) then
				Uw = sim_getux(xi, yi, zi, Omega, zmin, rmin)
			else if ( dir == 1 ) then
				Uw = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
			else if ( dir == 2 ) then
				Uw = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
			endif
			f14 = Uw
			f4 = 0.5*(f00 + f14)
		endif
		if( psi(i, j, k+1) >= 0 ) then
			xi = x
			yi = y
			zi = z + dpsi5
			if( dir == 0 ) then
				Uw = sim_getux(xi, yi, zi, Omega, zmin, rmin)
			else if ( dir == 1 ) then
				Uw = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
			else if ( dir == 2 ) then
				Uw = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
			endif
			f15 = Uw
			f5 = 0.5*(f00 + f15)
		endif
		if( psi(i, j, k-1) >= 0 ) then
			xi = x
			yi = y
			zi = z - dpsi6
			if( dir == 0 ) then
				Uw = sim_getux(xi, yi, zi, Omega, zmin, rmin)
			else if ( dir == 1 ) then
				Uw = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
			else if ( dir == 2 ) then
				Uw = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
			endif
			f16 = Uw
			f6 = 0.5*(f00 + f16)
		endif

    q1 = f1*vx1
    q3 = f3*vx3
    q2 = f2*vy2
    q4 = f4*vy4
    q5 = f5*vz5
    q6 = f6*vz6

    if( psi1 >= 0 .and. psi3 >= 0 ) then
      q1 = 0.0
      q3 = 0.0
    endif
    if( psi2 >= 0 .and. psi4 >= 0 ) then
      q2 = 0.0
      q4 = 0.0
    endif
    if( psi5 >= 0 .and. psi6 >= 0 ) then
      q5 = 0.0
      q6 = 0.0
    endif

    fc(i, j, k) = q1 - q3 + q2 - q4 + q5 - q6

  end do
  end do
  end do
!$omp end do
!$omp end parallel
end subroutine sim_calc_c_u_test

subroutine sim_calc_abd_u( &
									A &
								, b &
								, u1_ &
								, u0_ &
								, uc0_, ucp_ &
								, ud0_, udp_ &
								, ux_, uy_, uz_ &
								, phi &
								, psi &
								, rhol, rhog &
								, mul, mug &
								, nl_viscosity &
								, nl_viscosity_d &
								, dir &
								, Omega &
								, zmin &
								, zmax &
								, rmin &
								, rmax &
								, org &
								, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:6)	:: A
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: u1_, u0_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: uc0_, ucp_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: ud0_, udp_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: ux_, uy_, uz_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: psi
	real										:: rhol, rhog
	real										:: mul, mug
	real										:: nl_viscosity
	real										:: nl_viscosity_d
	integer									:: dir
	real										:: Omega, zmin, zmax, rmin, rmax
	real, dimension(3)			:: org
	real										:: phi0, phi1, phi3, phi2, phi4, phi5, phi6
	real										:: Hphi0, Hphi1, Hphi3, Hphi2, Hphi4, Hphi5, Hphi6
	real										:: dphi1, dphi3, dphi2, dphi4, dphi5, dphi6
	real										:: psi0, psi1, psi3, psi2, psi4, psi5, psi6
	real										:: Hpsi0, Hpsi1, Hpsi3, Hpsi2, Hpsi4, Hpsi5, Hpsi6
	real										:: dpsi1, dpsi3, dpsi2, dpsi4, dpsi5, dpsi6
	real										:: mask1, mask3, mask2, mask4, mask5, mask6
	real										:: muf0, muf1, muf3, muf2, muf4, muf5, muf6
	real										:: mu0, mu1, mu3, mu2, mu4, mu5, mu6
	real										:: mu1m, mu3m, mu2m, mu4m, mu5m, mu6m
	real										:: rhof0, rho0
	real										:: a1, a3, a2, a4, a5, a6
	real										:: u0, u1, u3, u2, u4, u5, u6
	real										:: ux0, ux1, ux3, ux2, ux4, ux5, ux6
	real										:: uy0, uy1, uy3, uy2, uy4, uy5, uy6
	real										:: uz0, uz1, uz3, uz2, uz4, uz5, uz6
	real										:: duxdx1, duxdx3, duydx2, duydx4, duzdx5, duzdx6
	real										:: duxdy1, duxdy3, duydy2, duydy4, duzdy5, duzdy6
	real										:: duxdz1, duxdz3, duydz2, duydz4, duzdz5, duzdz6
	real										:: Uw
	real										:: x, y, z, r2, rl, xi, yi, zi, theta
	real										:: alpha
	real										:: sim_geth
	real										:: sim_geth_s
	real										:: sim_geth_2
	real										:: sim_getd
	real										:: sim_getux
	real										:: sim_getuy
	real										:: sim_getuz
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	alpha = 1.0
	alpha = 0.5
!$omp parallel private(i, j, k) &
!$omp					,private(phi0, phi1, phi3, phi2, phi4, phi5, phi6) &
!$omp					,private(Hphi0, Hphi1, Hphi3, Hphi2, Hphi4, Hphi5, Hphi6) &
!$omp					,private(dphi1, dphi3, dphi2, dphi4, dphi5, dphi6) &
!$omp					,private(psi0, psi1, psi3, psi2, psi4, psi5, psi6) &
!$omp					,private(Hpsi0, Hpsi1, Hpsi3, Hpsi2, Hpsi4, Hpsi5, Hpsi6) &
!$omp					,private(dpsi1, dpsi3, dpsi2, dpsi4, dpsi5, dpsi6) &
!$omp					,private(mask1, mask3, mask2, mask4, mask5, mask6) &
!$omp					,private(muf0, muf1, muf3, muf2, muf4, muf5, muf6) &
!$omp					,private(mu0, mu1, mu3, mu2, mu4, mu5, mu6) &
!$omp					,private(mu1m, mu3m, mu2m, mu4m, mu5m, mu6m) &
!$omp					,private(rhof0, rho0) &
!$omp					,private(a1, a3, a2, a4, a5, a6) &
!$omp					,private(u0, u1, u3, u2, u4, u5, u6) &
!$omp					,private(ux0, ux1, ux3, ux2, ux4, ux5, ux6) &
!$omp					,private(uy0, uy1, uy3, uy2, uy4, uy5, uy6) &
!$omp					,private(uz0, uz1, uz3, uz2, uz4, uz5, uz6) &
!$omp					,private(duxdx1, duxdx3, duydx2, duydx4, duzdx5, duzdx6) &
!$omp					,private(duxdy1, duxdy3, duydy2, duydy4, duzdy5, duzdy6) &
!$omp					,private(duxdz1, duxdz3, duydz2, duydz4, duzdz5, duzdz6) &
!$omp					,private(Uw) &
!$omp					,private(x, y, z, r2, rl, xi, yi, zi, theta)
!$omp do schedule(dynamic, 1)
	do k=1, kx
	do j=1, jx
!ocl nouxsimd
	do i=1, ix
		x = org(1) + float(i) - 0.5
		y = org(2) + float(j) - 0.5
		z = org(3) + float(k) - 0.5

		phi0 = phi(i, j, k)
		phi1 = phi(i+1, j, k)
		phi3 = phi(i-1, j, k)
		phi2 = phi(i, j+1, k)
		phi4 = phi(i, j-1, k)
		phi5 = phi(i, j, k+1)
		phi6 = phi(i, j, k-1)

		Hphi0 = sim_geth(phi0)
		Hphi1 = sim_geth(phi1)
		Hphi3 = sim_geth(phi3)
		Hphi2 = sim_geth(phi2)
		Hphi4 = sim_geth(phi4)
		Hphi5 = sim_geth(phi5)
		Hphi6 = sim_geth(phi6)

		dphi1 = sim_getd(phi0, phi1)
		dphi3 = sim_getd(phi0, phi3)
		dphi2 = sim_getd(phi0, phi2)
		dphi4 = sim_getd(phi0, phi4)
		dphi5 = sim_getd(phi0, phi5)
		dphi6 = sim_getd(phi0, phi6)

		mu0 = mug + (mul - mug)*Hphi0
		mu1 = mug + (mul - mug)*Hphi1
		mu3 = mug + (mul - mug)*Hphi3
		mu2 = mug + (mul - mug)*Hphi2
		mu4 = mug + (mul - mug)*Hphi4
		mu5 = mug + (mul - mug)*Hphi5
		mu6 = mug + (mul - mug)*Hphi6

		mu1m = mu0*dphi1 + mu1*(1.0 - dphi1)
		mu3m = mu0*dphi3 + mu3*(1.0 - dphi3)
		mu2m = mu0*dphi2 + mu2*(1.0 - dphi2)
		mu4m = mu0*dphi4 + mu4*(1.0 - dphi4)
		mu5m = mu0*dphi5 + mu5*(1.0 - dphi5)
		mu6m = mu0*dphi6 + mu6*(1.0 - dphi6)

		mu1m = mu0*mu1/(mu0*(1.0 - dphi1) + mu1*dphi1)
		mu3m = mu0*mu3/(mu0*(1.0 - dphi3) + mu3*dphi3)
		mu2m = mu0*mu2/(mu0*(1.0 - dphi2) + mu2*dphi2)
		mu4m = mu0*mu4/(mu0*(1.0 - dphi4) + mu4*dphi4)
		mu5m = mu0*mu5/(mu0*(1.0 - dphi5) + mu5*dphi5)
		mu6m = mu0*mu6/(mu0*(1.0 - dphi6) + mu6*dphi6)

		psi0 = psi(i, j, k)
		psi1 = psi(i+1, j, k)
		psi3 = psi(i-1, j, k)
		psi2 = psi(i, j+1, k)
		psi4 = psi(i, j-1, k)
		psi5 = psi(i, j, k+1)
		psi6 = psi(i, j, k-1)

		Hpsi0 = sim_geth(psi0)
		Hpsi1 = sim_geth(psi1)
		Hpsi3 = sim_geth(psi3)
		Hpsi2 = sim_geth(psi2)
		Hpsi4 = sim_geth(psi4)
		Hpsi5 = sim_geth(psi5)
		Hpsi6 = sim_geth(psi6)

		dpsi1 = sim_getd(psi0, psi1)
		dpsi3 = sim_getd(psi0, psi3)
		dpsi2 = sim_getd(psi0, psi2)
		dpsi4 = sim_getd(psi0, psi4)
		dpsi5 = sim_getd(psi0, psi5)
		dpsi6 = sim_getd(psi0, psi6)

		mask1 = 0.0
		mask3 = 0.0
		mask2 = 0.0
		mask4 = 0.0
		mask5 = 0.0
		mask6 = 0.0
		if( psi(i, j, k) >= 0 ) then
		else
			if( psi(i+1, j, k) >= 0 ) then
				mask1 = 1.0
			endif
			if( psi(i-1, j, k) >= 0 ) then
				mask3 = 1.0
			endif
			if( psi(i, j+1, k) >= 0 ) then
				mask2 = 1.0
			endif
			if( psi(i, j-1, k) >= 0 ) then
				mask4 = 1.0
			endif
			if( psi(i, j, k+1) >= 0 ) then
				mask5 = 1.0
			endif
			if( psi(i, j, k-1) >= 0 ) then
				mask6 = 1.0
			endif
		endif

		if( psi(i, j, k) >= 0 ) then
		else
			if( psi(i+1, j, k) >= 0 ) then
				mu1m = mu0/dpsi1*2.0/(dpsi1 + dpsi3)
				mu3m = mu0/dpsi3*2.0/(dpsi1 + dpsi3)
			endif
			if( psi(i-1, j, k) >= 0 ) then
				mu1m = mu0/dpsi1*2.0/(dpsi1 + dpsi3)
				mu3m = mu0/dpsi3*2.0/(dpsi1 + dpsi3)
			endif
			if( psi(i, j+1, k) >= 0 ) then
				mu2m = mu0/dpsi2*2.0/(dpsi2 + dpsi4)
				mu4m = mu0/dpsi4*2.0/(dpsi2 + dpsi4)
			endif
			if( psi(i, j-1, k) >= 0 ) then
				mu2m = mu0/dpsi2*2.0/(dpsi2 + dpsi4)
				mu4m = mu0/dpsi4*2.0/(dpsi2 + dpsi4)
			endif
			if( psi(i, j, k+1) >= 0 ) then
				mu5m = mu0/dpsi5*2.0/(dpsi5 + dpsi6)
				mu6m = mu0/dpsi6*2.0/(dpsi5 + dpsi6)
			endif
			if( psi(i, j, k-1) >= 0 ) then
				mu5m = mu0/dpsi5*2.0/(dpsi5 + dpsi6)
				mu6m = mu0/dpsi6*2.0/(dpsi5 + dpsi6)
			endif
		endif

		u0 = u0_(i, j, k)
		u1 = u0_(i+1, j, k)
		u3 = u0_(i-1, j, k)
		u2 = u0_(i, j+1, k)
		u4 = u0_(i, j-1, k)
		u5 = u0_(i, j, k+1)
		u6 = u0_(i, j, k-1)
		if( psi(i, j, k) >= 0 ) then
		else
			if( psi(i+1, j, k) >= 0 ) then
				xi = x + dpsi1
				yi = y
				zi = z
				if( dir == 0 ) then
					u1 = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 1 ) then
					u1 = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 2 ) then
					u1 = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
				endif
			endif
			if( psi(i-1, j, k) >= 0 ) then
				xi = x - dpsi3
				yi = y
				zi = z
				if( dir == 0 ) then
					u3 = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 1 ) then
					u3 = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 2 ) then
					u3 = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
				endif
			endif
			if( psi(i, j+1, k) >= 0 ) then
				xi = x
				yi = y + dpsi2
				zi = z
				if( dir == 0 ) then
					u2 = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 1 ) then
					u2 = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 2 ) then
					u2 = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
				endif
			endif
			if( psi(i, j-1, k) >= 0 ) then
				xi = x
				yi = y - dpsi4
				zi = z
				if( dir == 0 ) then
					u4 = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 1 ) then
					u4 = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 2 ) then
					u4 = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
				endif
			endif
			if( psi(i, j, k+1) >= 0 ) then
				xi = x
				yi = y
				zi = z + dpsi5
				if( dir == 0 ) then
					u5 = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 1 ) then
					u5 = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 2 ) then
					u5 = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
				endif
			endif
			if( psi(i, j, k-1) >= 0 ) then
				xi = x
				yi = y
				zi = z - dpsi6
				if( dir == 0 ) then
					u6 = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 1 ) then
					u6 = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 2 ) then
					u6 = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
				endif
			endif
		endif

		Hphi0 = sim_geth(phi0)
		rho0 = rhog + (rhol - rhog)*Hphi0

		a1 = mu1m/(rho0)
		a3 = mu3m/(rho0)
		a2 = mu2m/(rho0)
		a4 = mu4m/(rho0)
		a5 = mu5m/(rho0)
		a6 = mu6m/(rho0)

		ux0 = ux_(i, j, k)
		ux1 = ux_(i+1, j, k)
		ux3 = ux_(i-1, j, k)
		ux2 = ux_(i, j+1, k)
		ux4 = ux_(i, j-1, k)
		ux5 = ux_(i, j, k+1)
		ux6 = ux_(i, j, k-1)

		uy0 = uy_(i, j, k)
		uy1 = uy_(i+1, j, k)
		uy3 = uy_(i-1, j, k)
		uy2 = uy_(i, j+1, k)
		uy4 = uy_(i, j-1, k)
		uy5 = uy_(i, j, k+1)
		uy6 = uy_(i, j, k-1)

		uz0 = uz_(i, j, k)
		uz1 = uz_(i+1, j, k)
		uz3 = uz_(i-1, j, k)
		uz2 = uz_(i, j+1, k)
		uz4 = uz_(i, j-1, k)
		uz5 = uz_(i, j, k+1)
		uz6 = uz_(i, j, k-1)

		duxdx1 = ux_(i+1, j, k) - ux_(i, j, k)
		duxdx3 = ux_(i, j, k) - ux_(i-1, j, k)
		duydx2 = 0.25*(uy_(i+1, j, k) - uy_(i-1, j, k) + uy_(i+1, j+1, k) - uy_(i-1, j+1, k))
		duydx4 = 0.25*(uy_(i+1, j, k) - uy_(i-1, j, k) + uy_(i+1, j-1, k) - uy_(i-1, j-1, k))
		duzdx5 = 0.25*(uz_(i+1, j, k) - uz_(i-1, j, k) + uz_(i+1, j, k+1) - uz_(i-1, j, k+1))
		duzdx6 = 0.25*(uz_(i+1, j, k) - uz_(i-1, j, k) + uz_(i+1, j, k-1) - uz_(i-1, j, k-1))

		duxdy1 = 0.25*(ux_(i, j+1, k) - ux_(i, j-1, k) + ux_(i+1, j+1, k) - ux_(i+1, j-1, k))
		duxdy3 = 0.25*(ux_(i, j+1, k) - ux_(i, j-1, k) + ux_(i-1, j+1, k) - ux_(i-1, j-1, k))
		duydy2 = uy_(i, j+1, k) - uy_(i, j, k)
		duydy4 = uy_(i, j, k) - uy_(i, j-1, k)
		duzdy5 = 0.25*(uz_(i, j+1, k) - uz_(i, j-1, k) + uz_(i, j+1, k+1) - uz_(i, j-1, k+1))
		duzdy6 = 0.25*(uz_(i, j+1, k) - uz_(i, j-1, k) + uz_(i, j+1, k-1) - uz_(i, j-1, k-1))

		duxdz1 = 0.25*(ux_(i, j, k+1) - ux_(i, j, k-1) + ux_(i+1, j, k+1) - ux_(i+1, j, k-1))
		duxdz3 = 0.25*(ux_(i, j, k+1) - ux_(i, j, k-1) + ux_(i-1, j, k+1) - ux_(i-1, j, k-1))
		duydz2 = 0.25*(uy_(i, j, k+1) - uy_(i, j, k-1) + uy_(i, j+1, k+1) - uy_(i, j+1, k-1))
		duydz4 = 0.25*(uy_(i, j, k+1) - uy_(i, j, k-1) + uy_(i, j-1, k+1) - uy_(i, j-1, k-1))
		duzdz5 = uz_(i, j, k+1) - uz_(i, j, k)
		duzdz6 = uz_(i, j, k) - uz_(i, j, k-1)

		Hphi0 = sim_geth(phi0)
		if( nl_viscosity_d > 0 ) then
			Hphi0 = sim_geth_2(phi0, nl_viscosity_d)
		end if
		rho0 = rhog + (rhol - rhog)*Hphi0

		ud0_(i, j, k) = 0.0
		if( dir == 0 ) then
			a1 = (1.0 + nl_viscosity)*mu1m/(rho0)
			a3 = (1.0 + nl_viscosity)*mu3m/(rho0)
			ud0_(i, j, k) = ( &
											+ mu2m*duydx2 - mu4m*duydx4 &
											+ mu5m*duzdx5 - mu6m*duzdx6 &
											)/rho0*nl_viscosity
		else if( dir == 1 ) then
			a2 = (1.0 + nl_viscosity)*mu2m/(rho0)
			a4 = (1.0 + nl_viscosity)*mu4m/(rho0)
			ud0_(i, j, k) = ( &
												mu1m*duxdy1 - mu3m*duxdy3 &
											+ mu5m*duzdy5 - mu6m*duzdy6 &
											)/rho0*nl_viscosity
		else if( dir == 2 ) then
			a5 = (1.0 + nl_viscosity)*mu5m/(rho0)
			a6 = (1.0 + nl_viscosity)*mu6m/(rho0)
			ud0_(i, j, k) = ( &
												mu1m*duxdz1 - mu3m*duxdz3 &
											+ mu2m*duydz2 - mu4m*duydz4 &
											)/rho0*nl_viscosity
		end if

		A(i, j, k, 0) = 1.0 + alpha*(a1 + a3) &
												+ alpha*(a2 + a4) &
												+ alpha*(a5 + a6)
		A(i, j, k, 1) = - alpha*a1*(1.0 - mask1)
		A(i, j, k, 3) = - alpha*a3*(1.0 - mask3)
		A(i, j, k, 2) = - alpha*a2*(1.0 - mask2)
		A(i, j, k, 4) = - alpha*a4*(1.0 - mask4)
		A(i, j, k, 5) = - alpha*a5*(1.0 - mask5)
		A(i, j, k, 6) = - alpha*a6*(1.0 - mask6)
		b(i, j, k) = u0_(i, j, k) &
!									- (1.5*uc0_(i, j, k) - 0.5*ucp_(i, j, k)) &
!									+ (1.5*ud0_(i, j, k) - 0.5*udp_(i, j, k)) &
									- (uc0_(i, j, k)) &
									+ (ud0_(i, j, k)) &
									+ (1.0 - alpha)*( &
										+ a1*(u1 - u0) + a3*(u3 - u0) &
										+ a2*(u2 - u0) + a4*(u4 - u0) &
										+ a5*(u5 - u0) + a6*(u6 - u0) &
									) &
									+ alpha*a1*u1*mask1 &
									+ alpha*a3*u3*mask3 &
									+ alpha*a2*u2*mask2 &
									+ alpha*a4*u4*mask4 &
									+ alpha*a5*u5*mask5 &
									+ alpha*a6*u6*mask6 

		u1_(i, j, k) = u0_(i, j, k)

		if( psi(i, j, k) >= 0 ) then
			A(i, j, k, 0) = 1.0
			A(i, j, k, 1) = 0.0
			A(i, j, k, 3) = 0.0
			A(i, j, k, 2) = 0.0
			A(i, j, k, 4) = 0.0
			A(i, j, k, 5) = 0.0
			A(i, j, k, 6) = 0.0
			b(i, j, k) = 0.0
			u1_(i, j, k) = 0.0
		endif

	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine sim_calc_abd_u

subroutine sim_calc_abd_u2( &
									A &
								, b &
								, u1_ &
								, u0_ &
								, uc0_, ucp_ &
								, ud0_, udp_ &
								, ux_, uy_, uz_ &
								, phi &
								, psi &
								, work0_ &
								, work1_ &
								, work2_ &
								, rhol, rhog &
								, mul, mug &
								, nl_viscosity &
								, nl_viscosity_d &
								, dir &
								, Omega &
								, zmin &
								, zmax &
								, rmin &
								, rmax &
								, org &
								, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:6)	:: A
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: u1_, u0_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: uc0_, ucp_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: ud0_, udp_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: ux_, uy_, uz_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: psi
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: work0_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: work1_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: work2_
	real										:: rhol, rhog
	real										:: mul, mug
	real										:: nl_viscosity
	real										:: nl_viscosity_d
	integer									:: dir
	real										:: Omega, zmin, zmax, rmin, rmax
	real, dimension(3)			:: org
	real										:: phi0, phi1, phi3, phi2, phi4, phi5, phi6
	real										:: Hphi0, Hphi1, Hphi3, Hphi2, Hphi4, Hphi5, Hphi6
	real										:: dphi1, dphi3, dphi2, dphi4, dphi5, dphi6
	real										:: psi0, psi1, psi3, psi2, psi4, psi5, psi6
	real										:: Hpsi0, Hpsi1, Hpsi3, Hpsi2, Hpsi4, Hpsi5, Hpsi6
	real										:: dpsi1, dpsi3, dpsi2, dpsi4, dpsi5, dpsi6
	real										:: mask1, mask3, mask2, mask4, mask5, mask6
	real										:: muf0, muf1, muf3, muf2, muf4, muf5, muf6
	real										:: mu0, mu1, mu3, mu2, mu4, mu5, mu6
	real										:: mu1m, mu3m, mu2m, mu4m, mu5m, mu6m
	real										:: rhof0, rho0
	real										:: a1, a3, a2, a4, a5, a6
	real										:: u0, u1, u3, u2, u4, u5, u6
	real										:: ux0, ux1, ux3, ux2, ux4, ux5, ux6
	real										:: uy0, uy1, uy3, uy2, uy4, uy5, uy6
	real										:: uz0, uz1, uz3, uz2, uz4, uz5, uz6
	real										:: duxdx1, duxdx3, duydx2, duydx4, duzdx5, duzdx6
	real										:: duxdy1, duxdy3, duydy2, duydy4, duzdy5, duzdy6
	real										:: duxdz1, duxdz3, duydz2, duydz4, duzdz5, duzdz6
	real										:: Uw
	real										:: x, y, z, r2, rl, xi, yi, zi, theta
	real										:: alpha
	real										:: sim_geth
	real										:: sim_geth_s
	real										:: sim_geth_2
	real										:: sim_getd
	real										:: sim_getux
	real										:: sim_getuy
	real										:: sim_getuz
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	alpha = 1.0
	alpha = 0.5
!$omp parallel private(i, j, k) &
!$omp					,private(phi0, phi1, phi3, phi2, phi4, phi5, phi6) &
!$omp					,private(Hphi0, Hphi1, Hphi3, Hphi2, Hphi4, Hphi5, Hphi6) &
!$omp					,private(dphi1, dphi3, dphi2, dphi4, dphi5, dphi6) &
!$omp					,private(psi0, psi1, psi3, psi2, psi4, psi5, psi6) &
!$omp					,private(Hpsi0, Hpsi1, Hpsi3, Hpsi2, Hpsi4, Hpsi5, Hpsi6) &
!$omp					,private(dpsi1, dpsi3, dpsi2, dpsi4, dpsi5, dpsi6) &
!$omp					,private(mask1, mask3, mask2, mask4, mask5, mask6) &
!$omp					,private(muf0, muf1, muf3, muf2, muf4, muf5, muf6) &
!$omp					,private(mu0, mu1, mu3, mu2, mu4, mu5, mu6) &
!$omp					,private(mu1m, mu3m, mu2m, mu4m, mu5m, mu6m) &
!$omp					,private(rhof0, rho0) &
!$omp					,private(a1, a3, a2, a4, a5, a6) &
!$omp					,private(u0, u1, u3, u2, u4, u5, u6) &
!$omp					,private(ux0, ux1, ux3, ux2, ux4, ux5, ux6) &
!$omp					,private(uy0, uy1, uy3, uy2, uy4, uy5, uy6) &
!$omp					,private(uz0, uz1, uz3, uz2, uz4, uz5, uz6) &
!$omp					,private(duxdx1, duxdx3, duydx2, duydx4, duzdx5, duzdx6) &
!$omp					,private(duxdy1, duxdy3, duydy2, duydy4, duzdy5, duzdy6) &
!$omp					,private(duxdz1, duxdz3, duydz2, duydz4, duzdz5, duzdz6) &
!$omp					,private(Uw) &
!$omp					,private(x, y, z, r2, rl, xi, yi, zi, theta)

!$omp do schedule(static, 1)
	do k=1, kx
	do j=1, jx
!ocl nouxsimd
	do i=1, ix
		x = org(1) + float(i) - 0.5
		y = org(2) + float(j) - 0.5
		z = org(3) + float(k) - 0.5

		psi0 = psi(i, j, k)
		psi1 = psi(i+1, j, k)
		psi3 = psi(i-1, j, k)
		psi2 = psi(i, j+1, k)
		psi4 = psi(i, j-1, k)
		psi5 = psi(i, j, k+1)
		psi6 = psi(i, j, k-1)

		Hpsi0 = sim_geth(psi0)
		Hpsi1 = sim_geth(psi1)
		Hpsi3 = sim_geth(psi3)
		Hpsi2 = sim_geth(psi2)
		Hpsi4 = sim_geth(psi4)
		Hpsi5 = sim_geth(psi5)
		Hpsi6 = sim_geth(psi6)

		dpsi1 = sim_getd(psi0, psi1)
		dpsi3 = sim_getd(psi0, psi3)
		dpsi2 = sim_getd(psi0, psi2)
		dpsi4 = sim_getd(psi0, psi4)
		dpsi5 = sim_getd(psi0, psi5)
		dpsi6 = sim_getd(psi0, psi6)

		ux0 = ux_(i, j, k)
		ux1 = ux_(i+1, j, k)
		ux3 = ux_(i-1, j, k)
		ux2 = ux_(i, j+1, k)
		ux4 = ux_(i, j-1, k)
		ux5 = ux_(i, j, k+1)
		ux6 = ux_(i, j, k-1)
		uy0 = uy_(i, j, k)
		uy1 = uy_(i+1, j, k)
		uy3 = uy_(i-1, j, k)
		uy2 = uy_(i, j+1, k)
		uy4 = uy_(i, j-1, k)
		uy5 = uy_(i, j, k+1)
		uy6 = uy_(i, j, k-1)
		uz0 = uz_(i, j, k)
		uz1 = uz_(i+1, j, k)
		uz3 = uz_(i-1, j, k)
		uz2 = uz_(i, j+1, k)
		uz4 = uz_(i, j-1, k)
		uz5 = uz_(i, j, k+1)
		uz6 = uz_(i, j, k-1)
		if( psi(i, j, k) >= 0 ) then
		else
			if( psi(i+1, j, k) >= 0 ) then
				xi = x + dpsi1
				yi = y
				zi = z
				ux1 = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				uy1 = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				uz1 = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
			endif
			if( psi(i-1, j, k) >= 0 ) then
				xi = x - dpsi3
				yi = y
				zi = z
				ux3 = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				uy3 = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				uz3 = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
			endif
			if( psi(i, j+1, k) >= 0 ) then
				xi = x
				yi = y + dpsi2
				zi = z
				ux2 = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				uy2 = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				uz2 = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
			endif
			if( psi(i, j-1, k) >= 0 ) then
				xi = x
				yi = y - dpsi4
				zi = z
				ux4 = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				uy4 = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				uz4 = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
			endif
			if( psi(i, j, k+1) >= 0 ) then
				xi = x
				yi = y
				zi = z + dpsi5
				ux5 = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				uy5 = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				uz5 = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
			endif
			if( psi(i, j, k-1) >= 0 ) then
				xi = x
				yi = y
				zi = z - dpsi6
				ux6 = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				uy6 = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				uz6 = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
			endif
		endif

		if( dir == 0 ) then
			work0_(i, j, k)   = (ux0 - ux3)/dpsi3
			work1_(i, j, k)   = (uy0 - uy3)/dpsi3
			work2_(i, j, k)   = (uz0 - uz3)/dpsi3
			work0_(i+1, j, k) = (ux1 - ux0)/dpsi1
			work1_(i+1, j, k) = (uy1 - uy0)/dpsi1
			work2_(i+1, j, k) = (uz1 - uz0)/dpsi1
		else if( dir == 1 ) then
			work0_(i, j, k)   = (ux0 - ux4)/dpsi4
			work1_(i, j, k)   = (uy0 - uy4)/dpsi4
			work2_(i, j, k)   = (uz0 - uz4)/dpsi4
			work0_(i, j+1, k) = (ux2 - ux0)/dpsi2
			work1_(i, j+1, k) = (uy2 - uy0)/dpsi2
			work2_(i, j+1, k) = (uz2 - uz0)/dpsi2
		else if( dir == 2 ) then
			work0_(i, j, k)   = (ux0 - ux6)/dpsi6
			work1_(i, j, k)   = (uy0 - uy6)/dpsi6
			work2_(i, j, k)   = (uz0 - uz6)/dpsi6
			work0_(i, j, k+1) = (ux5 - ux0)/dpsi5
			work1_(i, j, k+1) = (uy5 - uy0)/dpsi5
			work2_(i, j, k+1) = (uz5 - uz0)/dpsi5
		end if
	end do
	end do
	end do

!$omp do schedule(static, 1)
	do k=1, kx
	do j=1, jx
!ocl nouxsimd
	do i=1, ix
		x = org(1) + float(i) - 0.5
		y = org(2) + float(j) - 0.5
		z = org(3) + float(k) - 0.5

		phi0 = phi(i, j, k)
		phi1 = phi(i+1, j, k)
		phi3 = phi(i-1, j, k)
		phi2 = phi(i, j+1, k)
		phi4 = phi(i, j-1, k)
		phi5 = phi(i, j, k+1)
		phi6 = phi(i, j, k-1)

		Hphi0 = sim_geth(phi0)
		Hphi1 = sim_geth(phi1)
		Hphi3 = sim_geth(phi3)
		Hphi2 = sim_geth(phi2)
		Hphi4 = sim_geth(phi4)
		Hphi5 = sim_geth(phi5)
		Hphi6 = sim_geth(phi6)

		dphi1 = sim_getd(phi0, phi1)
		dphi3 = sim_getd(phi0, phi3)
		dphi2 = sim_getd(phi0, phi2)
		dphi4 = sim_getd(phi0, phi4)
		dphi5 = sim_getd(phi0, phi5)
		dphi6 = sim_getd(phi0, phi6)

		mu0 = mug + (mul - mug)*Hphi0
		mu1 = mug + (mul - mug)*Hphi1
		mu3 = mug + (mul - mug)*Hphi3
		mu2 = mug + (mul - mug)*Hphi2
		mu4 = mug + (mul - mug)*Hphi4
		mu5 = mug + (mul - mug)*Hphi5
		mu6 = mug + (mul - mug)*Hphi6

		mu1m = mu0*dphi1 + mu1*(1.0 - dphi1)
		mu3m = mu0*dphi3 + mu3*(1.0 - dphi3)
		mu2m = mu0*dphi2 + mu2*(1.0 - dphi2)
		mu4m = mu0*dphi4 + mu4*(1.0 - dphi4)
		mu5m = mu0*dphi5 + mu5*(1.0 - dphi5)
		mu6m = mu0*dphi6 + mu6*(1.0 - dphi6)

		mu1m = mu0*mu1/(mu0*(1.0 - dphi1) + mu1*dphi1)
		mu3m = mu0*mu3/(mu0*(1.0 - dphi3) + mu3*dphi3)
		mu2m = mu0*mu2/(mu0*(1.0 - dphi2) + mu2*dphi2)
		mu4m = mu0*mu4/(mu0*(1.0 - dphi4) + mu4*dphi4)
		mu5m = mu0*mu5/(mu0*(1.0 - dphi5) + mu5*dphi5)
		mu6m = mu0*mu6/(mu0*(1.0 - dphi6) + mu6*dphi6)

		psi0 = psi(i, j, k)
		psi1 = psi(i+1, j, k)
		psi3 = psi(i-1, j, k)
		psi2 = psi(i, j+1, k)
		psi4 = psi(i, j-1, k)
		psi5 = psi(i, j, k+1)
		psi6 = psi(i, j, k-1)

		Hpsi0 = sim_geth(psi0)
		Hpsi1 = sim_geth(psi1)
		Hpsi3 = sim_geth(psi3)
		Hpsi2 = sim_geth(psi2)
		Hpsi4 = sim_geth(psi4)
		Hpsi5 = sim_geth(psi5)
		Hpsi6 = sim_geth(psi6)

		dpsi1 = sim_getd(psi0, psi1)
		dpsi3 = sim_getd(psi0, psi3)
		dpsi2 = sim_getd(psi0, psi2)
		dpsi4 = sim_getd(psi0, psi4)
		dpsi5 = sim_getd(psi0, psi5)
		dpsi6 = sim_getd(psi0, psi6)

		mask1 = 0.0
		mask3 = 0.0
		mask2 = 0.0
		mask4 = 0.0
		mask5 = 0.0
		mask6 = 0.0
		if( psi(i, j, k) >= 0 ) then
		else
			if( psi(i+1, j, k) >= 0 ) then
				mask1 = 1.0
			endif
			if( psi(i-1, j, k) >= 0 ) then
				mask3 = 1.0
			endif
			if( psi(i, j+1, k) >= 0 ) then
				mask2 = 1.0
			endif
			if( psi(i, j-1, k) >= 0 ) then
				mask4 = 1.0
			endif
			if( psi(i, j, k+1) >= 0 ) then
				mask5 = 1.0
			endif
			if( psi(i, j, k-1) >= 0 ) then
				mask6 = 1.0
			endif
		endif

		if( psi(i, j, k) >= 0 ) then
		else
			if( psi(i+1, j, k) >= 0 ) then
				mu1m = mu0/dpsi1*2.0/(dpsi1 + dpsi3)
				mu3m = mu0/dpsi3*2.0/(dpsi1 + dpsi3)
			endif
			if( psi(i-1, j, k) >= 0 ) then
				mu1m = mu0/dpsi1*2.0/(dpsi1 + dpsi3)
				mu3m = mu0/dpsi3*2.0/(dpsi1 + dpsi3)
			endif
			if( psi(i, j+1, k) >= 0 ) then
				mu2m = mu0/dpsi2*2.0/(dpsi2 + dpsi4)
				mu4m = mu0/dpsi4*2.0/(dpsi2 + dpsi4)
			endif
			if( psi(i, j-1, k) >= 0 ) then
				mu2m = mu0/dpsi2*2.0/(dpsi2 + dpsi4)
				mu4m = mu0/dpsi4*2.0/(dpsi2 + dpsi4)
			endif
			if( psi(i, j, k+1) >= 0 ) then
				mu5m = mu0/dpsi5*2.0/(dpsi5 + dpsi6)
				mu6m = mu0/dpsi6*2.0/(dpsi5 + dpsi6)
			endif
			if( psi(i, j, k-1) >= 0 ) then
				mu5m = mu0/dpsi5*2.0/(dpsi5 + dpsi6)
				mu6m = mu0/dpsi6*2.0/(dpsi5 + dpsi6)
			endif
		endif

		u0 = u0_(i, j, k)
		u1 = u0_(i+1, j, k)
		u3 = u0_(i-1, j, k)
		u2 = u0_(i, j+1, k)
		u4 = u0_(i, j-1, k)
		u5 = u0_(i, j, k+1)
		u6 = u0_(i, j, k-1)
		if( psi(i, j, k) >= 0 ) then
		else
			if( psi(i+1, j, k) >= 0 ) then
				xi = x + dpsi1
				yi = y
				zi = z
				if( dir == 0 ) then
					u1 = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 1 ) then
					u1 = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 2 ) then
					u1 = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
				endif
			endif
			if( psi(i-1, j, k) >= 0 ) then
				xi = x - dpsi3
				yi = y
				zi = z
				if( dir == 0 ) then
					u3 = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 1 ) then
					u3 = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 2 ) then
					u3 = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
				endif
			endif
			if( psi(i, j+1, k) >= 0 ) then
				xi = x
				yi = y + dpsi2
				zi = z
				if( dir == 0 ) then
					u2 = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 1 ) then
					u2 = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 2 ) then
					u2 = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
				endif
			endif
			if( psi(i, j-1, k) >= 0 ) then
				xi = x
				yi = y - dpsi4
				zi = z
				if( dir == 0 ) then
					u4 = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 1 ) then
					u4 = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 2 ) then
					u4 = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
				endif
			endif
			if( psi(i, j, k+1) >= 0 ) then
				xi = x
				yi = y
				zi = z + dpsi5
				if( dir == 0 ) then
					u5 = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 1 ) then
					u5 = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 2 ) then
					u5 = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
				endif
			endif
			if( psi(i, j, k-1) >= 0 ) then
				xi = x
				yi = y
				zi = z - dpsi6
				if( dir == 0 ) then
					u6 = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 1 ) then
					u6 = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				else if( dir == 2 ) then
					u6 = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
				endif
			endif
		endif

		Hphi0 = sim_geth(phi0)
		rho0 = rhog + (rhol - rhog)*Hphi0

		a1 = mu1m/(rho0)
		a3 = mu3m/(rho0)
		a2 = mu2m/(rho0)
		a4 = mu4m/(rho0)
		a5 = mu5m/(rho0)
		a6 = mu6m/(rho0)

		ux0 = ux_(i, j, k)
		ux1 = ux_(i+1, j, k)
		ux3 = ux_(i-1, j, k)
		ux2 = ux_(i, j+1, k)
		ux4 = ux_(i, j-1, k)
		ux5 = ux_(i, j, k+1)
		ux6 = ux_(i, j, k-1)

		uy0 = uy_(i, j, k)
		uy1 = uy_(i+1, j, k)
		uy3 = uy_(i-1, j, k)
		uy2 = uy_(i, j+1, k)
		uy4 = uy_(i, j-1, k)
		uy5 = uy_(i, j, k+1)
		uy6 = uy_(i, j, k-1)

		uz0 = uz_(i, j, k)
		uz1 = uz_(i+1, j, k)
		uz3 = uz_(i-1, j, k)
		uz2 = uz_(i, j+1, k)
		uz4 = uz_(i, j-1, k)
		uz5 = uz_(i, j, k+1)
		uz6 = uz_(i, j, k-1)

		duxdx1 = ux_(i+1, j, k) - ux_(i, j, k)
		duxdx3 = ux_(i, j, k) - ux_(i-1, j, k)
		duydx2 = 0.25*(uy_(i+1, j, k) - uy_(i-1, j, k) + uy_(i+1, j+1, k) - uy_(i-1, j+1, k))
		duydx4 = 0.25*(uy_(i+1, j, k) - uy_(i-1, j, k) + uy_(i+1, j-1, k) - uy_(i-1, j-1, k))
		duzdx5 = 0.25*(uz_(i+1, j, k) - uz_(i-1, j, k) + uz_(i+1, j, k+1) - uz_(i-1, j, k+1))
		duzdx6 = 0.25*(uz_(i+1, j, k) - uz_(i-1, j, k) + uz_(i+1, j, k-1) - uz_(i-1, j, k-1))

		duxdy1 = 0.25*(ux_(i, j+1, k) - ux_(i, j-1, k) + ux_(i+1, j+1, k) - ux_(i+1, j-1, k))
		duxdy3 = 0.25*(ux_(i, j+1, k) - ux_(i, j-1, k) + ux_(i-1, j+1, k) - ux_(i-1, j-1, k))
		duydy2 = uy_(i, j+1, k) - uy_(i, j, k)
		duydy4 = uy_(i, j, k) - uy_(i, j-1, k)
		duzdy5 = 0.25*(uz_(i, j+1, k) - uz_(i, j-1, k) + uz_(i, j+1, k+1) - uz_(i, j-1, k+1))
		duzdy6 = 0.25*(uz_(i, j+1, k) - uz_(i, j-1, k) + uz_(i, j+1, k-1) - uz_(i, j-1, k-1))

		duxdz1 = 0.25*(ux_(i, j, k+1) - ux_(i, j, k-1) + ux_(i+1, j, k+1) - ux_(i+1, j, k-1))
		duxdz3 = 0.25*(ux_(i, j, k+1) - ux_(i, j, k-1) + ux_(i-1, j, k+1) - ux_(i-1, j, k-1))
		duydz2 = 0.25*(uy_(i, j, k+1) - uy_(i, j, k-1) + uy_(i, j+1, k+1) - uy_(i, j+1, k-1))
		duydz4 = 0.25*(uy_(i, j, k+1) - uy_(i, j, k-1) + uy_(i, j-1, k+1) - uy_(i, j-1, k-1))
		duzdz5 = uz_(i, j, k+1) - uz_(i, j, k)
		duzdz6 = uz_(i, j, k) - uz_(i, j, k-1)

		Hphi0 = sim_geth(phi0)
		if( nl_viscosity_d > 0 ) then
			Hphi0 = sim_geth_2(phi0, nl_viscosity_d)
		end if
		rho0 = rhog + (rhol - rhog)*Hphi0

		if( dir == 0 ) then
			a1 = (1.0 + nl_viscosity)*mu1m/(rho0)
			a3 = (1.0 + nl_viscosity)*mu3m/(rho0)
!			duydx2 = 0.25*( work0_(i, j, k) + work0_(i+1, j, k) + work0_(i, j+1, k) + work0_(i+1, j+1, k) )
!			duydx4 = 0.25*( work0_(i, j, k) + work0_(i+1, j, k) + work0_(i, j-1, k) + work0_(i+1, j-1, k) )
!			duzdx5 = 0.25*( work1_(i, j, k) + work1_(i+1, j, k) + work1_(i, j, k+1) + work1_(i+1, j, k+1) )
!			duzdx6 = 0.25*( work1_(i, j, k) + work1_(i+1, j, k) + work1_(i, j, k-1) + work1_(i+1, j, k-1) )
			ud0_(i, j, k) = ( &
											+ mu2m*duydx2 - mu4m*duydx4 &
											+ mu5m*duzdx5 - mu6m*duzdx6 &
											)/rho0*nl_viscosity
		else if( dir == 1 ) then
			a2 = (1.0 + nl_viscosity)*mu2m/(rho0)
			a4 = (1.0 + nl_viscosity)*mu4m/(rho0)
!			duxdy1 = 0.25*( work0_(i, j, k) + work0_(i, j+1, k) + work0_(i+1, j, k) + work0_(i+1, j+1, k) )
!			duxdy3 = 0.25*( work0_(i, j, k) + work0_(i, j+1, k) + work0_(i-1, j, k) + work0_(i-1, j+1, k) )
!			duzdy5 = 0.25*( work1_(i, j, k) + work1_(i, j+1, k) + work1_(i, j, k+1) + work1_(i, j+1, k+1) )
!			duzdy6 = 0.25*( work1_(i, j, k) + work1_(i, j+1, k) + work1_(i, j, k-1) + work1_(i, j+1, k-1) )
			ud0_(i, j, k) = ( &
												mu1m*duxdy1 - mu3m*duxdy3 &
											+ mu5m*duzdy5 - mu6m*duzdy6 &
											)/rho0*nl_viscosity
		else if( dir == 2 ) then
			a5 = (1.0 + nl_viscosity)*mu5m/(rho0)
			a6 = (1.0 + nl_viscosity)*mu6m/(rho0)
!			duxdz1 = 0.25*( work0_(i, j, k) + work0_(i, j, k+1) + work0_(i+1, j, k) + work0_(i+1, j, k+1) )
!			duxdz3 = 0.25*( work0_(i, j, k) + work0_(i, j, k+1) + work0_(i-1, j, k) + work0_(i-1, j, k+1) )
!			duydz2 = 0.25*( work1_(i, j, k) + work1_(i, j, k+1) + work1_(i, j+1, k) + work1_(i, j+1, k+1) )
!			duydz4 = 0.25*( work1_(i, j, k) + work1_(i, j, k+1) + work1_(i, j-1, k) + work1_(i, j-1, k+1) )
!			duzdz5 = work2_(i, j, k+1)
!			duzdz6 = work2_(i, j, k)
			ud0_(i, j, k) = ( &
												mu1m*duxdz1 - mu3m*duxdz3 &
											+ mu2m*duydz2 - mu4m*duydz4 &
											)/rho0*nl_viscosity
		end if

		A(i, j, k, 0) = 1.0 + alpha*(a1 + a3) &
												+ alpha*(a2 + a4) &
												+ alpha*(a5 + a6)
		A(i, j, k, 1) = - alpha*a1*(1.0 - mask1)
		A(i, j, k, 3) = - alpha*a3*(1.0 - mask3)
		A(i, j, k, 2) = - alpha*a2*(1.0 - mask2)
		A(i, j, k, 4) = - alpha*a4*(1.0 - mask4)
		A(i, j, k, 5) = - alpha*a5*(1.0 - mask5)
		A(i, j, k, 6) = - alpha*a6*(1.0 - mask6)
		b(i, j, k) = u0_(i, j, k) &
!									- (1.5*uc0_(i, j, k) - 0.5*ucp_(i, j, k)) &
!									+ (1.5*ud0_(i, j, k) - 0.5*udp_(i, j, k)) &
									- (uc0_(i, j, k)) &
									+ (ud0_(i, j, k)) &
									+ (1.0 - alpha)*( &
										+ a1*(u1 - u0) + a3*(u3 - u0) &
										+ a2*(u2 - u0) + a4*(u4 - u0) &
										+ a5*(u5 - u0) + a6*(u6 - u0) &
									) &
									+ alpha*a1*u1*mask1 &
									+ alpha*a3*u3*mask3 &
									+ alpha*a2*u2*mask2 &
									+ alpha*a4*u4*mask4 &
									+ alpha*a5*u5*mask5 &
									+ alpha*a6*u6*mask6 

		u1_(i, j, k) = u0_(i, j, k)

		if( psi(i, j, k) >= 0 ) then
			A(i, j, k, 0) = 1.0
			A(i, j, k, 1) = 0.0
			A(i, j, k, 3) = 0.0
			A(i, j, k, 2) = 0.0
			A(i, j, k, 4) = 0.0
			A(i, j, k, 5) = 0.0
			A(i, j, k, 6) = 0.0
			b(i, j, k) = 0.0
			u1_(i, j, k) = 0.0
		endif

	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine sim_calc_abd_u2

subroutine sim_calc_abd_t( &
									A &
								, b &
								, t1_ &
								, t0_ &
								, tc0_, tcp_ &
								, td0_, tdp_ &
								, phi &
								, psi &
								, rhol, rhog, rhos &
								, cpl, cpg, cps &
								, kl, kg, ks &
								, org &
								, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:6)	:: A
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: t1_, t0_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: tc0_, tcp_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: td0_, tdp_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: psi
	real										:: rhol, rhog, rhos
	real										:: cpl, cpg, cps
	real										:: kl, kg, ks
	real, dimension(3)			:: org
	real										:: phi0, phi1, phi3, phi2, phi4, phi5, phi6
	real										:: Hphi0, Hphi1, Hphi3, Hphi2, Hphi4, Hphi5, Hphi6
	real										:: dphi1, dphi3, dphi2, dphi4, dphi5, dphi6
	real										:: psi0, psi1, psi3, psi2, psi4, psi5, psi6
	real										:: Hpsi0, Hpsi1, Hpsi3, Hpsi2, Hpsi4, Hpsi5, Hpsi6
	real										:: dpsi1, dpsi3, dpsi2, dpsi4, dpsi5, dpsi6
	real										:: mask1, mask3, mask2, mask4, mask5, mask6
	real										:: kf0, kf1, kf3, kf2, kf4, kf5, kf6
	real										:: k0, k1, k3, k2, k4, k5, k6
	real										:: k1m, k3m, k2m, k4m, k5m, k6m
	real										:: rhof0, rho0
	real										:: cpf0, cp0
	real										:: a1, a3, a2, a4, a5, a6
	real										:: t0, t1, t3, t2, t4, t5, t6
	real										:: Tw
	real										:: x, y, z, r2, rl, xi, yi, zi, theta
	real										:: alpha
	real										:: sim_geth
	real										:: sim_geth_s
	real										:: sim_getd
	real										:: sim_getux
	real										:: sim_getuy
	real										:: sim_getuz
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	alpha = 0.5
!$omp parallel private(i, j, k) &
!$omp					,private(phi0, phi1, phi3, phi2, phi4, phi5, phi6) &
!$omp					,private(Hphi0, Hphi1, Hphi3, Hphi2, Hphi4, Hphi5, Hphi6) &
!$omp					,private(dphi1, dphi3, dphi2, dphi4, dphi5, dphi6) &
!$omp					,private(psi0, psi1, psi3, psi2, psi4, psi5, psi6) &
!$omp					,private(Hpsi0, Hpsi1, Hpsi3, Hpsi2, Hpsi4, Hpsi5, Hpsi6) &
!$omp					,private(dpsi1, dpsi3, dpsi2, dpsi4, dpsi5, dpsi6) &
!$omp					,private(mask1, mask3, mask2, mask4, mask5, mask6) &
!$omp					,private(kf0, kf1, kf3, kf2, kf4, kf5, kf6) &
!$omp					,private(k0, k1, k3, k2, k4, k5, k6) &
!$omp					,private(k1m, k3m, k2m, k4m, k5m, k6m) &
!$omp					,private(rhof0, rho0) &
!$omp					,private(cpf0, cp0) &
!$omp					,private(a1, a3, a2, a4, a5, a6) &
!$omp					,private(t0, t1, t3, t2, t4, t5, t6) &
!$omp					,private(Tw) &
!$omp					,private(x, y, z, r2, rl, xi, yi, zi, theta)
!$omp do schedule(dynamic, 1)
	do k=1, kx
	do j=1, jx
!ocl nouxsimd
	do i=1, ix
		x = org(1) + float(i) - 0.5
		y = org(2) + float(j) - 0.5
		z = org(3) + float(k) - 0.5

		phi0 = phi(i, j, k)
		phi1 = phi(i+1, j, k)
		phi3 = phi(i-1, j, k)
		phi2 = phi(i, j+1, k)
		phi4 = phi(i, j-1, k)
		phi5 = phi(i, j, k+1)
		phi6 = phi(i, j, k-1)

		Hphi0 = sim_geth(phi0)
		Hphi1 = sim_geth(phi1)
		Hphi3 = sim_geth(phi3)
		Hphi2 = sim_geth(phi2)
		Hphi4 = sim_geth(phi4)
		Hphi5 = sim_geth(phi5)
		Hphi6 = sim_geth(phi6)

		dphi1 = sim_getd(phi0, phi1)
		dphi3 = sim_getd(phi0, phi3)
		dphi2 = sim_getd(phi0, phi2)
		dphi4 = sim_getd(phi0, phi4)
		dphi5 = sim_getd(phi0, phi5)
		dphi6 = sim_getd(phi0, phi6)

		kf0 = kg + (kl - kg)*Hphi0
		kf1 = kg + (kl - kg)*Hphi1
		kf3 = kg + (kl - kg)*Hphi3
		kf2 = kg + (kl - kg)*Hphi2
		kf4 = kg + (kl - kg)*Hphi4
		kf5 = kg + (kl - kg)*Hphi5
		kf6 = kg + (kl - kg)*Hphi6

		k0 = kg + (kl - kg)*Hphi0
		k1 = kg + (kl - kg)*Hphi1
		k3 = kg + (kl - kg)*Hphi3
		k2 = kg + (kl - kg)*Hphi2
		k4 = kg + (kl - kg)*Hphi4
		k5 = kg + (kl - kg)*Hphi5
		k6 = kg + (kl - kg)*Hphi6

		k1m = kf0*kf1/(kf0*(1.0 - dphi1) + kf1*dphi1)
		k3m = kf0*kf3/(kf0*(1.0 - dphi3) + kf3*dphi3)
		k2m = kf0*kf2/(kf0*(1.0 - dphi2) + kf2*dphi2)
		k4m = kf0*kf4/(kf0*(1.0 - dphi4) + kf4*dphi4)
		k5m = kf0*kf5/(kf0*(1.0 - dphi5) + kf5*dphi5)
		k6m = kf0*kf6/(kf0*(1.0 - dphi6) + kf6*dphi6)

		psi0 = psi(i, j, k)
		psi1 = psi(i+1, j, k)
		psi3 = psi(i-1, j, k)
		psi2 = psi(i, j+1, k)
		psi4 = psi(i, j-1, k)
		psi5 = psi(i, j, k+1)
		psi6 = psi(i, j, k-1)

		Hpsi0 = sim_geth(psi0)
		Hpsi1 = sim_geth(psi1)
		Hpsi3 = sim_geth(psi3)
		Hpsi2 = sim_geth(psi2)
		Hpsi4 = sim_geth(psi4)
		Hpsi5 = sim_geth(psi5)
		Hpsi6 = sim_geth(psi6)

		dpsi1 = sim_getd(psi0, psi1)
		dpsi3 = sim_getd(psi0, psi3)
		dpsi2 = sim_getd(psi0, psi2)
		dpsi4 = sim_getd(psi0, psi4)
		dpsi5 = sim_getd(psi0, psi5)
		dpsi6 = sim_getd(psi0, psi6)

		mask1 = 0.0
		mask3 = 0.0
		mask2 = 0.0
		mask4 = 0.0
		mask5 = 0.0
		mask6 = 0.0
		if( psi(i, j, k) >= 0 ) then
		else
			if( psi(i+1, j, k) >= 0 ) then
				mask1 = 1.0
			endif
			if( psi(i-1, j, k) >= 0 ) then
				mask3 = 1.0
			endif
			if( psi(i, j+1, k) >= 0 ) then
				mask2 = 1.0
			endif
			if( psi(i, j-1, k) >= 0 ) then
				mask4 = 1.0
			endif
			if( psi(i, j, k+1) >= 0 ) then
				mask5 = 1.0
			endif
			if( psi(i, j, k-1) >= 0 ) then
				mask6 = 1.0
			endif
		endif

		if( psi(i, j, k) >= 0 ) then
		else
			if( psi(i+1, j, k) >= 0 ) then
				k1m = k0/dpsi1*2.0/(dpsi1 + dpsi3)
				k3m = k0/dpsi3*2.0/(dpsi1 + dpsi3)
			endif
			if( psi(i-1, j, k) >= 0 ) then
				k1m = k0/dpsi1*2.0/(dpsi1 + dpsi3)
				k3m = k0/dpsi3*2.0/(dpsi1 + dpsi3)
			endif
			if( psi(i, j+1, k) >= 0 ) then
				k2m = k0/dpsi2*2.0/(dpsi2 + dpsi4)
				k4m = k0/dpsi4*2.0/(dpsi2 + dpsi4)
			endif
			if( psi(i, j-1, k) >= 0 ) then
				k2m = k0/dpsi2*2.0/(dpsi2 + dpsi4)
				k4m = k0/dpsi4*2.0/(dpsi2 + dpsi4)
			endif
			if( psi(i, j, k+1) >= 0 ) then
				k5m = k0/dpsi5*2.0/(dpsi5 + dpsi6)
				k6m = k0/dpsi6*2.0/(dpsi5 + dpsi6)
			endif
			if( psi(i, j, k-1) >= 0 ) then
				k5m = k0/dpsi5*2.0/(dpsi5 + dpsi6)
				k6m = k0/dpsi6*2.0/(dpsi5 + dpsi6)
			endif
		endif

		t0 = t0_(i, j, k)
		t1 = t0_(i+1, j, k)
		t3 = t0_(i-1, j, k)
		t2 = t0_(i, j+1, k)
		t4 = t0_(i, j-1, k)
		t5 = t0_(i, j, k+1)
		t6 = t0_(i, j, k-1)
		if( psi(i, j, k) >= 0 ) then
		else
			if( psi(i+1, j, k) >= 0 ) then
				xi = x + dpsi1
				yi = y
				zi = z
				t1 = 0.0
			endif
			if( psi(i-1, j, k) >= 0 ) then
				xi = x - dpsi3
				yi = y
				zi = z
				t3 = 0.0
			endif
			if( psi(i, j+1, k) >= 0 ) then
				xi = x
				yi = y + dpsi2
				zi = z
				t2 = 0.0
			endif
			if( psi(i, j-1, k) >= 0 ) then
				xi = x
				yi = y - dpsi4
				zi = z
				t4 = 0.0
			endif
			if( psi(i, j, k+1) >= 0 ) then
				xi = x
				yi = y
				zi = z + dpsi5
				t5 = 0.0
			endif
			if( psi(i, j, k-1) >= 0 ) then
				xi = x
				yi = y
				zi = z - dpsi6
				t6 = 0.0
			endif
		endif

		rho0 = rhog + (rhol - rhog)*Hphi0
		cp0 = cpg + (cpl - cpg)*Hphi0

		a1 = k1m/(rho0*cp0)
		a3 = k3m/(rho0*cp0)
		a2 = k2m/(rho0*cp0)
		a4 = k4m/(rho0*cp0)
		a5 = k5m/(rho0*cp0)
		a6 = k6m/(rho0*cp0)

		td0_(i, j, k) = 0.0

		A(i, j, k, 0) = 1.0 + alpha*(a1 + a3) &
												+ alpha*(a2 + a4) &
												+ alpha*(a5 + a6)
		A(i, j, k, 1) = - alpha*a1*(1.0 - mask1)
		A(i, j, k, 3) = - alpha*a3*(1.0 - mask3)
		A(i, j, k, 2) = - alpha*a2*(1.0 - mask2)
		A(i, j, k, 4) = - alpha*a4*(1.0 - mask4)
		A(i, j, k, 5) = - alpha*a5*(1.0 - mask5)
		A(i, j, k, 6) = - alpha*a6*(1.0 - mask6)
		b(i, j, k) = t0_(i, j, k) &
!									- (1.5*tc0_(i, j, k) - 0.5*tcp_(i, j, k)) &
!									+ (1.5*td0_(i, j, k) - 0.5*tdp_(i, j, k)) &
									- tc0_(i, j, k) &
									+ (1.0 - alpha)*( &
										+ a1*(t1 - t0) + a3*(t3 - t0) &
										+ a2*(t2 - t0) + a4*(t4 - t0) &
										+ a5*(t5 - t0) + a6*(t6 - t0) &
									) &
									+ td0_(i, j, k) &
									+ alpha*a1*t1*mask1 &
									+ alpha*a3*t3*mask3 &
									+ alpha*a2*t2*mask2 &
									+ alpha*a4*t4*mask4 &
									+ alpha*a5*t5*mask5 &
									+ alpha*a6*t6*mask6 

		t1_(i, j, k) = t0_(i, j, k)

		if( psi(i, j, k) >= 0 ) then
			A(i, j, k, 0) = 1.0
			A(i, j, k, 1) = 0.0
			A(i, j, k, 3) = 0.0
			A(i, j, k, 2) = 0.0
			A(i, j, k, 4) = 0.0
			A(i, j, k, 5) = 0.0
			A(i, j, k, 6) = 0.0
			b(i, j, k) = 0.0
			t1_(i, j, k) = 0.0
		endif

	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine sim_calc_abd_t

subroutine sim_add_g(ux, uy, uz, phi, psi, gx, gy, gz, OmegaCoriolis, OmegaCentrifugal, org, sz, g)
  implicit none
  integer                  :: i, j, k
  integer                  :: ix, jx, kx
  integer                  :: g
  integer, dimension(3)    :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: ux, uy, uz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: phi, psi
  real                    :: gx, gy, gz
	real										:: OmegaCoriolis
	real										:: OmegaCentrifugal
	real, dimension(3)			:: org
	real										:: x, y, z, r2, rl
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(x, y, z, r2, rl)
!$omp do
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
		x = org(1) + float(i) - 0.5
		y = org(2) + float(j) - 0.5
		z = org(3) + float(k) - 0.5
    ux(i, j, k) = ux(i, j, k) + gx + 2.0*OmegaCoriolis*uy(i, j, k) + OmegaCentrifugal*OmegaCentrifugal*x
    uy(i, j, k) = uy(i, j, k) + gy - 2.0*OmegaCoriolis*ux(i, j, k) + OmegaCentrifugal*OmegaCentrifugal*y
    uz(i, j, k) = uz(i, j, k) + gz
    if( psi(i, j, k) >= 0 ) then
      ux(i, j, k) = 0.0
      uy(i, j, k) = 0.0
      uz(i, j, k) = 0.0
    endif
  end do
  end do
  end do
!$omp end do
!$omp end parallel
end subroutine sim_add_g

subroutine sim_add_g_nc(ux, uy, uz, t, phi, psi, gx, gy, gz, betag, tr, OmegaCoriolis, OmegaCentrifugal, org, sz, g)
  implicit none
  integer                  :: i, j, k
  integer                  :: ix, jx, kx
  integer                  :: g
  integer, dimension(3)    :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: ux, uy, uz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: t
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: phi, psi
  real                    :: gx, gy, gz
  real                    :: betag, tr
	real										:: OmegaCoriolis
	real										:: OmegaCentrifugal
	real, dimension(3)			:: org
	real										:: x, y, z, r2, rl
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(x, y, z, r2, rl)
!$omp do
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
		x = org(1) + float(i) - 0.5
		y = org(2) + float(j) - 0.5
		z = org(3) + float(k) - 0.5
    ux(i, j, k) = ux(i, j, k) + gx + 2.0*OmegaCoriolis*uy(i, j, k) + OmegaCentrifugal*OmegaCentrifugal*x
    uy(i, j, k) = uy(i, j, k) + gy - 2.0*OmegaCoriolis*ux(i, j, k) + OmegaCentrifugal*OmegaCentrifugal*y
    uz(i, j, k) = uz(i, j, k) + gz + betag*(t(i, j, k) - tr)
    if( psi(i, j, k) >= 0 ) then
      ux(i, j, k) = 0.0
      uy(i, j, k) = 0.0
      uz(i, j, k) = 0.0
    endif
  end do
  end do
  end do
!$omp end do
!$omp end parallel
end subroutine sim_add_g_nc

subroutine sim_calc_ab_p( &
                        A &
                      , b &
                      ,  p1_ &
                      , v &
                      , p0_ &
                      , ux, uy, uz &
                      , phi &
                      , psi &
                      , kappa &
                      , rhol, rhog &
                      , sigma &
											, Omega &
											, zmin &
											, zmax &
											, rmin &
											, rmax &
											, org &
                      , sz, g)
  implicit none
  integer                  :: i, j, k
  integer                  :: ix, jx, kx
  integer                  :: g
  integer, dimension(3)    :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:6)  :: A
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: b
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: p1_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:6)  :: v
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: p0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: ux, uy, uz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: phi
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: psi
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: kappa
  real                    :: rhol, rhog
  real                    :: sigma
	real										:: Omega, zmin, zmax, rmin, rmax
	real, dimension(3)			:: org
  real                    :: ux0, uy0, uz0
  real                    :: ux1, ux3, uy2, uy4, uz5, uz6
  real                    :: vx1, vx3, vy2, vy4, vz5, vz6
  real                    :: phi0, phi1, phi3, phi2, phi4, phi5, phi6
  real                    :: Hphi0, Hphi1, Hphi3, Hphi2, Hphi4, Hphi5, Hphi6
  real                    :: dphi1, dphi3, dphi2, dphi4, dphi5, dphi6
  real                    :: psi0, psi1, psi3, psi2, psi4, psi5, psi6
  real                    :: Hpsi0, Hpsi1, Hpsi3, Hpsi2, Hpsi4, Hpsi5, Hpsi6
  real                    :: dpsi1, dpsi3, dpsi2, dpsi4, dpsi5, dpsi6
	real										:: mask1, mask3, mask2, mask4, mask5, mask6
  real                    :: rho0, rho1, rho3, rho2, rho4, rho5, rho6
  real                    :: rho1m, rho3m, rho2m, rho4m, rho5m, rho6m
  real                    :: c1, c3, c2, c4, c5, c6
  real                    :: kappa0, kappa1, kappa3, kappa2, kappa4, kappa5, kappa6
  real                    :: kappa1m, kappa3m, kappa2m, kappa4m, kappa5m, kappa6m
  real                    :: dp1, dp3, dp2, dp4, dp5, dp6
  real                    :: divv
	real										:: Uw, Vw, Ww
	real										:: x, y, z, r2, rl, xi, yi, zi, theta
  real                    :: sim_geth
  real                    :: sim_getd
	real										:: sim_getux
	real										:: sim_getuy
	real										:: sim_getuz
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
!$omp parallel private(i, j, k) &
!$omp          ,private(ux0, uy0, uz0) &
!$omp          ,private(ux1, ux3, uy2, uy4, uz5, uz6) &
!$omp          ,private(vx1, vx3, vy2, vy4, vz5, vz6) &
!$omp          ,private(phi0, phi1, phi3, phi2, phi4, phi5, phi6) &
!$omp          ,private(Hphi0, Hphi1, Hphi3, Hphi2, Hphi4, Hphi5, Hphi6) &
!$omp          ,private(dphi1, dphi3, dphi2, dphi4, dphi5, dphi6) &
!$omp          ,private(psi0, psi1, psi3, psi2, psi4, psi5, psi6) &
!$omp          ,private(Hpsi0, Hpsi1, Hpsi3, Hpsi2, Hpsi4, Hpsi5, Hpsi6) &
!$omp          ,private(dpsi1, dpsi3, dpsi2, dpsi4, dpsi5, dpsi6) &
!$omp					 ,private(mask1, mask3, mask2, mask4, mask5, mask6) &
!$omp          ,private(rho0, rho1, rho3, rho2, rho4, rho5, rho6) &
!$omp          ,private(rho1m, rho3m, rho2m, rho4m, rho5m, rho6m) &
!$omp          ,private(c1, c3, c2, c4, c5, c6) &
!$omp          ,private(kappa0, kappa1, kappa3, kappa2, kappa4, kappa5, kappa6) &
!$omp          ,private(kappa1m, kappa3m, kappa2m, kappa4m, kappa5m, kappa6m) &
!$omp          ,private(dp1, dp3, dp2, dp4, dp5, dp6) &
!$omp          ,private(divv) &
!$omp          ,private(Uw, Vw, Ww) &
!$omp					 ,private(x, y, z, r2, rl, xi, yi, zi, theta)
!$omp do schedule(dynamic, 1)
  do k=1, kx
  do j=1, jx
  do i=1, ix
		x = org(1) + float(i) - 0.5
		y = org(2) + float(j) - 0.5
		z = org(3) + float(k) - 0.5

    phi0 = phi(i, j, k)
    phi1 = phi(i+1, j, k)
    phi3 = phi(i-1, j, k)
    phi2 = phi(i, j+1, k)
    phi4 = phi(i, j-1, k)
    phi5 = phi(i, j, k+1)
    phi6 = phi(i, j, k-1)

    Hphi0 = sim_geth(phi0)
    Hphi1 = sim_geth(phi1)
    Hphi3 = sim_geth(phi3)
    Hphi2 = sim_geth(phi2)
    Hphi4 = sim_geth(phi4)
    Hphi5 = sim_geth(phi5)
    Hphi6 = sim_geth(phi6)

    dphi1 = sim_getd(phi0, phi1)
    dphi3 = sim_getd(phi0, phi3)
    dphi2 = sim_getd(phi0, phi2)
    dphi4 = sim_getd(phi0, phi4)
    dphi5 = sim_getd(phi0, phi5)
    dphi6 = sim_getd(phi0, phi6)

    rho0 = rhog + (rhol - rhog)*Hphi0
    rho1 = rhog + (rhol - rhog)*Hphi1
    rho3 = rhog + (rhol - rhog)*Hphi3
    rho2 = rhog + (rhol - rhog)*Hphi2
    rho4 = rhog + (rhol - rhog)*Hphi4
    rho5 = rhog + (rhol - rhog)*Hphi5
    rho6 = rhog + (rhol - rhog)*Hphi6

    rho1m = rho0*dphi1 + rho1*(1.0 - dphi1)
    rho3m = rho0*dphi3 + rho3*(1.0 - dphi3)
    rho2m = rho0*dphi2 + rho2*(1.0 - dphi2)
    rho4m = rho0*dphi4 + rho4*(1.0 - dphi4)
    rho5m = rho0*dphi5 + rho5*(1.0 - dphi5)
    rho6m = rho0*dphi6 + rho6*(1.0 - dphi6)

    c1 = 1.0/rho1m
    c3 = 1.0/rho3m
    c2 = 1.0/rho2m
    c4 = 1.0/rho4m
    c5 = 1.0/rho5m
    c6 = 1.0/rho6m

    kappa0 = kappa(i, j, k)
    kappa1 = kappa(i+1, j, k)
    kappa3 = kappa(i-1, j, k)
    kappa2 = kappa(i, j+1, k)
    kappa4 = kappa(i, j-1, k)
    kappa5 = kappa(i, j, k+1)
    kappa6 = kappa(i, j, k-1)

    kappa1m = kappa0*(1.0 - dphi1) + kappa1*dphi1
    kappa3m = kappa0*(1.0 - dphi3) + kappa3*dphi3
    kappa2m = kappa0*(1.0 - dphi2) + kappa2*dphi2
    kappa4m = kappa0*(1.0 - dphi4) + kappa4*dphi4
    kappa5m = kappa0*(1.0 - dphi5) + kappa5*dphi5
    kappa6m = kappa0*(1.0 - dphi6) + kappa6*dphi6

    dp1 = -sigma*kappa1m
    dp3 = -sigma*kappa3m
    dp2 = -sigma*kappa2m
    dp4 = -sigma*kappa4m
    dp5 = -sigma*kappa5m
    dp6 = -sigma*kappa6m

    psi0 = psi(i, j, k)
    psi1 = psi(i+1, j, k)
    psi3 = psi(i-1, j, k)
    psi2 = psi(i, j+1, k)
    psi4 = psi(i, j-1, k)
    psi5 = psi(i, j, k+1)
    psi6 = psi(i, j, k-1)

    dpsi1 = sim_getd(psi0, psi1)
    dpsi3 = sim_getd(psi0, psi3)
    dpsi2 = sim_getd(psi0, psi2)
    dpsi4 = sim_getd(psi0, psi4)
    dpsi5 = sim_getd(psi0, psi5)
    dpsi6 = sim_getd(psi0, psi6)

		mask1 = 0.0
		mask3 = 0.0
		mask2 = 0.0
		mask4 = 0.0
		mask5 = 0.0
		mask6 = 0.0
		if( psi(i, j, k) >= 0 ) then
		else
			if( psi(i+1, j, k) >= 0 ) then
				mask1 = 1.0
			endif
			if( psi(i-1, j, k) >= 0 ) then
				mask3 = 1.0
			endif
			if( psi(i, j+1, k) >= 0 ) then
				mask2 = 1.0
			endif
			if( psi(i, j-1, k) >= 0 ) then
				mask4 = 1.0
			endif
			if( psi(i, j, k+1) >= 0 ) then
				mask5 = 1.0
			endif
			if( psi(i, j, k-1) >= 0 ) then
				mask6 = 1.0
			endif
		endif

    ux0 = ux(i, j, k)
    uy0 = uy(i, j, k)
    uz0 = uz(i, j, k)

    ux1 = ux(i+1, j, k)
    ux3 = ux(i-1, j, k)
    uy2 = uy(i, j+1, k)
    uy4 = uy(i, j-1, k)
    uz5 = uz(i, j, k+1)
    uz6 = uz(i, j, k-1)

    if( psi0 >= 0 ) then
      ux0 = 0.0
      uy0 = 0.0
      uz0 = 0.0

      ux1 = 0.0
      ux3 = 0.0
      uy2 = 0.0
      uy4 = 0.0
      uz5 = 0.0
      uz6 = 0.0
    else
      if( psi1 >= 0 ) then
        c1 = 0.0
        c3 = 1.0/rho3m/(dpsi1 + 0.5)
        dp1 = 0.0
      endif
      if( psi3 >= 0 ) then
        c3 = 0.0
        c1 = 1.0/rho1m/(dpsi3 + 0.5)
        dp3 = 0.0
      endif
      if( psi2 >= 0 ) then
        c2 = 0.0
        c4 = 1.0/rho4m/(dpsi2 + 0.5)
        dp2 = 0.0
      endif
      if( psi4 >= 0 ) then
        c4 = 0.0
        c2 = 1.0/rho2m/(dpsi4 + 0.5)
        dp4 = 0.0
      endif
      if( psi5 >= 0 ) then
        c5 = 0.0
        c6 = 1.0/rho6m/(dpsi5 + 0.5)
        dp5 = 0.0
      endif
      if( psi6 >= 0 ) then
        c6 = 0.0
        c5 = 1.0/rho5m/(dpsi6 + 0.5)
        dp6 = 0.0
      endif

			if( psi1 >= 0 .and. psi3 >= 0 ) then
				c1 = 0.0
				c3 = 0.0
				dp1 = 0.0
				dp3 = 0.0
			endif
			if( psi2 >= 0 .and. psi4 >= 0 ) then
				c2 = 0.0
				c4 = 0.0
				dp2 = 0.0
				dp4 = 0.0
			endif
			if( psi5 >= 0 .and. psi6 >= 0 ) then
				c5 = 0.0
				c6 = 0.0
				dp5 = 0.0
				dp6 = 0.0
			endif
    endif

    vx1 = 0.5*(ux0 + ux1)
    vx3 = 0.5*(ux0 + ux3)
    vy2 = 0.5*(uy0 + uy2)
    vy4 = 0.5*(uy0 + uy4)
    vz5 = 0.5*(uz0 + uz5) 
    vz6 = 0.5*(uz0 + uz6) 

		if( phi(i, j, k)*phi(i+1, j, k) < 0 ) then
!			vx1 = vx3
		end if
		if( phi(i, j, k)*phi(i-1, j, k) < 0 ) then
!			vx3 = vx1
		end if
		if( phi(i, j, k)*phi(i, j+1, k) < 0 ) then
!			vy2 = vy4
		end if
		if( phi(i, j, k)*phi(i, j-1, k) < 0 ) then
!			vy4 = vy2
		end if
		if( phi(i, j, k)*phi(i, j, k+1) < 0 ) then
!			vz5 = vz6
		end if
		if( phi(i, j, k)*phi(i, j, k-1) < 0 ) then
!			vz6 = vz5
		end if

		Uw = 0.0
		Vw = 0.0
		Ww = 0.0
		if( psi(i, j, k) >= 0 ) then
      vx1 = 0.0
      vx3 = 0.0
      vy2 = 0.0
      vy4 = 0.0
      vz5 = 0.0
      vz6 = 0.0
		else
			if( psi(i+1, j, k) >= 0 ) then
				xi = x + dpsi1
				yi = y
				zi = z
				Uw = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				Vw = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				Ww = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
				vx1 = vx3*(dpsi1 - 0.5)/(dpsi1 + 0.5) + Uw/(dpsi1 + 0.5)
			endif
			if( psi(i-1, j, k) >= 0 ) then
				xi = x - dpsi3
				yi = y
				zi = z
				Uw = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				Vw = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				Ww = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
				vx3 = vx1*(dpsi3 - 0.5)/(dpsi3 + 0.5) + Uw/(dpsi3 + 0.5)
			endif
			if( psi(i, j+1, k) >= 0 ) then
				xi = x
				yi = y + dpsi2
				zi = z
				Uw = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				Vw = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				Ww = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
				vy2 = vy4*(dpsi2 - 0.5)/(dpsi2 + 0.5) + Vw/(dpsi2 + 0.5)
			endif
			if( psi(i, j-1, k) >= 0 ) then
				xi = x
				yi = y - dpsi4
				zi = z
				Uw = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				Vw = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				Ww = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
				vy4 = vy2*(dpsi4 - 0.5)/(dpsi4 + 0.5) + Vw/(dpsi4 + 0.5)
			endif
			if( psi(i, j, k+1) >= 0 ) then
				xi = x
				yi = y
				zi = z + dpsi5
				Uw = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				Vw = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				Ww = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
				vz5 = vz6*(dpsi5 - 0.5)/(dpsi5 + 0.5) + Ww/(dpsi5 + 0.5)
			endif
			if( psi(i, j, k-1) >= 0 ) then
				xi = x
				yi = y
				zi = z - dpsi6
				Uw = sim_getux(xi, yi, zi, Omega, zmin, rmin)
				Vw = sim_getuy(xi, yi, zi, Omega, zmin, rmin)
				Ww = sim_getuz(xi, yi, zi, Omega, zmin, rmin)
				vz6 = vz5*(dpsi6 - 0.5)/(dpsi6 + 0.5) + Ww/(dpsi6 + 0.5)
			endif

			if( psi1 >= 0 .and. psi3 >= 0 ) then
				vx1 = 0.0
				vx3 = 0.0
			endif
			if( psi2 >= 0 .and. psi4 >= 0 ) then
				vy2 = 0.0
				vy4 = 0.0
			endif
			if( psi5 >= 0 .and. psi6 >= 0 ) then
				vz5 = 0.0
				vz6 = 0.0
			endif
		endif

    divv = vx1 - vx3 + vy2 - vy4 + vz5 - vz6

    A(i, j, k, 0) = -(c1 + c3) - (c2 + c4) - (c5 + c6)
    A(i, j, k, 1) = c1
    A(i, j, k, 3) = c3
    A(i, j, k, 2) = c2
    A(i, j, k, 4) = c4
    A(i, j, k, 5) = c5
    A(i, j, k, 6) = c6

    b(i, j, k)   = divv &
                  + c1*dp1*(Hphi1 - Hphi0) &
                  + c3*dp3*(Hphi3 - Hphi0) &
                  + c2*dp2*(Hphi2 - Hphi0) &
                  + c4*dp4*(Hphi4 - Hphi0) &
                  + c5*dp5*(Hphi5 - Hphi0) &
                  + c6*dp6*(Hphi6 - Hphi0) 

    p1_(i, j, k) = p0_(i, j, k)

    v(i, j, k, 1) = vx1
    v(i, j, k, 3) = vx3
    v(i, j, k, 2) = vy2
    v(i, j, k, 4) = vy4
    v(i, j, k, 5) = vz5
    v(i, j, k, 6) = vz6

    if( psi0 >= 0 ) then
      A(i, j, k, 0) = 1.0
      A(i, j, k, 1) = 0.0
      A(i, j, k, 3) = 0.0
      A(i, j, k, 2) = 0.0
      A(i, j, k, 4) = 0.0
      A(i, j, k, 5) = 0.0
      A(i, j, k, 6) = 0.0
      b(i, j, k) = 0.0
      p1_(i, j, k) = 0.0
      v(i, j, k, 1) = 0.0
      v(i, j, k, 3) = 0.0
      v(i, j, k, 2) = 0.0
      v(i, j, k, 4) = 0.0
      v(i, j, k, 5) = 0.0
      v(i, j, k, 6) = 0.0
    endif

  end do
  end do
  end do
!$omp end do
!$omp end parallel
end subroutine sim_calc_ab_p

subroutine sim_set_refbox_p( &
									A &
								, b &
								, x &
								, xmin &
								, ymin &
								, zmin &
								, xmax &
								, ymax &
								, zmax &
								, pr &
								, dx &
								, org &
								, sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:6)	:: A
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: b
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: x
	real										:: xmin, ymin, zmin
	real										:: xmax, ymax, zmax
	real										:: pr
	real										:: dx
	real, dimension(3)			:: org
	real										:: x0, y0, z0
	real										:: x1, y1, z1
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
!$omp parallel private(i, j, k) &
!$omp					 private(x0, y0, z0) &
!$omp					 private(x1, y1, z1) 
!$omp do schedule(static, 1)
  do k=1, kx
  do j=1, jx
  do i=1, ix
		x0 = org(1)*dx + (real(i-1))*dx
		y0 = org(2)*dx + (real(j-1))*dx
		z0 = org(3)*dx + (real(k-1))*dx
		x1 = org(1)*dx + (real(i))*dx
		y1 = org(2)*dx + (real(j))*dx
		z1 = org(3)*dx + (real(k))*dx

		if( (x1 >= xmin .and. x0 <= xmax) .and. &
				(y1 >= ymin .and. y0 <= ymax) .and. &
				(z1 >= zmin .and. z0 <= zmax) ) then
			A(i, j, k, 0) = 1.0d0
			A(i, j, k, 3) = 0.0d0
			A(i, j, k, 1) = 0.0d0
			A(i, j, k, 4) = 0.0d0
			A(i, j, k, 2) = 0.0d0
			A(i, j, k, 6) = 0.0d0
			A(i, j, k, 5) = 0.0d0
			b(i, j, k) = pr
		endif
  end do
  end do
  end do
!$omp end do
!$omp end parallel
end subroutine sim_set_refbox_p

subroutine sim_corr_u( &
                    ux1_, uy1_, uz1_ &
                  , v1_ &
                  , ux0_, uy0_, uz0_ &
                  , v0_ &
                  , p0_ &
                  , phi &
                  , psi &
                  , kappa &
                  , rhol, rhog &
                  , sigma &
                  , gx, gy, gz &
									, Omega &
									, zmin &
									, zmax &
									, rmin &
									, rmax &
									, org &
                  , sz, g)
  implicit none
  integer                  :: i, j, k
  integer                  :: ix, jx, kx
  integer                  :: g
  integer, dimension(3)    :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: ux1_, uy1_, uz1_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:6)  :: v1_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: ux0_, uy0_, uz0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:6)  :: v0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: p0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: phi
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: psi
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: kappa
  real                    :: rhol, rhog
  real                    :: sigma
  real                    :: gx, gy, gz
	real										:: Omega, zmin, zmax, rmin, rmax
	real, dimension(3)			:: org
  real                    :: phi0, phi1, phi3, phi2, phi4, phi5, phi6
  real                    :: Hphi0, Hphi1, Hphi3, Hphi2, Hphi4, Hphi5, Hphi6
  real                    :: dphi1, dphi3, dphi2, dphi4, dphi5, dphi6
  real                    :: psi0, psi1, psi3, psi2, psi4, psi5, psi6
  real                    :: Hpsi0, Hpsi1, Hpsi3, Hpsi2, Hpsi4, Hpsi5, Hpsi6
  real                    :: dpsi1, dpsi3, dpsi2, dpsi4, dpsi5, dpsi6
  real                    :: mask1, mask3, mask2, mask4, mask5, mask6
  real                    :: rho0, rho1, rho3, rho2, rho4, rho5, rho6
  real                    :: rho1m, rho3m, rho2m, rho4m, rho5m, rho6m
  real                    :: c1, c3, c2, c4, c5, c6
  real                    :: kappa0, kappa1, kappa3, kappa2, kappa4, kappa5, kappa6
  real                    :: kappa1m, kappa3m, kappa2m, kappa4m, kappa5m, kappa6m
  real                    :: dp1, dp3, dp2, dp4, dp5, dp6
  real                    :: dpx1, dpx3, dpy2, dpy4, dpz5, dpz6
  real                    :: cdpx1, cdpx3, cdpy2, cdpy4, cdpz5, cdpz6
  real                    :: cdpx, cdpy, cdpz
	real										:: Uw, Vw, Ww
	real										:: x, y, z, r2, rl, xi, yi, zi, theta
  real                    :: sim_geth
  real                    :: sim_getd
	real										:: sim_getux
	real										:: sim_getuy
	real										:: sim_getuz
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
!$omp parallel private(i, j, k) &
!$omp          ,private(phi0, phi1, phi3, phi2, phi4, phi5, phi6) &
!$omp          ,private(Hphi0, Hphi1, Hphi3, Hphi2, Hphi4, Hphi5, Hphi6) &
!$omp          ,private(dphi1, dphi3, dphi2, dphi4, dphi5, dphi6) &
!$omp          ,private(psi0, psi1, psi3, psi2, psi4, psi5, psi6) &
!$omp          ,private(Hpsi0, Hpsi1, Hpsi3, Hpsi2, Hpsi4, Hpsi5, Hpsi6) &
!$omp          ,private(dpsi1, dpsi3, dpsi2, dpsi4, dpsi5, dpsi6) &
!$omp          ,private(mask1, mask3, mask2, mask4, mask5, mask6) &
!$omp          ,private(rho0, rho1, rho3, rho2, rho4, rho5, rho6) &
!$omp          ,private(rho1m, rho3m, rho2m, rho4m, rho5m, rho6m) &
!$omp          ,private(c1, c3, c2, c4, c5, c6) &
!$omp          ,private(kappa0, kappa1, kappa3, kappa2, kappa4, kappa5, kappa6) &
!$omp          ,private(kappa1m, kappa3m, kappa2m, kappa4m, kappa5m, kappa6m) &
!$omp          ,private(dp1, dp3, dp2, dp4, dp5, dp6) &
!$omp          ,private(dpx1, dpx3, dpy2, dpy4, dpz5, dpz6) &
!$omp          ,private(cdpx1, cdpx3, cdpy2, cdpy4, cdpz5, cdpz6) &
!$omp          ,private(cdpx, cdpy, cdpz) &
!$omp          ,private(Uw, Vw, Ww) &
!$omp					 ,private(x, y, z, r2, rl, xi, yi, zi, theta)
!$omp do schedule(dynamic, 1)
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
		x = org(1) + float(i) - 0.5
		y = org(2) + float(j) - 0.5
		z = org(3) + float(k) - 0.5

    phi0 = phi(i, j, k)
    phi1 = phi(i+1, j, k)
    phi3 = phi(i-1, j, k)
    phi2 = phi(i, j+1, k)
    phi4 = phi(i, j-1, k)
    phi5 = phi(i, j, k+1)
    phi6 = phi(i, j, k-1)

    Hphi0 = sim_geth(phi0)
    Hphi1 = sim_geth(phi1)
    Hphi3 = sim_geth(phi3)
    Hphi2 = sim_geth(phi2)
    Hphi4 = sim_geth(phi4)
    Hphi5 = sim_geth(phi5)
    Hphi6 = sim_geth(phi6)

    dphi1 = sim_getd(phi0, phi1)
    dphi3 = sim_getd(phi0, phi3)
    dphi2 = sim_getd(phi0, phi2)
    dphi4 = sim_getd(phi0, phi4)
    dphi5 = sim_getd(phi0, phi5)
    dphi6 = sim_getd(phi0, phi6)

    rho0 = rhog + (rhol - rhog)*Hphi0
    rho1 = rhog + (rhol - rhog)*Hphi1
    rho3 = rhog + (rhol - rhog)*Hphi3
    rho2 = rhog + (rhol - rhog)*Hphi2
    rho4 = rhog + (rhol - rhog)*Hphi4
    rho5 = rhog + (rhol - rhog)*Hphi5
    rho6 = rhog + (rhol - rhog)*Hphi6

    rho1m = rho0*dphi1 + rho1*(1.0 - dphi1)
    rho3m = rho0*dphi3 + rho3*(1.0 - dphi3)
    rho2m = rho0*dphi2 + rho2*(1.0 - dphi2)
    rho4m = rho0*dphi4 + rho4*(1.0 - dphi4)
    rho5m = rho0*dphi5 + rho5*(1.0 - dphi5)
    rho6m = rho0*dphi6 + rho6*(1.0 - dphi6)

    c1 = 1.0/rho1m
    c3 = 1.0/rho3m
    c2 = 1.0/rho2m
    c4 = 1.0/rho4m
    c5 = 1.0/rho5m
    c6 = 1.0/rho6m

    kappa0 = kappa(i, j, k)
    kappa1 = kappa(i+1, j, k)
    kappa3 = kappa(i-1, j, k)
    kappa2 = kappa(i, j+1, k)
    kappa4 = kappa(i, j-1, k)
    kappa5 = kappa(i, j, k+1)
    kappa6 = kappa(i, j, k-1)

    kappa1m = kappa0*(1.0 - dphi1) + kappa1*dphi1
    kappa3m = kappa0*(1.0 - dphi3) + kappa3*dphi3
    kappa2m = kappa0*(1.0 - dphi2) + kappa2*dphi2
    kappa4m = kappa0*(1.0 - dphi4) + kappa4*dphi4
    kappa5m = kappa0*(1.0 - dphi5) + kappa5*dphi5
    kappa6m = kappa0*(1.0 - dphi6) + kappa6*dphi6

    dp1 = -sigma*kappa1m
    dp3 = -sigma*kappa3m
    dp2 = -sigma*kappa2m
    dp4 = -sigma*kappa4m
    dp5 = -sigma*kappa5m
    dp6 = -sigma*kappa6m

    psi0 = psi(i, j, k)
    psi1 = psi(i+1, j, k)
    psi3 = psi(i-1, j, k)
    psi2 = psi(i, j+1, k)
    psi4 = psi(i, j-1, k)
    psi5 = psi(i, j, k+1)
    psi6 = psi(i, j, k-1)

    dpsi1 = sim_getd(psi0, psi1)
    dpsi3 = sim_getd(psi0, psi3)
    dpsi2 = sim_getd(psi0, psi2)
    dpsi4 = sim_getd(psi0, psi4)
    dpsi5 = sim_getd(psi0, psi5)
    dpsi6 = sim_getd(psi0, psi6)

    dpx1 = p0_(i+1, j, k) - p0_(i, j, k)
    dpx3 = p0_(i, j, k) - p0_(i-1, j, k)
    dpy2 = p0_(i, j+1, k) - p0_(i, j, k)
    dpy4 = p0_(i, j, k) - p0_(i, j-1, k)
    dpz5 = p0_(i, j, k+1) - p0_(i, j, k)
    dpz6 = p0_(i, j, k) - p0_(i, j, k-1)

		mask1 = 0.0
		mask3 = 0.0
		mask2 = 0.0
		mask4 = 0.0
		mask5 = 0.0
		mask6 = 0.0
		if( psi(i, j, k) >= 0 ) then
		else
			if( psi(i+1, j, k) >= 0 ) then
				mask1 = 1.0
			endif
			if( psi(i-1, j, k) >= 0 ) then
				mask3 = 1.0
			endif
			if( psi(i, j+1, k) >= 0 ) then
				mask2 = 1.0
			endif
			if( psi(i, j-1, k) >= 0 ) then
				mask4 = 1.0
			endif
			if( psi(i, j, k+1) >= 0 ) then
				mask5 = 1.0
			endif
			if( psi(i, j, k-1) >= 0 ) then
				mask6 = 1.0
			endif
		endif

    if( psi0 >= 0 ) then
      dpx1 = 0.0
      dpx3 = 0.0
      dpy2 = 0.0
      dpy4 = 0.0
      dpz5 = 0.0
      dpz6 = 0.0
    else
      if( psi1 >= 0 ) then
        c1 = c3
        dpx1 = dpx3*(dpsi1 - 0.5)/(dpsi1 + 0.5) 
        dp1 = 0.0
      endif
      if( psi3 >= 0 ) then
        c3 = c1
        dpx3 = dpx1*(dpsi3 - 0.5)/(dpsi3 + 0.5) 
        dp3 = 0.0
      endif
      if( psi2 >= 0 ) then
        c2 = c4
        dpy2 = dpy4*(dpsi2 - 0.5)/(dpsi2 + 0.5)
        dp2 = 0.0
      endif
      if( psi4 >= 0 ) then
        c4 = c2
        dpy4 = dpy2*(dpsi4 - 0.5)/(dpsi4 + 0.5) 
        dp4 = 0.0
      endif
      if( psi5 >= 0 ) then
        c5 = c6
        dpz5 = dpz6*(dpsi5 - 0.5)/(dpsi5 + 0.5)
        dp5 = 0.0
      endif
      if( psi6 >= 0 ) then
        c6 = c5
        dpz6 = dpz5*(dpsi6 - 0.5)/(dpsi6 + 0.5)
        dp6 = 0.0
      endif

			if( psi1 >= 0 .and. psi3 >= 0 ) then
				c1 = 0.0
				c3 = 0.0
				dpx1 = 0.0
				dpx3 = 0.0
				dp1 = 0.0
				dp3 = 0.0
			endif
			if( psi2 >= 0 .and. psi4 >= 0 ) then
				c2 = 0.0
				c4 = 0.0
				dpy2 = 0.0
				dpy4 = 0.0
				dp2 = 0.0
				dp4 = 0.0
			endif
			if( psi5 >= 0 .and. psi6 >= 0 ) then
				c5 = 0.0
				c6 = 0.0
				dpz5 = 0.0
				dpz6 = 0.0
				dp5 = 0.0
				dp6 = 0.0
			endif
    endif

    cdpx1 = c1*dpx1 - c1*dp1*(Hphi1 - Hphi0)
    cdpx3 = c3*dpx3 - c3*dp3*(Hphi0 - Hphi3)
    cdpy2 = c2*dpy2 - c2*dp2*(Hphi2 - Hphi0)
    cdpy4 = c4*dpy4 - c4*dp4*(Hphi0 - Hphi4)
    cdpz5 = c5*dpz5 - c5*dp5*(Hphi5 - Hphi0)
    cdpz6 = c6*dpz6 - c6*dp6*(Hphi0 - Hphi6)

    cdpx = 0.5*(cdpx1 + cdpx3)
    cdpy = 0.5*(cdpy2 + cdpy4)
    cdpz = 0.5*(cdpz5 + cdpz6)

    ux1_(i, j, k)   = ux0_(i, j, k) - cdpx - 0.5*gx/(dpsi1 + 0.5)*mask1 - 0.5*gx/(dpsi3 + 0.5)*mask3
    uy1_(i, j, k)   = uy0_(i, j, k) - cdpy - 0.5*gy/(dpsi2 + 0.5)*mask2 - 0.5*gy/(dpsi4 + 0.5)*mask4
    uz1_(i, j, k)   = uz0_(i, j, k) - cdpz - 0.5*gz/(dpsi5 + 0.5)*mask5 - 0.5*gz/(dpsi6 + 0.5)*mask6
    v1_(i, j, k, 1) = v0_(i, j, k, 1) - cdpx1
    v1_(i, j, k, 3) = v0_(i, j, k, 3) - cdpx3
    v1_(i, j, k, 2) = v0_(i, j, k, 2) - cdpy2
    v1_(i, j, k, 4) = v0_(i, j, k, 4) - cdpy4
    v1_(i, j, k, 5) = v0_(i, j, k, 5) - cdpz5
    v1_(i, j, k, 6) = v0_(i, j, k, 6) - cdpz6

!      if( psi6 >= 0 ) then
!        write(*, *) v0_(i, j, k, 6), v1_(i, j, k, 6), uz0_(i, j, k), uz1_(i, j, k)
!        write(*,*) dpsi6
!      end if

    if( psi0 >= 0 ) then
      ux1_(i, j, k) = 0.0
      uy1_(i, j, k) = 0.0
      uz1_(i, j, k) = 0.0
      v1_(i, j, k, 1) = 0.0
      v1_(i, j, k, 3) = 0.0
      v1_(i, j, k, 2) = 0.0
      v1_(i, j, k, 4) = 0.0
      v1_(i, j, k, 5) = 0.0
      v1_(i, j, k, 6) = 0.0
    endif

  end do
  end do
  end do
!$omp end do
!$omp end parallel
end subroutine sim_corr_u

subroutine sim_set_ic_u( &
									ux_, uy_, uz_ &
								, phi &
								, psi &
								, Omega &
								, zmin &
								, zmax &
								, rmin &
								, rmax &
								, org &
								, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: ux_, uy_, uz_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: psi
	real										:: Omega, zmin, zmax, rmin, rmax
	real, dimension(3)			:: org
	real										:: phi0, phi1, phi3, phi2, phi4, phi5, phi6
	real										:: Hphi0, Hphi1, Hphi3, Hphi2, Hphi4, Hphi5, Hphi6
	real										:: dphi1, dphi3, dphi2, dphi4, dphi5, dphi6
	real										:: psi0, psi1, psi3, psi2, psi4, psi5, psi6
	real										:: Hpsi0, Hpsi1, Hpsi3, Hpsi2, Hpsi4, Hpsi5, Hpsi6
	real										:: dpsi1, dpsi3, dpsi2, dpsi4, dpsi5, dpsi6
	real										:: mask1, mask3, mask2, mask4, mask5, mask6
	real										:: muf0, muf1, muf3, muf2, muf4, muf5, muf6
	real										:: mu0, mu1, mu3, mu2, mu4, mu5, mu6
	real										:: mu1m, mu3m, mu2m, mu4m, mu5m, mu6m
	real										:: rhof0, rho0
	real										:: a1, a3, a2, a4, a5, a6
	real										:: u0, u1, u3, u2, u4, u5, u6
	real										:: ux0, ux1, ux3, ux2, ux4, ux5, ux6
	real										:: uy0, uy1, uy3, uy2, uy4, uy5, uy6
	real										:: uz0, uz1, uz3, uz2, uz4, uz5, uz6
	real										:: duxdx1, duxdx3, duydx2, duydx4, duzdx5, duzdx6
	real										:: duxdy1, duxdy3, duydy2, duydy4, duzdy5, duzdy6
	real										:: duxdz1, duxdz3, duydz2, duydz4, duzdz5, duzdz6
	real										:: Uw
	real										:: x, y, z, r2, rl, xi, yi, zi, theta
	real										:: alpha
	real										:: sim_geth
	real										:: sim_geth_s
	real										:: sim_geth_2
	real										:: sim_getd
	real										:: sim_getux
	real										:: sim_getuy
	real										:: sim_getuz
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	alpha = 0.5
	alpha = 1.0
!$omp parallel private(i, j, k) &
!$omp					,private(phi0, phi1, phi3, phi2, phi4, phi5, phi6) &
!$omp					,private(Hphi0, Hphi1, Hphi3, Hphi2, Hphi4, Hphi5, Hphi6) &
!$omp					,private(dphi1, dphi3, dphi2, dphi4, dphi5, dphi6) &
!$omp					,private(psi0, psi1, psi3, psi2, psi4, psi5, psi6) &
!$omp					,private(Hpsi0, Hpsi1, Hpsi3, Hpsi2, Hpsi4, Hpsi5, Hpsi6) &
!$omp					,private(dpsi1, dpsi3, dpsi2, dpsi4, dpsi5, dpsi6) &
!$omp					,private(mask1, mask3, mask2, mask4, mask5, mask6) &
!$omp					,private(muf0, muf1, muf3, muf2, muf4, muf5, muf6) &
!$omp					,private(mu0, mu1, mu3, mu2, mu4, mu5, mu6) &
!$omp					,private(mu1m, mu3m, mu2m, mu4m, mu5m, mu6m) &
!$omp					,private(rhof0, rho0) &
!$omp					,private(a1, a3, a2, a4, a5, a6) &
!$omp					,private(u0, u1, u3, u2, u4, u5, u6) &
!$omp					,private(ux0, ux1, ux3, ux2, ux4, ux5, ux6) &
!$omp					,private(uy0, uy1, uy3, uy2, uy4, uy5, uy6) &
!$omp					,private(uz0, uz1, uz3, uz2, uz4, uz5, uz6) &
!$omp					,private(duxdx1, duxdx3, duydx2, duydx4, duzdx5, duzdx6) &
!$omp					,private(duxdy1, duxdy3, duydy2, duydy4, duzdy5, duzdy6) &
!$omp					,private(duxdz1, duxdz3, duydz2, duydz4, duzdz5, duzdz6) &
!$omp					,private(Uw) &
!$omp					,private(x, y, z, r2, rl, xi, yi, zi, theta)
!$omp do schedule(dynamic, 1)
	do k=1, kx
	do j=1, jx
!ocl nouxsimd
	do i=1, ix
		x = org(1) + float(i) - 0.5
		y = org(2) + float(j) - 0.5
		z = org(3) + float(k) - 0.5

		ux_(i, j, k) = sim_getux(x, y, z, Omega, 0.0, 0.0)
		uy_(i, j, k) = sim_getuy(x, y, z, Omega, 0.0, 0.0)
		uz_(i, j, k) = sim_getuz(x, y, z, Omega, 0.0, 0.0)

		if( psi(i, j, k) >= 0 ) then
			ux_(i, j, k) = 0.0
			uy_(i, j, k) = 0.0
			uz_(i, j, k) = 0.0
		endif

	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine sim_set_ic_u

