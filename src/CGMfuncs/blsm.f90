!>  @file  blsm.f90
!!  @brief Basic subprograms for the level-set method for Cartesian grid data structure
!<

function lsm_getd(phi0, phi1)
  implicit none
  real          :: lsm_getd
  real          :: phi0, phi1
  lsm_getd = abs(phi0)/(abs(phi0) + abs(phi1))
  if( phi0 >= 0 .and. phi1 >= 0 ) then
    lsm_getd = 1.0
  endif
  if( phi0 < 0 .and. phi1 < 0 ) then
    lsm_getd = 1.0
  endif
  return
end function lsm_getd

function lsm_getsign(a)
	implicit none
	real					:: lsm_getsign
	real					:: a
	lsm_getsign = 1.0
	if( a < 0 ) then
		lsm_getsign = -1.0
	end if
	return
end function lsm_getsign

function lsm_getminmod(a, b)
	implicit none
	real					:: lsm_getminmod
	real					:: a, b
!	lsm_getminmod = 0.5*(sign(a, a) + sign(b, b))*min(abs(a), abs(b))
	lsm_getminmod = 0.0
	if( abs(a) .lt. abs(b) ) then
		lsm_getminmod = a
	else if( abs(a) .gt. abs(b) ) then
		lsm_getminmod = b
	endif
	return
end function lsm_getminmod

function lsm_getupwind(u, n, p)
	implicit none
	real					:: lsm_getupwind
	real					:: u
	real					:: n, p
	lsm_getupwind = n
	if( u .lt. 0 ) then
		lsm_getupwind = p
	endif
	return
end function lsm_getupwind

function lsm_getupwind1(u, n, p)
	implicit none
	real					:: lsm_getupwind1
	real					:: u
	real					:: n, p
	lsm_getupwind1 = 0.0
	if( u .gt. 0 ) then
		lsm_getupwind1 = n
	else if( u .lt. 0 ) then
		lsm_getupwind1 = p
	endif
	return
end function lsm_getupwind1

function lsm_getupwind2(u, n, p)
	implicit none
	real					:: lsm_getupwind2
	real					:: u
	real					:: n, p
	lsm_getupwind2 = 0.5*(p + n)
	if( u .gt. 0 ) then
		lsm_getupwind2 = n
	else if( u .lt. 0 ) then
		lsm_getupwind2 = p
	endif
	return
end function lsm_getupwind2

function lsm_getweno3(v1, v2, v3)
	implicit none
	real					:: lsm_getweno3
	real					:: v1, v2, v3
	real					:: r
	real					:: w
	real					:: eps
	eps = 0.000001
	r = (eps + (v2-v1)**2)/(eps + (v3-v2)**2)
	w = 1.0/(1.0 + 2*r*r)
	lsm_getweno3 = 0.5*(v2 + v3) - 0.5*w*(v1 - 2.0*v2 + v3)
	return
end function lsm_getweno3

function lsm_getintval(f, xi, yi, zi, sz, g)
	implicit none
	real					:: lsm_getintval
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: f
	real					:: xi, yi, zi
	integer				:: i0, j0, k0
	real					:: x0, y0, z0
	real					:: dx, dy, dz
	real					:: f_ccc, f_pcc, f_cpc, f_ccp, f_ppc, f_cpp, f_pcp, f_ppp
	i0 = floor(xi + 0.5) 
	j0 = floor(yi + 0.5) 
	k0 = floor(zi + 0.5) 
	x0 = real(i0) - 0.5
	y0 = real(j0) - 0.5
	z0 = real(k0) - 0.5
	dx = xi - x0
	dy = yi - y0
	dz = zi - z0
	f_ccc = f(i0, j0, k0)
	f_pcc = f(i0+1, j0, k0)
	f_cpc = f(i0, j0+1, k0)
	f_ccp = f(i0, j0, k0+1)
	f_ppc = f(i0+1, j0+1, k0)
	f_cpp = f(i0, j0+1, k0+1)
	f_pcp = f(i0+1, j0, k0+1)
	f_ppp = f(i0+1, j0+1, k0+1)
	lsm_getintval = &
						+ f_ccc*(1.0 - dx)*(1.0 - dy)*(1.0 - dz) &
						+ f_pcc*       dx *(1.0 - dy)*(1.0 - dz) &
						+ f_cpc*(1.0 - dx)*       dy *(1.0 - dz) &
						+ f_ccp*(1.0 - dx)*(1.0 - dy)*       dz  &
						+ f_ppc*       dx*        dy *(1.0 - dz) &
						+ f_cpp*(1.0 - dx)*       dy *       dz  &
						+ f_pcp*       dx *(1.0 - dy)*       dz  &
						+ f_ppp*       dx *       dy *       dz  
	return
end function lsm_getintval

subroutine lsm_calc_n(nx, ny, nz, divn, phi, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: nx, ny, nz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: divn
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi
	real										:: dpx_c, dpy_c, dpz_c
	real										:: dpx_p, dpy_p, dpz_p
	real										:: dpx_n, dpy_n, dpz_n
	real										:: dp2_ppp, dp2_ppn, dp2_pnp, dp2_npp
	real										:: dp2_nnn, dp2_nnp, dp2_npn, dp2_pnn
	real										:: nx_tmp, ny_tmp, nz_tmp, n2_tmp, nl_tmp
	real										:: dpl_c, dp2_c
	real										:: dpxx_c, dpyy_c, dpzz_c
	real										:: dpxy_c, dpyz_c, dpzx_c
	real										:: eps
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	eps = 1.0e-16
!$omp parallel private(i, j, k) &
!$omp					,private(dpx_c, dpy_c, dpz_c) &
!$omp					,private(dpx_p, dpy_p, dpz_p) &
!$omp					,private(dpx_n, dpy_n, dpz_n) &
!$omp					,private(dp2_ppp, dp2_ppn, dp2_pnp, dp2_npp) &
!$omp					,private(dp2_nnn, dp2_nnp, dp2_npn, dp2_pnn) &
!$omp					,private(nx_tmp, ny_tmp, nz_tmp, n2_tmp, nl_tmp) &
!$omp					,private(dpl_c, dp2_c) &
!$omp					,private(dpxx_c, dpyy_c, dpzz_c) &
!$omp					,private(dpxy_c, dpyz_c, dpzx_c)
!$omp do schedule(dynamic, 1)
	do k=1, kx
	do j=1, jx
!ocl nouxsimd
	do i=1, ix
		dpx_p = phi(i+1, j, k) - phi(i, j, k)
		dpx_n = phi(i, j, k) - phi(i-1, j, k)
		dpy_p = phi(i, j+1, k) - phi(i, j, k)
		dpy_n = phi(i, j, k) - phi(i, j-1, k)
		dpz_p = phi(i, j, k+1) - phi(i, j, k)
		dpz_n = phi(i, j, k) - phi(i, j, k-1)
		dp2_ppp = dpx_p*dpx_p + dpy_p*dpy_p + dpz_p*dpz_p
		dp2_ppn = dpx_p*dpx_p + dpy_p*dpy_p + dpz_n*dpz_n
		dp2_pnp = dpx_p*dpx_p + dpy_n*dpy_n + dpz_p*dpz_p
		dp2_npp = dpx_n*dpx_n + dpy_p*dpy_p + dpz_p*dpz_p
		dp2_nnp = dpx_n*dpx_n + dpy_n*dpy_n + dpz_p*dpz_p
		dp2_npn = dpx_n*dpx_n + dpy_p*dpy_p + dpz_n*dpz_n
		dp2_pnn = dpx_p*dpx_p + dpy_n*dpy_n + dpz_n*dpz_n
		dp2_nnn = dpx_n*dpx_n + dpy_n*dpy_n + dpz_n*dpz_n
		nx_tmp = dpx_p/sqrt(dp2_ppp) &
					 + dpx_p/sqrt(dp2_ppn) &
					 + dpx_p/sqrt(dp2_pnp) &
					 + dpx_n/sqrt(dp2_npp) &
					 + dpx_n/sqrt(dp2_nnp) &
					 + dpx_n/sqrt(dp2_npn) &
					 + dpx_p/sqrt(dp2_pnn) &
					 + dpx_n/sqrt(dp2_nnn) 
		ny_tmp = dpy_p/sqrt(dp2_ppp) &
					 + dpy_p/sqrt(dp2_ppn) &
					 + dpy_n/sqrt(dp2_pnp) &
					 + dpy_p/sqrt(dp2_npp) &
					 + dpy_n/sqrt(dp2_nnp) &
					 + dpy_p/sqrt(dp2_npn) &
					 + dpy_n/sqrt(dp2_pnn) &
					 + dpy_n/sqrt(dp2_nnn) 
		nz_tmp = dpz_p/sqrt(dp2_ppp) &
					 + dpz_n/sqrt(dp2_ppn) &
					 + dpz_p/sqrt(dp2_pnp) &
					 + dpz_p/sqrt(dp2_npp) &
					 + dpz_p/sqrt(dp2_nnp) &
					 + dpz_n/sqrt(dp2_npn) &
					 + dpz_n/sqrt(dp2_pnn) &
					 + dpz_n/sqrt(dp2_nnn) 
		n2_tmp = nx_tmp*nx_tmp + ny_tmp*ny_tmp + nz_tmp*nz_tmp
		nl_tmp = sqrt(n2_tmp)
		dpx_c = 0.5*(phi(i+1, j, k) - phi(i-1, j, k))
		dpy_c = 0.5*(phi(i, j+1, k) - phi(i, j-1, k))
		dpz_c = 0.5*(phi(i, j, k+1) - phi(i, j, k-1))
		dp2_c = dpx_c*dpx_c + dpy_c*dpy_c + dpz_c*dpz_c
		dpl_c = sqrt(dp2_c);
		dpxx_c = phi(i+1, j, k) + phi(i-1, j, k) - 2.0*phi(i, j, k)
		dpyy_c = phi(i, j+1, k) + phi(i, j-1, k) - 2.0*phi(i, j, k)
		dpzz_c = phi(i, j, k+1) + phi(i, j, k-1) - 2.0*phi(i, j, k)
		dpxy_c = (phi(i+1, j+1, k) + phi(i-1, j-1, k) - phi(i+1, j-1, k) - phi(i-1, j+1, k))*0.25
		dpyz_c = (phi(i, j+1, k+1) + phi(i, j-1, k-1) - phi(i, j-1, k+1) - phi(i, j+1, k-1))*0.25
		dpzx_c = (phi(i+1, j, k+1) + phi(i-1, j, k-1) - phi(i+1, j, k-1) - phi(i-1, j, k+1))*0.25
		nx(i, j, k) = dpx_c/dpl_c
		ny(i, j, k) = dpy_c/dpl_c
		nz(i, j, k) = dpz_c/dpl_c
		nx(i, j, k) = nx_tmp/nl_tmp
		ny(i, j, k) = ny_tmp/nl_tmp
		nz(i, j, k) = nz_tmp/nl_tmp
		divn(i, j, k) = ( &
										+ (dpyy_c + dpzz_c)*dpx_c*dpx_c &
										+ (dpzz_c + dpxx_c)*dpy_c*dpy_c &
										+ (dpxx_c + dpyy_c)*dpz_c*dpz_c &
										- 2.0*dpx_c*dpy_c*dpxy_c &
										- 2.0*dpy_c*dpz_c*dpyz_c &
										- 2.0*dpz_c*dpx_c*dpzx_c &
									)/(dpl_c**3 + eps);
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine lsm_calc_n

subroutine lsm_calc_vin(vin, phi, psi, v, nx, ny, nz, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: vin
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: psi
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:6)	:: v
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: nx, ny, nz
	real										:: vx, vy, vz
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(vx, vy, vz)
!$omp do 
	do k=1, kx
	do j=1, jx
	do i=1, ix
		vx = 0.5*(v(i, j, k, 1) + v(i, j, k, 3))
		vy = 0.5*(v(i, j, k, 2) + v(i, j, k, 4))
		vz = 0.5*(v(i, j, k, 5) + v(i, j, k, 6))
		vin(i, j, k) = vx*nx(i, j, k) + vy*ny(i, j, k) + vz*nz(i, j, k)
		if( psi(i, j, k) >= 0.0 ) then
			vin(i, j, k) = 0.0d0
		end if
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine lsm_calc_vin

subroutine lsm_calc_i(fi, f0, phi, phinx, phiny, phinz, psi, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: fi
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: f0
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phinx, phiny, phinz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: psi
	real										:: phi0, phi1, phi3, phi2, phi4, phi5, phi6
	real										:: x, y, z
	real										:: xi, yi, zi
	real										:: nx, ny, nz, n2, nl
	integer									:: i0, j0, k0
	real										:: x0, y0, z0
	real										:: dx, dy, dz
	real										:: f_ccc, f_pcc, f_cpc, f_ccp, f_ppc, f_cpp, f_pcp, f_ppp
	real										:: lsm_getintval
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(phi0, phi1, phi3, phi2, phi4, phi5, phi6) &
!$omp					,private(x, y, z) &
!$omp					,private(xi, yi, zi) &
!$omp					,private(nx, ny, nz, n2, nl) &
!$omp					,private(i0, j0, k0) &
!$omp					,private(x0, y0, z0) &
!$omp					,private(dx, dy, dz) &
!$omp					,private(f_ccc, f_pcc, f_cpc, f_ccp, f_ppc, f_cpp, f_pcp, f_ppp)
!$omp do schedule(dynamic, 1)
	do k=1, kx
	do j=1, jx
	do i=1, ix
		fi(i, j, k) = 0.0
		if( psi(i, j, k) >= 0 ) then
			cycle
		endif
		phi0 = phi(i, j, k)		
		phi1 = phi(i+1, j, k)		
		phi3 = phi(i-1, j, k)		
		phi2 = phi(i, j+1, k)		
		phi4 = phi(i, j-1, k)		
		phi5 = phi(i, j, k+1)		
		phi6 = phi(i, j, k-1)		
		if( phi0*phi1 <= 0 .or. &
				phi0*phi3 <= 0 .or. &
				phi0*phi2 <= 0 .or. &
				phi0*phi4 <= 0 .or. &
				phi0*phi5 <= 0 .or. &
				phi0*phi6 <= 0 ) then
			if( abs(phi0) >= 1.0 ) then
				fi(i, j, k) = f0(i, j, k)
			else
				x = real(i) - 0.5
				y = real(j) - 0.5
				z = real(k) - 0.5
				xi = x - phi0*phinx(i, j, k)
				yi = y - phi0*phiny(i, j, k)
				zi = z - phi0*phinz(i, j, k)
				fi(i, j, k) = lsm_getintval(f0, xi, yi, zi, sz, g)
			endif
		endif
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine lsm_calc_i

subroutine lsm_update_phi_v(phi1_, phi0_, v0_, psi_, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi1_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi0_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:6)	:: v0_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: psi_
	real										:: dpx_p2, dpy_p2, dpz_p2
	real										:: dpx_p1, dpy_p1, dpz_p1
	real										:: dpx_n1, dpy_n1, dpz_n1
	real										:: dpx_n2, dpy_n2, dpz_n2
	real										:: dpx_p, dpy_p, dpz_p
	real										:: dpx_c, dpy_c, dpz_c
	real										:: dpx_n, dpy_n, dpz_n
	real										:: phi1p, phi3p, phi2p, phi4p, phi5p, phi6p
	real										:: phi1n, phi3n, phi2n, phi4n, phi5n, phi6n
	real										:: phi0, phi1, phi3, phi2, phi4, phi5, phi6
	real										:: psi0, psi1, psi3, psi2, psi4, psi5, psi6
	real										:: dpsi0, dpsi1, dpsi3, dpsi2, dpsi4, dpsi5, dpsi6
	real										:: q1, q3, q2, q4, q5, q6
	real										:: lsm_getminmod
	real										:: lsm_getupwind
	real										:: lsm_getd
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(dpx_p2, dpy_p2, dpz_p2) &
!$omp					,private(dpx_p1, dpy_p1, dpz_p1) &
!$omp					,private(dpx_n1, dpy_n1, dpz_n1) &
!$omp					,private(dpx_n2, dpy_n2, dpz_n2) &
!$omp					,private(dpx_p, dpy_p, dpz_p) &
!$omp					,private(dpx_c, dpy_c, dpz_c) &
!$omp					,private(dpx_n, dpy_n, dpz_n) &
!$omp					,private(phi1p, phi3p, phi2p, phi4p, phi5p, phi6p) &
!$omp					,private(phi1n, phi3n, phi2n, phi4n, phi5n, phi6n) &
!$omp					,private(phi0, phi1, phi3, phi2, phi4, phi5, phi6) &
!$omp					,private(psi0, psi1, psi3, psi2, psi4, psi5, psi6) &
!$omp					,private(dpsi0, dpsi1, dpsi3, dpsi2, dpsi4, dpsi5, dpsi6) &
!$omp					,private(q1, q3, q2, q4, q5, q6)
!$omp do schedule(dynamic, 1)
	do k=1, kx
	do j=1, jx
!ocl nouxsimd
	do i=1, ix
		dpx_p2 = phi0_(i+2, j, k) - phi0_(i+1, j, k)
		dpx_p1 = phi0_(i+1, j, k) - phi0_(i  , j, k)
		dpx_n1 = phi0_(i  , j, k) - phi0_(i-1, j, k)
		dpx_n2 = phi0_(i-1, j, k) - phi0_(i-2, j, k)

		dpy_p2 = phi0_(i, j+2, k) - phi0_(i, j+1, k)
		dpy_p1 = phi0_(i, j+1, k) - phi0_(i, j  , k)
		dpy_n1 = phi0_(i, j  , k) - phi0_(i, j-1, k)
		dpy_n2 = phi0_(i, j-1, k) - phi0_(i, j-2, k)

		dpz_p2 = phi0_(i, j, k+2) - phi0_(i, j, k+1)
		dpz_p1 = phi0_(i, j, k+1) - phi0_(i, j, k  )
		dpz_n1 = phi0_(i, j, k  ) - phi0_(i, j, k-1)
		dpz_n2 = phi0_(i, j, k-1) - phi0_(i, j, k-2)

		dpx_p = lsm_getminmod(dpx_p2, dpx_p1)
		dpx_c = lsm_getminmod(dpx_p1, dpx_n1)
		dpx_n = lsm_getminmod(dpx_n1, dpx_n2)
		dpy_p = lsm_getminmod(dpy_p2, dpy_p1)
		dpy_c = lsm_getminmod(dpy_p1, dpy_n1)
		dpy_n = lsm_getminmod(dpy_n1, dpy_n2)
		dpz_p = lsm_getminmod(dpz_p2, dpz_p1)
		dpz_c = lsm_getminmod(dpz_p1, dpz_n1)
		dpz_n = lsm_getminmod(dpz_n1, dpz_n2)

		phi1p = phi0_(i+1, j, k) - 0.5*dpx_p
		phi1n = phi0_(i  , j, k) + 0.5*dpx_c
		phi3p = phi0_(i  , j, k) - 0.5*dpx_c
		phi3n = phi0_(i-1, j, k) + 0.5*dpx_n

		phi2p = phi0_(i, j+1, k) - 0.5*dpy_p
		phi2n = phi0_(i, j  , k) + 0.5*dpy_c
		phi4p = phi0_(i, j  , k) - 0.5*dpy_c
		phi4n = phi0_(i, j-1, k) + 0.5*dpy_n

		phi5p = phi0_(i, j, k+1) - 0.5*dpz_p
		phi5n = phi0_(i, j, k  ) + 0.5*dpz_c
		phi6p = phi0_(i, j, k  ) - 0.5*dpz_c
		phi6n = phi0_(i, j, k-1) + 0.5*dpz_n

		phi1 = lsm_getupwind(v0_(i, j, k, 1), phi1n, phi1p)
		phi3 = lsm_getupwind(v0_(i, j, k, 3), phi3n, phi3p)
		phi2 = lsm_getupwind(v0_(i, j, k, 2), phi2n, phi2p)
		phi4 = lsm_getupwind(v0_(i, j, k, 4), phi4n, phi4p)
		phi5 = lsm_getupwind(v0_(i, j, k, 5), phi5n, phi5p)
		phi6 = lsm_getupwind(v0_(i, j, k, 6), phi6n, phi6p)

		if( psi_(i+1, j, k) >= 0 ) then
			phi1 = 0.0
			phi1 = phi0_(i, j, k) + (dpsi1 - 0.5)/(dpsi1 + 0.5)*(phi0_(i, j, k) - phi0_(i-1, j, k))
			phi1 = phi0_(i, j, k) + 0.5*(dpsi1 - 0.5)/(dpsi1 + 0.5)*(phi0_(i, j, k) - phi0_(i-1, j, k))
			phi1 = phi0_(i, j, k)
		end if
		if( psi_(i-1, j, k) >= 0 ) then
			phi3 = 0.0
			phi3 = phi0_(i, j, k) - (dpsi3 - 0.5)/(dpsi3 + 0.5)*(phi0_(i+1, j, k) - phi0_(i, j, k))
			phi3 = phi0_(i, j, k) - 0.5*(dpsi3 - 0.5)/(dpsi3 + 0.5)*(phi0_(i+1, j, k) - phi0_(i, j, k))
			phi3 = phi0_(i, j, k)
		end if
		if( psi_(i, j+1, k) >= 0 ) then
			phi2 = 0.0
			phi2 = phi0_(i, j, k) + (dpsi2 - 0.5)/(dpsi2 + 0.5)*(phi0_(i, j, k) - phi0_(i, j-1, k))
			phi2 = phi0_(i, j, k) + 0.5*(dpsi2 - 0.5)/(dpsi2 + 0.5)*(phi0_(i, j, k) - phi0_(i, j-1, k))
			phi2 = phi0_(i, j, k)
		end if
		if( psi_(i, j-1, k) >= 0 ) then
			phi4 = 0.0
			phi4 = phi0_(i, j, k) - (dpsi4 - 0.5)/(dpsi4 + 0.5)*(phi0_(i, j+1, k) - phi0_(i, j, k))
			phi4 = phi0_(i, j, k) - 0.5*(dpsi4 - 0.5)/(dpsi4 + 0.5)*(phi0_(i, j+1, k) - phi0_(i, j, k))
			phi4 = phi0_(i, j, k)
		end if
		if( psi_(i, j, k+1) >= 0 ) then
			phi5 = 0.0
			phi5 = phi0_(i, j, k) + (dpsi5 - 0.5)/(dpsi5 + 0.5)*(phi0_(i, j, k) - phi0_(i, j, k-1))
			phi5 = phi0_(i, j, k) + 0.5*(dpsi5 - 0.5)/(dpsi5 + 0.5)*(phi0_(i, j, k) - phi0_(i, j, k-1))
			phi5 = phi0_(i, j, k)
		end if
		if( psi_(i, j, k-1) >= 0 ) then
			phi6 = 0.0
			phi6 = phi0_(i, j, k) - (dpsi6 - 0.5)/(dpsi6 + 0.5)*(phi0_(i, j, k+1) - phi0_(i, j, k))
			phi6 = phi0_(i, j, k) - 0.5*(dpsi6 - 0.5)/(dpsi6 + 0.5)*(phi0_(i, j, k+1) - phi0_(i, j, k))
			phi6 = phi0_(i, j, k)
		end if

		q1 = phi1*v0_(i, j, k, 1)
		q3 = phi3*v0_(i, j, k, 3)
		q2 = phi2*v0_(i, j, k, 2)
		q4 = phi4*v0_(i, j, k, 4)
		q5 = phi5*v0_(i, j, k, 5)
		q6 = phi6*v0_(i, j, k, 6)

		phi1_(i, j, k) = phi0_(i, j, k) &
											- (q1 - q3) &
											- (q2 - q4) &
											- (q5 - q6)
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine lsm_update_phi_v

subroutine lsm_update_phi_v_2(phi1_, phi0_, v0_, psi, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi1_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi0_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:6)	:: v0_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: psi
	real										:: dpx_p2, dpy_p2, dpz_p2
	real										:: dpx_p1, dpy_p1, dpz_p1
	real										:: dpx_n1, dpy_n1, dpz_n1
	real										:: dpx_n2, dpy_n2, dpz_n2
	real										:: dpx_p, dpy_p, dpz_p
	real										:: dpx_c, dpy_c, dpz_c
	real										:: dpx_n, dpy_n, dpz_n
	real										:: phi1p, phi3p, phi2p, phi4p, phi5p, phi6p
	real										:: phi1n, phi3n, phi2n, phi4n, phi5n, phi6n
	real										:: phi0, phi1, phi3, phi2, phi4, phi5, phi6
	real										:: lsm_getminmod
	real										:: lsm_getupwind
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(dpx_p2, dpy_p2, dpz_p2) &
!$omp					,private(dpx_p1, dpy_p1, dpz_p1) &
!$omp					,private(dpx_n1, dpy_n1, dpz_n1) &
!$omp					,private(dpx_n2, dpy_n2, dpz_n2) &
!$omp					,private(dpx_p, dpy_p, dpz_p) &
!$omp					,private(dpx_c, dpy_c, dpz_c) &
!$omp					,private(dpx_n, dpy_n, dpz_n) &
!$omp					,private(phi1p, phi3p, phi2p, phi4p, phi5p, phi6p) &
!$omp					,private(phi1n, phi3n, phi2n, phi4n, phi5n, phi6n) &
!$omp					,private(phi0, phi1, phi3, phi2, phi4, phi5, phi6)
!$omp do schedule(dynamic, 1)
	do k=1, kx
	do j=1, jx
!ocl nouxsimd
	do i=1, ix
		phi1p = 0.125*(3.0*phi0_(i, j, k) + 6.0*phi0_(i+1, j, k) - phi0_(i+2, j, k))
		phi1n = 0.125*(6.0*phi0_(i, j, k) + 3.0*phi0_(i+1, j, k) - phi0_(i-1, j, k))
		phi3p = 0.125*(6.0*phi0_(i, j, k) + 3.0*phi0_(i-1, j, k) - phi0_(i+1, j, k))
		phi3n = 0.125*(3.0*phi0_(i, j, k) + 6.0*phi0_(i-1, j, k) - phi0_(i-2, j, k))

		phi2p = 0.125*(3.0*phi0_(i, j, k) + 6.0*phi0_(i, j+1, k) - phi0_(i, j+2, k))
		phi2n = 0.125*(6.0*phi0_(i, j, k) + 3.0*phi0_(i, j+1, k) - phi0_(i, j-1, k))
		phi4p = 0.125*(6.0*phi0_(i, j, k) + 3.0*phi0_(i, j-1, k) - phi0_(i, j+1, k))
		phi4n = 0.125*(3.0*phi0_(i, j, k) + 6.0*phi0_(i, j-1, k) - phi0_(i, j-2, k))

		phi5p = 0.125*(3.0*phi0_(i, j, k) + 6.0*phi0_(i, j, k+1) - phi0_(i, j, k+2))
		phi5n = 0.125*(6.0*phi0_(i, j, k) + 3.0*phi0_(i, j, k+1) - phi0_(i, j, k-1))
		phi6p = 0.125*(6.0*phi0_(i, j, k) + 3.0*phi0_(i, j, k-1) - phi0_(i, j, k+1))
		phi6n = 0.125*(3.0*phi0_(i, j, k) + 6.0*phi0_(i, j, k-1) - phi0_(i, j, k-2))

		phi1 = lsm_getupwind(v0_(i, j, k, 1), phi1n, phi1p)
		phi3 = lsm_getupwind(v0_(i, j, k, 3), phi3n, phi3p)
		phi2 = lsm_getupwind(v0_(i, j, k, 2), phi2n, phi2p)
		phi4 = lsm_getupwind(v0_(i, j, k, 4), phi4n, phi4p)
		phi5 = lsm_getupwind(v0_(i, j, k, 5), phi5n, phi5p)
		phi6 = lsm_getupwind(v0_(i, j, k, 6), phi6n, phi6p)

		phi1_(i, j, k) = phi0_(i, j, k) &
											- (v0_(i, j, k, 1)*phi1 - v0_(i, j, k, 3)*phi3) &
											- (v0_(i, j, k, 2)*phi2 - v0_(i, j, k, 4)*phi4) &
											- (v0_(i, j, k, 5)*phi5 - v0_(i, j, k, 6)*phi6)
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine lsm_update_phi_v_2

subroutine lsm_update_phi_v_quick(phi1_, phi0_, v0_, psi_, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi1_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi0_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:6)	:: v0_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: psi_
	real										:: dpx_p2, dpy_p2, dpz_p2
	real										:: dpx_p1, dpy_p1, dpz_p1
	real										:: dpx_n1, dpy_n1, dpz_n1
	real										:: dpx_n2, dpy_n2, dpz_n2
	real										:: dpx_p, dpy_p, dpz_p
	real										:: dpx_c, dpy_c, dpz_c
	real										:: dpx_n, dpy_n, dpz_n
	real										:: phi1p, phi3p, phi2p, phi4p, phi5p, phi6p
	real										:: phi1n, phi3n, phi2n, phi4n, phi5n, phi6n
	real										:: phi0, phi1, phi3, phi2, phi4, phi5, phi6
	real										:: psi0, psi1, psi3, psi2, psi4, psi5, psi6
	real										:: dpsi0, dpsi1, dpsi3, dpsi2, dpsi4, dpsi5, dpsi6
	real										:: lsm_getminmod
	real										:: lsm_getupwind
	real										:: lsm_getd
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(dpx_p2, dpy_p2, dpz_p2) &
!$omp					,private(dpx_p1, dpy_p1, dpz_p1) &
!$omp					,private(dpx_n1, dpy_n1, dpz_n1) &
!$omp					,private(dpx_n2, dpy_n2, dpz_n2) &
!$omp					,private(dpx_p, dpy_p, dpz_p) &
!$omp					,private(dpx_c, dpy_c, dpz_c) &
!$omp					,private(dpx_n, dpy_n, dpz_n) &
!$omp					,private(phi1p, phi3p, phi2p, phi4p, phi5p, phi6p) &
!$omp					,private(phi1n, phi3n, phi2n, phi4n, phi5n, phi6n) &
!$omp					,private(phi0, phi1, phi3, phi2, phi4, phi5, phi6) &
!$omp					,private(psi0, psi1, psi3, psi2, psi4, psi5, psi6) &
!$omp					,private(dpsi0, dpsi1, dpsi3, dpsi2, dpsi4, dpsi5, dpsi6)
!$omp do schedule(dynamic, 1)
	do k=1, kx
	do j=1, jx
!ocl nouxsimd
	do i=1, ix
		phi1p = 0.125*(3.0*phi0_(i, j, k) + 6.0*phi0_(i+1, j, k) - phi0_(i+2, j, k))
		phi1n = 0.125*(6.0*phi0_(i, j, k) + 3.0*phi0_(i+1, j, k) - phi0_(i-1, j, k))
		phi3p = 0.125*(6.0*phi0_(i, j, k) + 3.0*phi0_(i-1, j, k) - phi0_(i+1, j, k))
		phi3n = 0.125*(3.0*phi0_(i, j, k) + 6.0*phi0_(i-1, j, k) - phi0_(i-2, j, k))

		phi2p = 0.125*(3.0*phi0_(i, j, k) + 6.0*phi0_(i, j+1, k) - phi0_(i, j+2, k))
		phi2n = 0.125*(6.0*phi0_(i, j, k) + 3.0*phi0_(i, j+1, k) - phi0_(i, j-1, k))
		phi4p = 0.125*(6.0*phi0_(i, j, k) + 3.0*phi0_(i, j-1, k) - phi0_(i, j+1, k))
		phi4n = 0.125*(3.0*phi0_(i, j, k) + 6.0*phi0_(i, j-1, k) - phi0_(i, j-2, k))

		phi5p = 0.125*(3.0*phi0_(i, j, k) + 6.0*phi0_(i, j, k+1) - phi0_(i, j, k+2))
		phi5n = 0.125*(6.0*phi0_(i, j, k) + 3.0*phi0_(i, j, k+1) - phi0_(i, j, k-1))
		phi6p = 0.125*(6.0*phi0_(i, j, k) + 3.0*phi0_(i, j, k-1) - phi0_(i, j, k+1))
		phi6n = 0.125*(3.0*phi0_(i, j, k) + 6.0*phi0_(i, j, k-1) - phi0_(i, j, k-2))

		phi1 = lsm_getupwind(v0_(i, j, k, 1), phi1n, phi1p)
		phi3 = lsm_getupwind(v0_(i, j, k, 3), phi3n, phi3p)
		phi2 = lsm_getupwind(v0_(i, j, k, 2), phi2n, phi2p)
		phi4 = lsm_getupwind(v0_(i, j, k, 4), phi4n, phi4p)
		phi5 = lsm_getupwind(v0_(i, j, k, 5), phi5n, phi5p)
		phi6 = lsm_getupwind(v0_(i, j, k, 6), phi6n, phi6p)

		psi0 = psi_(i, j, k)
		psi1 = psi_(i+1, j, k)
		psi3 = psi_(i-1, j, k)
		psi2 = psi_(i, j+1, k)
		psi4 = psi_(i, j-1, k)
		psi5 = psi_(i, j, k+1)
		psi6 = psi_(i, j, k-1)

    dpsi1 = lsm_getd(psi0, psi1)
    dpsi3 = lsm_getd(psi0, psi3)
    dpsi2 = lsm_getd(psi0, psi2)
    dpsi4 = lsm_getd(psi0, psi4)
    dpsi5 = lsm_getd(psi0, psi5)
    dpsi6 = lsm_getd(psi0, psi6)

		if( psi_(i+1, j, k) >= 0 ) then
			phi1 = 0.0
			phi1 = phi0_(i, j, k) + (dpsi1 - 0.5)/(dpsi1 + 0.5)*(phi0_(i, j, k) - phi0_(i-1, j, k))
			phi1 = phi0_(i, j, k)
		end if
		if( psi_(i-1, j, k) >= 0 ) then
			phi3 = 0.0
			phi3 = phi0_(i, j, k) - (dpsi3 - 0.5)/(dpsi3 + 0.5)*(phi0_(i+1, j, k) - phi0_(i, j, k))
			phi3 = phi0_(i, j, k)
		end if
		if( psi_(i, j+1, k) >= 0 ) then
			phi2 = 0.0
			phi2 = phi0_(i, j, k) + (dpsi2 - 0.5)/(dpsi2 + 0.5)*(phi0_(i, j, k) - phi0_(i, j-1, k))
			phi2 = phi0_(i, j, k)
		end if
		if( psi_(i, j-1, k) >= 0 ) then
			phi4 = 0.0
			phi4 = phi0_(i, j, k) - (dpsi4 - 0.5)/(dpsi4 + 0.5)*(phi0_(i, j+1, k) - phi0_(i, j, k))
			phi4 = phi0_(i, j, k)
		end if
		if( psi_(i, j, k+1) >= 0 ) then
			phi5 = 0.0
			phi5 = phi0_(i, j, k) + (dpsi5 - 0.5)/(dpsi5 + 0.5)*(phi0_(i, j, k) - phi0_(i, j, k-1))
			phi5 = phi0_(i, j, k)
		end if
		if( psi_(i, j, k-1) >= 0 ) then
			phi6 = 0.0
			phi6 = phi0_(i, j, k) - (dpsi6 - 0.5)/(dpsi6 + 0.5)*(phi0_(i, j, k+1) - phi0_(i, j, k))
			phi6 = phi0_(i, j, k)
		end if

		phi1_(i, j, k) = phi0_(i, j, k) &
											- (v0_(i, j, k, 1)*phi1 - v0_(i, j, k, 3)*phi3) &
											- (v0_(i, j, k, 2)*phi2 - v0_(i, j, k, 4)*phi4) &
											- (v0_(i, j, k, 5)*phi5 - v0_(i, j, k, 6)*phi6)
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine lsm_update_phi_v_quick

subroutine lsm_update_phi_v_weno3(phi1_, phi0_, v0_, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi1_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi0_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:6)	:: v0_
	real										:: dpx_p2, dpy_p2, dpz_p2
	real										:: dpx_p1, dpy_p1, dpz_p1
	real										:: dpx_n1, dpy_n1, dpz_n1
	real										:: dpx_n2, dpy_n2, dpz_n2
	real										:: dpx_p, dpy_p, dpz_p
	real										:: dpx_n, dpy_n, dpz_n
	real										:: vx, vy, vz
	real										:: lsm_getweno3
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(dpx_p2, dpy_p2, dpz_p2) &
!$omp					,private(dpx_p1, dpy_p1, dpz_p1) &
!$omp					,private(dpx_n1, dpy_n1, dpz_n1) &
!$omp					,private(dpx_n2, dpy_n2, dpz_n2) &
!$omp					,private(dpx_p, dpy_p, dpz_p) &
!$omp					,private(dpx_n, dpy_n, dpz_n) &
!$omp					,private(vx, vy, vz)
!$omp do schedule(dynamic, 1)
	do k=1, kx
	do j=1, jx
	do i=1, ix
		dpx_p2 = phi0_(i+2, j, k) - phi0_(i+1, j, k)
		dpx_p1 = phi0_(i+1, j, k) - phi0_(i  , j, k)
		dpx_n1 = phi0_(i  , j, k) - phi0_(i-1, j, k)
		dpx_n2 = phi0_(i-1, j, k) - phi0_(i-2, j, k)

		dpy_p2 = phi0_(i, j+2, k) - phi0_(i, j+1, k)
		dpy_p1 = phi0_(i, j+1, k) - phi0_(i, j  , k)
		dpy_n1 = phi0_(i, j  , k) - phi0_(i, j-1, k)
		dpy_n2 = phi0_(i, j-1, k) - phi0_(i, j-2, k)

		dpz_p2 = phi0_(i, j, k+2) - phi0_(i, j, k+1)
		dpz_p1 = phi0_(i, j, k+1) - phi0_(i, j, k  )
		dpz_n1 = phi0_(i, j, k  ) - phi0_(i, j, k-1)
		dpz_n2 = phi0_(i, j, k-1) - phi0_(i, j, k-2)

		dpx_p = lsm_getweno3(dpx_p2, dpx_p1, dpx_n1)
		dpx_n = lsm_getweno3(dpx_n2, dpx_n1, dpx_p1)
		dpy_p = lsm_getweno3(dpy_p2, dpy_p1, dpy_n1)
		dpy_n = lsm_getweno3(dpy_n2, dpy_n1, dpy_p1)
		dpz_p = lsm_getweno3(dpz_p2, dpz_p1, dpz_n1)
		dpz_n = lsm_getweno3(dpz_n2, dpz_n1, dpz_p1)

		vx = 0.5*(v0_(i, j, k, 1) + v0_(i, j, k, 3))
		vy = 0.5*(v0_(i, j, k, 2) + v0_(i, j, k, 4))
		vz = 0.5*(v0_(i, j, k, 5) + v0_(i, j, k, 6))

		phi1_(i, j, k) = phi0_(i, j, k) &
											- (max(vx, 0.0)*dpx_n + min(vx, 0.0)*dpx_p) &
											- (max(vy, 0.0)*dpy_n + min(vy, 0.0)*dpy_p) &
											- (max(vz, 0.0)*dpz_n + min(vz, 0.0)*dpz_p)
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine lsm_update_phi_v_weno3

subroutine lsm_update_phi_f(phi1_, phi0_, phif, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi1_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi0_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phif
	real										:: dpxx_n, dpxx_c, dpxx_p
	real										:: dpyy_n, dpyy_c, dpyy_p
	real										:: dpzz_n, dpzz_c, dpzz_p
	real										:: dpxx_nc, dpxx_pc
	real										:: dpyy_nc, dpyy_pc
	real										:: dpzz_nc, dpzz_pc
	real										:: dpx_p1, dpy_p1, dpz_p1
	real										:: dpx_n1, dpy_n1, dpz_n1
	real										:: dpx_p2, dpy_p2, dpz_p2
	real										:: dpx_n2, dpy_n2, dpz_n2
	real										:: dp_p, dp_n
	real										:: lsm_getminmod
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(dpxx_n, dpxx_c, dpxx_p) &
!$omp					,private(dpyy_n, dpyy_c, dpyy_p) &
!$omp					,private(dpzz_n, dpzz_c, dpzz_p) &
!$omp					,private(dpxx_nc, dpxx_pc) &
!$omp					,private(dpyy_nc, dpyy_pc) &
!$omp					,private(dpzz_nc, dpzz_pc) &
!$omp					,private(dpx_p1, dpy_p1, dpz_p1) &
!$omp					,private(dpx_n1, dpy_n1, dpz_n1) &
!$omp					,private(dpx_p2, dpy_p2, dpz_p2) &
!$omp					,private(dpx_n2, dpy_n2, dpz_n2) &
!$omp					,private(dp_p, dp_n)
!$omp do schedule(dynamic, 1)
	do k=1, kx
	do j=1, jx
!ocl nouxsimd
	do i=1, ix
		dpxx_n = phi0_(i+0, j, k) + phi0_(i-2, j, k) - 2.0*phi0_(i-1, j, k)
		dpxx_c = phi0_(i+1, j, k) + phi0_(i-1, j, k) - 2.0*phi0_(i  , j, k)
		dpxx_p = phi0_(i+2, j, k) + phi0_(i-0, j, k) - 2.0*phi0_(i+1, j, k)
		dpyy_n = phi0_(i, j+0, k) + phi0_(i, j-2, k) - 2.0*phi0_(i, j-1, k)
		dpyy_c = phi0_(i, j+1, k) + phi0_(i, j-1, k) - 2.0*phi0_(i, j  , k)
		dpyy_p = phi0_(i, j+2, k) + phi0_(i, j-0, k) - 2.0*phi0_(i, j+1, k)
		dpzz_n = phi0_(i, j, k+0) + phi0_(i, j, k-2) - 2.0*phi0_(i, j, k-1)
		dpzz_c = phi0_(i, j, k+1) + phi0_(i, j, k-1) - 2.0*phi0_(i, j, k  )
		dpzz_p = phi0_(i, j, k+2) + phi0_(i, j, k-0) - 2.0*phi0_(i, j, k+1)
		dpx_n1 = phi0_(i, j, k) - phi0_(i-1, j, k)
		dpy_n1 = phi0_(i, j, k) - phi0_(i, j-1, k)
		dpz_n1 = phi0_(i, j, k) - phi0_(i, j, k-1)
		dpx_p1 = phi0_(i+1, j, k) - phi0_(i, j, k)
		dpy_p1 = phi0_(i, j+1, k) - phi0_(i, j, k)
		dpz_p1 = phi0_(i, j, k+1) - phi0_(i, j, k)
		dpx_n2 = dpx_n1 + 0.5*lsm_getminmod(dpxx_c, dpxx_n)
		dpy_n2 = dpy_n1 + 0.5*lsm_getminmod(dpyy_c, dpyy_n)
		dpz_n2 = dpz_n1 + 0.5*lsm_getminmod(dpzz_c, dpzz_n)
		dpx_p2 = dpx_p1 - 0.5*lsm_getminmod(dpxx_c, dpxx_p)
		dpy_p2 = dpy_p1 - 0.5*lsm_getminmod(dpyy_c, dpyy_p)
		dpz_p2 = dpz_p1 - 0.5*lsm_getminmod(dpzz_c, dpzz_p)
		dp_p = sqrt( &
								max(dpx_n2, 0.0)*max(dpx_n2, 0.0) &
							+ min(dpx_p2, 0.0)*min(dpx_p2, 0.0) &
							+ max(dpy_n2, 0.0)*max(dpy_n2, 0.0) &
							+ min(dpy_p2, 0.0)*min(dpy_p2, 0.0) &
							+ max(dpz_n2, 0.0)*max(dpz_n2, 0.0) &
							+ min(dpz_p2, 0.0)*min(dpz_p2, 0.0) &
						)
		dp_n = sqrt( &
								min(dpx_n2, 0.0)*min(dpx_n2, 0.0) &
							+ max(dpx_p2, 0.0)*max(dpx_p2, 0.0) &
							+ min(dpy_n2, 0.0)*min(dpy_n2, 0.0) &
							+ max(dpy_p2, 0.0)*max(dpy_p2, 0.0) &
							+ min(dpz_n2, 0.0)*min(dpz_n2, 0.0) &
							+ max(dpz_p2, 0.0)*max(dpz_p2, 0.0) &
						)
		phi1_(i, j, k) = phi0_(i, j, k) &
											- ( max(phif(i, j, k), 0.0)*dp_p &
												+ min(phif(i, j, k), 0.0)*dp_n )
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine lsm_update_phi_f

subroutine lsm_update_phi_f_2(phi1_, phi0_, phif, psi, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi1_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi0_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phif
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: psi
	real										:: dpxx_n, dpxx_c, dpxx_p
	real										:: dpyy_n, dpyy_c, dpyy_p
	real										:: dpzz_n, dpzz_c, dpzz_p
	real										:: dpxx_nc, dpxx_pc
	real										:: dpyy_nc, dpyy_pc
	real										:: dpzz_nc, dpzz_pc
	real										:: dpx_p1, dpy_p1, dpz_p1
	real										:: dpx_n1, dpy_n1, dpz_n1
	real										:: dpx_p2, dpy_p2, dpz_p2
	real										:: dpx_n2, dpy_n2, dpz_n2
	real										:: dp_p, dp_n
	real										:: lsm_getminmod
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(dpxx_n, dpxx_c, dpxx_p) &
!$omp					,private(dpyy_n, dpyy_c, dpyy_p) &
!$omp					,private(dpzz_n, dpzz_c, dpzz_p) &
!$omp					,private(dpxx_nc, dpxx_pc) &
!$omp					,private(dpyy_nc, dpyy_pc) &
!$omp					,private(dpzz_nc, dpzz_pc) &
!$omp					,private(dpx_p1, dpy_p1, dpz_p1) &
!$omp					,private(dpx_n1, dpy_n1, dpz_n1) &
!$omp					,private(dpx_p2, dpy_p2, dpz_p2) &
!$omp					,private(dpx_n2, dpy_n2, dpz_n2) &
!$omp					,private(dp_p, dp_n)
!$omp do schedule(dynamic, 1)
	do k=1, kx
	do j=1, jx
!ocl nouxsimd
	do i=1, ix
		dpxx_n = phi0_(i+0, j, k) + phi0_(i-2, j, k) - 2.0*phi0_(i-1, j, k)
		dpxx_c = phi0_(i+1, j, k) + phi0_(i-1, j, k) - 2.0*phi0_(i  , j, k)
		dpxx_p = phi0_(i+2, j, k) + phi0_(i-0, j, k) - 2.0*phi0_(i+1, j, k)
		dpyy_n = phi0_(i, j+0, k) + phi0_(i, j-2, k) - 2.0*phi0_(i, j-1, k)
		dpyy_c = phi0_(i, j+1, k) + phi0_(i, j-1, k) - 2.0*phi0_(i, j  , k)
		dpyy_p = phi0_(i, j+2, k) + phi0_(i, j-0, k) - 2.0*phi0_(i, j+1, k)
		dpzz_n = phi0_(i, j, k+0) + phi0_(i, j, k-2) - 2.0*phi0_(i, j, k-1)
		dpzz_c = phi0_(i, j, k+1) + phi0_(i, j, k-1) - 2.0*phi0_(i, j, k  )
		dpzz_p = phi0_(i, j, k+2) + phi0_(i, j, k-0) - 2.0*phi0_(i, j, k+1)
		dpx_n1 = phi0_(i, j, k) - phi0_(i-1, j, k)
		dpy_n1 = phi0_(i, j, k) - phi0_(i, j-1, k)
		dpz_n1 = phi0_(i, j, k) - phi0_(i, j, k-1)
		dpx_p1 = phi0_(i+1, j, k) - phi0_(i, j, k)
		dpy_p1 = phi0_(i, j+1, k) - phi0_(i, j, k)
		dpz_p1 = phi0_(i, j, k+1) - phi0_(i, j, k)
		dpx_n2 = dpx_n1 + 0.5*lsm_getminmod(dpxx_c, dpxx_n)
		dpy_n2 = dpy_n1 + 0.5*lsm_getminmod(dpyy_c, dpyy_n)
		dpz_n2 = dpz_n1 + 0.5*lsm_getminmod(dpzz_c, dpzz_n)
		dpx_p2 = dpx_p1 - 0.5*lsm_getminmod(dpxx_c, dpxx_p)
		dpy_p2 = dpy_p1 - 0.5*lsm_getminmod(dpyy_c, dpyy_p)
		dpz_p2 = dpz_p1 - 0.5*lsm_getminmod(dpzz_c, dpzz_p)

		if( psi(i+1, j, k) >= 0 ) then
			dpx_p2 = dpx_n2
		endif
		if( psi(i-1, j, k) >= 0 ) then
			dpx_n2 = dpx_p2
		endif
		if( psi(i, j+1, k) >= 0 ) then
			dpy_p2 = dpy_n2
		endif
		if( psi(i, j-1, k) >= 0 ) then
			dpy_n2 = dpy_p2
		endif
		if( psi(i, j, k+1) >= 0 ) then
			dpz_p2 = dpz_n2
		endif
		if( psi(i, j, k-1) >= 0 ) then
			dpz_n2 = dpz_p2
		endif

		dp_p = sqrt( &
								max(dpx_n2, 0.0)*max(dpx_n2, 0.0) &
							+ min(dpx_p2, 0.0)*min(dpx_p2, 0.0) &
							+ max(dpy_n2, 0.0)*max(dpy_n2, 0.0) &
							+ min(dpy_p2, 0.0)*min(dpy_p2, 0.0) &
							+ max(dpz_n2, 0.0)*max(dpz_n2, 0.0) &
							+ min(dpz_p2, 0.0)*min(dpz_p2, 0.0) &
						)
		dp_n = sqrt( &
								min(dpx_n2, 0.0)*min(dpx_n2, 0.0) &
							+ max(dpx_p2, 0.0)*max(dpx_p2, 0.0) &
							+ min(dpy_n2, 0.0)*min(dpy_n2, 0.0) &
							+ max(dpy_p2, 0.0)*max(dpy_p2, 0.0) &
							+ min(dpz_n2, 0.0)*min(dpz_n2, 0.0) &
							+ max(dpz_p2, 0.0)*max(dpz_p2, 0.0) &
						)
		phi1_(i, j, k) = phi0_(i, j, k) &
											- ( max(phif(i, j, k), 0.0)*dp_p &
												+ min(phif(i, j, k), 0.0)*dp_n )
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine lsm_update_phi_f_2

subroutine lsm_calc_signeddistance(d, phi, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: d
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi
	real										:: dpx_p1, dpy_p1, dpz_p1
	real										:: dpx_n1, dpy_n1, dpz_n1
	real										:: dpx_p, dpy_p, dpz_p
	real										:: dpx_c, dpy_c, dpz_c
	real										:: dpx_n, dpy_n, dpz_n
	real										:: dphi1, dphi2, dphi3, dphi4, dphi5, dphi6, dphi7, dphi
	real										:: lsm_getsign
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(dpx_p1, dpy_p1, dpz_p1) &
!$omp					,private(dpx_n1, dpy_n1, dpz_n1) &
!$omp					,private(dpx_p, dpy_p, dpz_p) &
!$omp					,private(dpx_c, dpy_c, dpz_c) &
!$omp					,private(dpx_n, dpy_n, dpz_n) &
!$omp					,private(dphi1, dphi2, dphi3, dphi4, dphi5, dphi6, dphi7, dphi)
!$omp do schedule(dynamic, 1)
	do k=1, kx
	do j=1, jx
!ocl nouxsimd
	do i=1, ix
		dpx_p1 = phi(i+1, j, k) - phi(i  , j, k)
		dpx_n1 = phi(i  , j, k) - phi(i-1, j, k)
		dpy_p1 = phi(i, j+1, k) - phi(i, j  , k)
		dpy_n1 = phi(i, j  , k) - phi(i, j-1, k)
		dpz_p1 = phi(i, j, k+1) - phi(i, j, k  )
		dpz_n1 = phi(i, j, k  ) - phi(i, j, k-1)

		dpx_c = 0.0
		dpx_c = 0.0
		dpy_c = 0.0
		dpy_c = 0.0
		dpz_c = 0.0
		dpz_c = 0.0
		if( phi(i, j, k)*phi(i+1, j, k) < 0 ) then
			dpx_c = phi(i+1, j, k) - phi(i  , j, k)
		end if
		if( phi(i, j, k)*phi(i-1, j, k) < 0 ) then
			dpx_c = phi(i  , j, k) - phi(i-1, j, k)
		end if
		if( phi(i, j, k)*phi(i, j+1, k) < 0 ) then
			dpy_c = phi(i, j+1, k) - phi(i, j  , k)
		end if
		if( phi(i, j, k)*phi(i, j-1, k) < 0 ) then
			dpy_c = phi(i, j  , k) - phi(i, j-1, k)
		end if
		if( phi(i, j, k)*phi(i, j, k+1) < 0 ) then
			dpz_c = phi(i, j, k+1) - phi(i, j, k  )
		end if
		if( phi(i, j, k)*phi(i, j, k-1) < 0 ) then
			dpz_c = phi(i, j, k  ) - phi(i, j, k-1)
		end if

		if( phi(i, j, k)*phi(i+1, j, k) < 0 .and. phi(i, j, k)*phi(i-1, j, k) < 0 ) then
			dpx_c = 0.5*(phi(i+1, j, k) - phi(i-1, j, k))
		end if
		if( phi(i, j, k)*phi(i, j+1, k) < 0 .and. phi(i, j, k)*phi(i, j-1, k) < 0 ) then
			dpy_c = 0.5*(phi(i, j+1, k) - phi(i, j-1, k))
		end if
		if( phi(i, j, k)*phi(i, j, k+1) < 0 .and. phi(i, j, k)*phi(i, j, k-1) < 0 ) then
			dpz_c = 0.5*(phi(i, j, k+1) - phi(i, j, k-1))
		end if

		dpx_p = dpx_p1
		dpx_n = dpx_n1
		dpy_p = dpy_p1
		dpy_n = dpy_n1
		dpz_p = dpz_p1
		dpz_n = dpz_n1

!		dpx_c = 0.5*(dpx_p1 + dpx_n1)
!		dpy_c = 0.5*(dpy_p1 + dpy_n1)
!		dpz_c = 0.5*(dpz_p1 + dpz_n1)

		dphi1 = max(abs(dpx_p), abs(dpx_n))
		dphi2 = max(abs(dpy_p), abs(dpy_n))
		dphi3 = max(abs(dpz_p), abs(dpz_n))
		dphi4 = max(dphi1, abs(dpx_c))
		dphi5 = max(dphi2, abs(dpy_c))
		dphi6 = max(dphi3, abs(dpz_c))
		dphi7 = sqrt(dpx_c*dpx_c + dpy_c*dpy_c + dpz_c*dpz_c) 
		dphi1 = max(dphi4, dphi5)
		dphi2 = max(dphi6, dphi7)
		dphi3 = max(dphi1, dphi2)
		dphi = max(dphi7, 1.0)

		d(i, j, k) = phi(i, j, k)
		d(i, j, k) = phi(i, j, k)/(dphi7 + 1.0e-6)
		d(i, j, k) = phi(i, j, k)/dphi
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine lsm_calc_signeddistance

subroutine lsm_calc_signeddistance2(d, phi, cid, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: d
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: cid
	real										:: dpx_p1, dpy_p1, dpz_p1
	real										:: dpx_n1, dpy_n1, dpz_n1
	real										:: dpx_p, dpy_p, dpz_p
	real										:: dpx_c, dpy_c, dpz_c
	real										:: dpx_n, dpy_n, dpz_n
	real										:: dphi1, dphi2, dphi3, dphi4, dphi5, dphi6, dphi7, dphi
	real										:: lsm_getsign
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(dpx_p1, dpy_p1, dpz_p1) &
!$omp					,private(dpx_n1, dpy_n1, dpz_n1) &
!$omp					,private(dpx_p, dpy_p, dpz_p) &
!$omp					,private(dpx_c, dpy_c, dpz_c) &
!$omp					,private(dpx_n, dpy_n, dpz_n) &
!$omp					,private(dphi1, dphi2, dphi3, dphi4, dphi5, dphi6, dphi7, dphi)
!$omp do schedule(dynamic, 1)
	do k=1, kx
	do j=1, jx
!ocl nouxsimd
	do i=1, ix
		cid(i, j, k) = -1.0

		if( phi(i, j, k)*phi(i+1, j, k) < 0 .or. &
				phi(i, j, k)*phi(i-1, j, k) < 0 .or. &
				phi(i, j, k)*phi(i, j+1, k) < 0 .or. &
				phi(i, j, k)*phi(i, j-1, k) < 0 .or. &
				phi(i, j, k)*phi(i, j, k+1) < 0 .or. &
				phi(i, j, k)*phi(i, j, k-1) < 0 ) then
			cid(i, j, k) = 1.0
		end if
	end do
	end do
	end do
!$omp end do

!$omp do schedule(dynamic, 1)
	do k=1, kx
	do j=1, jx
!ocl nouxsimd
	do i=1, ix
		d(i, j, k) = 0.0

		if( cid(i, j, k) < 0 ) then
			cycle
		end if

		dpx_c = 0.0
		dpx_c = 0.0
		dpy_c = 0.0
		dpy_c = 0.0
		dpz_c = 0.0
		dpz_c = 0.0
		if( cid(i+1, j, k) > 0 ) then
			dpx_c = phi(i+1, j, k) - phi(i  , j, k)
		end if
		if( cid(i-1, j, k) > 0 ) then
			dpx_c = phi(i  , j, k) - phi(i-1, j, k)
		end if
		if( cid(i, j+1, k) > 0 ) then
			dpy_c = phi(i, j+1, k) - phi(i, j  , k)
		end if
		if( cid(i, j-1, k) > 0 ) then
			dpy_c = phi(i, j  , k) - phi(i, j-1, k)
		end if
		if( cid(i, j, k+1) > 0 ) then
			dpz_c = phi(i, j, k+1) - phi(i, j, k  )
		end if
		if( cid(i, j, k-1) > 0 ) then
			dpz_c = phi(i, j, k  ) - phi(i, j, k-1)
		end if

		if( cid(i+1, j, k) > 0 .and. cid(i-1, j, k) > 0 ) then
			dpx_c = 0.5*(phi(i+1, j, k) - phi(i-1, j, k))
		end if
		if( cid(i, j+1, k) > 0 .and. cid(i, j-1, k) > 0 ) then
			dpy_c = 0.5*(phi(i, j+1, k) - phi(i, j-1, k))
		end if
		if( cid(i, j, k+1) > 0 .and. cid(i, j, k-1) > 0 ) then
			dpz_c = 0.5*(phi(i, j, k+1) - phi(i, j, k-1))
		end if

		dphi7 = sqrt(dpx_c*dpx_c + dpy_c*dpy_c + dpz_c*dpz_c) 
		dphi4 = max(dphi7, abs(dpx_c))
		dphi5 = max(dphi7, abs(dpy_c))
		dphi6 = max(dphi7, abs(dpz_c))
		dphi1 = max(dphi4, dphi5)
		dphi2 = max(dphi1, dphi6)
		dphi3 = max(dphi2, dphi7)
		dphi = max(dphi3, 1.0)

		dphi7 = sqrt(dpx_c*dpx_c + dpy_c*dpy_c + dpz_c*dpz_c) 
		dphi = min(dphi7, 1.0)
		dphi = max(dphi7, 1.0)

		d(i, j, k) = phi(i, j, k)/dphi
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine lsm_calc_signeddistance2

subroutine lsm_reinit_core(phi1_, phi0_, dtau, phii_, phid_, psi, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi1_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi0_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phii_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phid_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: psi
	real										:: dtau
	real										:: dpxx_p, dpyy_p, dpzz_p
	real										:: dpxx_n, dpyy_n, dpzz_n
	real										:: dpxx_c, dpyy_c, dpzz_c
	real										:: dpx_p1, dpy_p1, dpz_p1
	real										:: dpx_n1, dpy_n1, dpz_n1
	real										:: dpx_p2, dpy_p2, dpz_p2
	real										:: dpx_n2, dpy_n2, dpz_n2
	real										:: dpx_p_max2, dpy_p_max2, dpz_p_max2
	real										:: dpx_n_max2, dpy_n_max2, dpz_n_max2
	real										:: dpx_p_min2, dpy_p_min2, dpz_p_min2
	real										:: dpx_n_min2, dpy_n_min2, dpz_n_min2
	real										:: dp_p, dp_n
	real										:: dpl
	real										:: phis
	real										:: lsm_getsign
	real										:: lsm_getminmod
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(dpxx_p, dpyy_p, dpzz_p) &
!$omp					,private(dpxx_n, dpyy_n, dpzz_n) &
!$omp					,private(dpxx_c, dpyy_c, dpzz_c) &
!$omp					,private(dpx_p1, dpy_p1, dpz_p1) &
!$omp					,private(dpx_n1, dpy_n1, dpz_n1) &
!$omp					,private(dpx_p2, dpy_p2, dpz_p2) &
!$omp					,private(dpx_n2, dpy_n2, dpz_n2) &
!$omp					,private(dpx_p_max2, dpy_p_max2, dpz_p_max2) &
!$omp					,private(dpx_n_max2, dpy_n_max2, dpz_n_max2) &
!$omp					,private(dpx_p_min2, dpy_p_min2, dpz_p_min2) &
!$omp					,private(dpx_n_min2, dpy_n_min2, dpz_n_min2) &
!$omp					,private(dp_p, dp_n) &
!$omp					,private(dpl) &
!$omp					,private(phis)
!$omp do schedule(dynamic, 1)
	do k=1, kx
	do j=1, jx
	do i=1, ix
!		if( psi(i, j, k) >= 0 ) then
!			phi1_(i, j, k) = phi0_(i, j, k)
!			cycle
!		endif

		dpx_n1 = phi0_(i, j, k) - phi0_(i-1, j, k)
		dpy_n1 = phi0_(i, j, k) - phi0_(i, j-1, k)
		dpz_n1 = phi0_(i, j, k) - phi0_(i, j, k-1)
		dpx_p1 = phi0_(i+1, j, k) - phi0_(i, j, k)
		dpy_p1 = phi0_(i, j+1, k) - phi0_(i, j, k)
		dpz_p1 = phi0_(i, j, k+1) - phi0_(i, j, k)
		dpxx_n = phi0_(i+0, j, k) + phi0_(i-2, j, k) - 2.0*phi0_(i-1, j, k)
		dpxx_c = phi0_(i+1, j, k) + phi0_(i-1, j, k) - 2.0*phi0_(i  , j, k)
		dpxx_p = phi0_(i+2, j, k) + phi0_(i-0, j, k) - 2.0*phi0_(i+1, j, k)
		dpyy_n = phi0_(i, j+0, k) + phi0_(i, j-2, k) - 2.0*phi0_(i, j-1, k)
		dpyy_c = phi0_(i, j+1, k) + phi0_(i, j-1, k) - 2.0*phi0_(i, j  , k)
		dpyy_p = phi0_(i, j+2, k) + phi0_(i, j-0, k) - 2.0*phi0_(i, j+1, k)
		dpzz_n = phi0_(i, j, k+0) + phi0_(i, j, k-2) - 2.0*phi0_(i, j, k-1)
		dpzz_c = phi0_(i, j, k+1) + phi0_(i, j, k-1) - 2.0*phi0_(i, j, k  )
		dpzz_p = phi0_(i, j, k+2) + phi0_(i, j, k-0) - 2.0*phi0_(i, j, k+1)
		dpx_n2 = dpx_n1 + 0.5*lsm_getminmod(dpxx_c, dpxx_n)
		dpy_n2 = dpy_n1 + 0.5*lsm_getminmod(dpyy_c, dpyy_n)
		dpz_n2 = dpz_n1 + 0.5*lsm_getminmod(dpzz_c, dpzz_n)
		dpx_p2 = dpx_p1 - 0.5*lsm_getminmod(dpxx_c, dpxx_p)
		dpy_p2 = dpy_p1 - 0.5*lsm_getminmod(dpyy_c, dpyy_p)
		dpz_p2 = dpz_p1 - 0.5*lsm_getminmod(dpzz_c, dpzz_p)

		dpx_n_max2 = max(dpx_n2, 0.0)*max(dpx_n2, 0.0)
		dpx_n_min2 = min(dpx_n2, 0.0)*min(dpx_n2, 0.0)
		dpx_p_max2 = max(dpx_p2, 0.0)*max(dpx_p2, 0.0)
		dpx_p_min2 = min(dpx_p2, 0.0)*min(dpx_p2, 0.0)
		dpy_n_max2 = max(dpy_n2, 0.0)*max(dpy_n2, 0.0)
		dpy_n_min2 = min(dpy_n2, 0.0)*min(dpy_n2, 0.0)
		dpy_p_max2 = max(dpy_p2, 0.0)*max(dpy_p2, 0.0)
		dpy_p_min2 = min(dpy_p2, 0.0)*min(dpy_p2, 0.0)
		dpz_n_max2 = max(dpz_n2, 0.0)*max(dpz_n2, 0.0)
		dpz_n_min2 = min(dpz_n2, 0.0)*min(dpz_n2, 0.0)
		dpz_p_max2 = max(dpz_p2, 0.0)*max(dpz_p2, 0.0)
		dpz_p_min2 = min(dpz_p2, 0.0)*min(dpz_p2, 0.0)

		dpx_p2 = max(dpx_n_max2, dpx_p_min2)
		dpy_p2 = max(dpy_n_max2, dpy_p_min2)
		dpz_p2 = max(dpz_n_max2, dpz_p_min2)
		dpx_n2 = max(dpx_n_min2, dpx_p_max2)
		dpy_n2 = max(dpy_n_min2, dpy_p_max2)
		dpz_n2 = max(dpz_n_min2, dpz_p_max2)

		dp_p = sqrt(dpx_p2 + dpy_p2 + dpz_p2)
		dp_n = sqrt(dpx_n2 + dpy_n2 + dpz_n2)

		phis = lsm_getsign( phii_(i, j, k) )

		if( phis > 0.0 ) then
			dpl = dp_p
		else if( phis < 0.0 ) then
			dpl = dp_n
		else
			dpl = 1.0
		end if

		if( phii_(i, j, k)*phii_(i+1, j, k) < 0 .or. &
				phii_(i, j, k)*phii_(i-1, j, k) < 0 .or. &
				phii_(i, j, k)*phii_(i, j+1, k) < 0 .or. &
				phii_(i, j, k)*phii_(i, j-1, k) < 0 .or. &
				phii_(i, j, k)*phii_(i, j, k+1) < 0 .or. &
				phii_(i, j, k)*phii_(i, j, k-1) < 0 ) then
			phi1_(i, j, k) = phid_(i, j, k)
		else
			phi1_(i, j, k) = phi0_(i, j, k) - phis*(dpl - 1.0)*dtau
		endif
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine lsm_reinit_core

subroutine lsm_reinit_s_core(phi1_, phi0_, dtau, psi, psinx, psiny, psinz, thetas, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi1_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi0_
	real										:: dtau
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: psi
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: psinx, psiny, psinz
	real										:: thetas
	real										:: wx, wy, wz
	real										:: dpxx_p, dpyy_p, dpzz_p
	real										:: dpxx_n, dpyy_n, dpzz_n
	real										:: dpxx_c, dpyy_c, dpzz_c
	real										:: dpx_p1, dpy_p1, dpz_p1
	real										:: dpx_n1, dpy_n1, dpz_n1
	real										:: dpx_p2, dpy_p2, dpz_p2
	real										:: dpx_n2, dpy_n2, dpz_n2
	real										:: dpx_p, dpy_p, dpz_p
	real										:: dpx_n, dpy_n, dpz_n
	real										:: dpx, dpy, dpz
	real										:: dpl
	real										:: lsm_getupwind2
	real										:: lsm_getupwind1
	real										:: lsm_getupwind
	real										:: lsm_getweno3
	real										:: lsm_getminmod
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(wx, wy, wz) &
!$omp					,private(dpxx_p, dpyy_p, dpzz_p) &
!$omp					,private(dpxx_n, dpyy_n, dpzz_n) &
!$omp					,private(dpxx_c, dpyy_c, dpzz_c) &
!$omp					,private(dpx_p1, dpy_p1, dpz_p1) &
!$omp					,private(dpx_n1, dpy_n1, dpz_n1) &
!$omp					,private(dpx_p2, dpy_p2, dpz_p2) &
!$omp					,private(dpx_n2, dpy_n2, dpz_n2) &
!$omp					,private(dpx_p, dpy_p, dpz_p) &
!$omp					,private(dpx_n, dpy_n, dpz_n) &
!$omp					,private(dpx, dpy, dpz) &
!$omp					,private(dpl) 
!$omp do schedule(static, 1)
	do k=1, kx
	do j=1, jx
	do i=1, ix
		if( psi(i, j, k) < 0 ) then
			phi1_(i, j, k) = phi0_(i, j, k)
			cycle
		endif

		if( abs(psi(i, j, k)) > 4.0*sqrt(3.0) ) then
!			cycle
		endif

		if( abs(phi0_(i, j, k)) > 4.0*sqrt(3.0) ) then
!			cycle
		endif

		wx = psinx(i, j, k)
		wy = psiny(i, j, k)
		wz = psinz(i, j, k)

		dpx_n1 = phi0_(i, j, k) - phi0_(i-1, j, k)
		dpy_n1 = phi0_(i, j, k) - phi0_(i, j-1, k)
		dpz_n1 = phi0_(i, j, k) - phi0_(i, j, k-1)
		dpx_p1 = phi0_(i+1, j, k) - phi0_(i, j, k)
		dpy_p1 = phi0_(i, j+1, k) - phi0_(i, j, k)
		dpz_p1 = phi0_(i, j, k+1) - phi0_(i, j, k)
		dpxx_n = phi0_(i+0, j, k) + phi0_(i-2, j, k) - 2.0*phi0_(i-1, j, k)
		dpxx_c = phi0_(i+1, j, k) + phi0_(i-1, j, k) - 2.0*phi0_(i  , j, k)
		dpxx_p = phi0_(i+2, j, k) + phi0_(i-0, j, k) - 2.0*phi0_(i+1, j, k)
		dpyy_n = phi0_(i, j+0, k) + phi0_(i, j-2, k) - 2.0*phi0_(i, j-1, k)
		dpyy_c = phi0_(i, j+1, k) + phi0_(i, j-1, k) - 2.0*phi0_(i, j  , k)
		dpyy_p = phi0_(i, j+2, k) + phi0_(i, j-0, k) - 2.0*phi0_(i, j+1, k)
		dpzz_n = phi0_(i, j, k+0) + phi0_(i, j, k-2) - 2.0*phi0_(i, j, k-1)
		dpzz_c = phi0_(i, j, k+1) + phi0_(i, j, k-1) - 2.0*phi0_(i, j, k  )
		dpzz_p = phi0_(i, j, k+2) + phi0_(i, j, k-0) - 2.0*phi0_(i, j, k+1)
		dpx_n2 = dpx_n1 + 0.5*lsm_getminmod(dpxx_c, dpxx_n)
		dpy_n2 = dpy_n1 + 0.5*lsm_getminmod(dpyy_c, dpyy_n)
		dpz_n2 = dpz_n1 + 0.5*lsm_getminmod(dpzz_c, dpzz_n)
		dpx_p2 = dpx_p1 - 0.5*lsm_getminmod(dpxx_c, dpxx_p)
		dpy_p2 = dpy_p1 - 0.5*lsm_getminmod(dpyy_c, dpyy_p)
		dpz_p2 = dpz_p1 - 0.5*lsm_getminmod(dpzz_c, dpzz_p)
		dpx_p = dpx_p2
		dpx_n = dpx_n2
		dpy_p = dpy_p2
		dpy_n = dpy_n2
		dpz_p = dpz_p2
		dpz_n = dpz_n2

		phi1_(i, j, k) = phi0_(i, j, k) -( &
																			+ (max(wx, 0.0)*dpx_n + min(wx, 0.0)*dpx_p) &
																			+ (max(wy, 0.0)*dpy_n + min(wy, 0.0)*dpy_p) &
																			+ (max(wz, 0.0)*dpz_n + min(wz, 0.0)*dpz_p) - cos(thetas) )*dtau
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine lsm_reinit_s_core

subroutine lsm_reinit_core_0(phi1_, phi0_, dtau, phii_, phid_, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi1_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi0_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phii_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phid_
	real										:: dtau
	real										:: dpx_p1, dpy_p1, dpz_p1
	real										:: dpx_n1, dpy_n1, dpz_n1
	real										:: dpx_p_max2, dpy_p_max2, dpz_p_max2
	real										:: dpx_n_max2, dpy_n_max2, dpz_n_max2
	real										:: dpx_p_min2, dpy_p_min2, dpz_p_min2
	real										:: dpx_n_min2, dpy_n_min2, dpz_n_min2
	real										:: dpx_p2, dpy_p2, dpz_p2
	real										:: dpx_n2, dpy_n2, dpz_n2
	real										:: dp_p, dp_n
	real										:: dpl
	real										:: phis
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(dpx_p1, dpy_p1, dpz_p1) &
!$omp					,private(dpx_n1, dpy_n1, dpz_n1) &
!$omp					,private(dpx_p_max2, dpy_p_max2, dpz_p_max2) &
!$omp					,private(dpx_n_max2, dpy_n_max2, dpz_n_max2) &
!$omp					,private(dpx_p_min2, dpy_p_min2, dpz_p_min2) &
!$omp					,private(dpx_n_min2, dpy_n_min2, dpz_n_min2) &
!$omp					,private(dpx_p2, dpy_p2, dpz_p2) &
!$omp					,private(dpx_n2, dpy_n2, dpz_n2) &
!$omp					,private(dp_p, dp_n) &
!$omp					,private(dpl) &
!$omp					,private(phis)
!$omp do schedule(dynamic, 1)
	do k=1, kx
	do j=1, jx
	do i=1, ix
		dpx_p1 = phi0_(i+1, j, k) - phi0_(i  , j, k)
		dpx_n1 = phi0_(i  , j, k) - phi0_(i-1, j, k)
		dpy_p1 = phi0_(i, j+1, k) - phi0_(i, j  , k)
		dpy_n1 = phi0_(i, j  , k) - phi0_(i, j-1, k)
		dpz_p1 = phi0_(i, j, k+1) - phi0_(i, j, k  )
		dpz_n1 = phi0_(i, j, k  ) - phi0_(i, j, k-1)

		dpx_n_max2 = max(dpx_n1, 0.0)*max(dpx_n1, 0.0)
		dpx_n_min2 = min(dpx_n1, 0.0)*min(dpx_n1, 0.0)
		dpx_p_max2 = max(dpx_p1, 0.0)*max(dpx_p1, 0.0)
		dpx_p_min2 = min(dpx_p1, 0.0)*min(dpx_p1, 0.0)
		dpy_n_max2 = max(dpy_n1, 0.0)*max(dpy_n1, 0.0)
		dpy_n_min2 = min(dpy_n1, 0.0)*min(dpy_n1, 0.0)
		dpy_p_max2 = max(dpy_p1, 0.0)*max(dpy_p1, 0.0)
		dpy_p_min2 = min(dpy_p1, 0.0)*min(dpy_p1, 0.0)
		dpz_n_max2 = max(dpz_n1, 0.0)*max(dpz_n1, 0.0)
		dpz_n_min2 = min(dpz_n1, 0.0)*min(dpz_n1, 0.0)
		dpz_p_max2 = max(dpz_p1, 0.0)*max(dpz_p1, 0.0)
		dpz_p_min2 = min(dpz_p1, 0.0)*min(dpz_p1, 0.0)

		dpx_p2 = max(dpx_n_max2, dpx_p_min2)
		dpy_p2 = max(dpy_n_max2, dpy_p_min2)
		dpz_p2 = max(dpz_n_max2, dpz_p_min2)
		dpx_n2 = max(dpx_n_min2, dpx_p_max2)
		dpy_n2 = max(dpy_n_min2, dpy_p_max2)
		dpz_n2 = max(dpz_n_min2, dpz_p_max2)

		dp_p = sqrt(dpx_p2 + dpy_p2 + dpz_p2)
		dp_n = sqrt(dpx_n2 + dpy_n2 + dpz_n2)

		if(phii_(i, j, k) > 0.0) then
			dpl = dp_p
		else if(phii_(i, j, k) < 0.0) then
			dpl = dp_n
		else
			dpl = 1.0
		endif

		if(phii_(i, j, k) > 0.0) then
			phis = 1.0
		else if(phii_(i, j, k) < 0.0) then
			phis = -1.0
		else
			phis = 0.0
		endif

		if( phii_(i, j, k)*phii_(i+1, j, k) < 0 .or. &
				phii_(i, j, k)*phii_(i-1, j, k) < 0 .or. &
				phii_(i, j, k)*phii_(i, j+1, k) < 0 .or. &
				phii_(i, j, k)*phii_(i, j-1, k) < 0 .or. &
				phii_(i, j, k)*phii_(i, j, k+1) < 0 .or. &
				phii_(i, j, k)*phii_(i, j, k-1) < 0 ) then
			phi1_(i, j, k) = phi0_(i, j, k)*(1.0 - dtau) + phid_(i, j, k)*dtau
		else
			phi1_(i, j, k) = phi0_(i, j, k) - phis*(dpl - 1.0)*dtau
		endif
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine lsm_reinit_core_0

subroutine lsm_smooth_core(phi1_, phi0_, dtau, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi1_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi0_
	real										:: dtau
	real										:: dpx_p1, dpy_p1, dpz_p1
	real										:: dpx_n1, dpy_n1, dpz_n1
	real										:: dpx_p_max2, dpy_p_max2, dpz_p_max2
	real										:: dpx_n_max2, dpy_n_max2, dpz_n_max2
	real										:: dpx_p_min2, dpy_p_min2, dpz_p_min2
	real										:: dpx_n_min2, dpy_n_min2, dpz_n_min2
	real										:: dpx_p2, dpy_p2, dpz_p2
	real										:: dpx_n2, dpy_n2, dpz_n2
	real										:: dp_p, dp_n
	real										:: dpl
	real										:: phis
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(dpx_p1, dpy_p1, dpz_p1) &
!$omp					,private(dpx_n1, dpy_n1, dpz_n1) &
!$omp					,private(dpx_p_max2, dpy_p_max2, dpz_p_max2) &
!$omp					,private(dpx_n_max2, dpy_n_max2, dpz_n_max2) &
!$omp					,private(dpx_p_min2, dpy_p_min2, dpz_p_min2) &
!$omp					,private(dpx_n_min2, dpy_n_min2, dpz_n_min2) &
!$omp					,private(dpx_p2, dpy_p2, dpz_p2) &
!$omp					,private(dpx_n2, dpy_n2, dpz_n2) &
!$omp					,private(dp_p, dp_n) &
!$omp					,private(dpl) &
!$omp					,private(phis)
!$omp do schedule(dynamic, 1)
	do k=1, kx
	do j=1, jx
	do i=1, ix
		dpx_p1 = phi0_(i+1, j, k) - phi0_(i  , j, k)
		dpx_n1 = phi0_(i  , j, k) - phi0_(i-1, j, k)
		dpy_p1 = phi0_(i, j+1, k) - phi0_(i, j  , k)
		dpy_n1 = phi0_(i, j  , k) - phi0_(i, j-1, k)
		dpz_p1 = phi0_(i, j, k+1) - phi0_(i, j, k  )
		dpz_n1 = phi0_(i, j, k  ) - phi0_(i, j, k-1)

		dpx_n_max2 = max(dpx_n1, 0.0)*max(dpx_n1, 0.0)
		dpx_n_min2 = min(dpx_n1, 0.0)*min(dpx_n1, 0.0)
		dpx_p_max2 = max(dpx_p1, 0.0)*max(dpx_p1, 0.0)
		dpx_p_min2 = min(dpx_p1, 0.0)*min(dpx_p1, 0.0)
		dpy_n_max2 = max(dpy_n1, 0.0)*max(dpy_n1, 0.0)
		dpy_n_min2 = min(dpy_n1, 0.0)*min(dpy_n1, 0.0)
		dpy_p_max2 = max(dpy_p1, 0.0)*max(dpy_p1, 0.0)
		dpy_p_min2 = min(dpy_p1, 0.0)*min(dpy_p1, 0.0)
		dpz_n_max2 = max(dpz_n1, 0.0)*max(dpz_n1, 0.0)
		dpz_n_min2 = min(dpz_n1, 0.0)*min(dpz_n1, 0.0)
		dpz_p_max2 = max(dpz_p1, 0.0)*max(dpz_p1, 0.0)
		dpz_p_min2 = min(dpz_p1, 0.0)*min(dpz_p1, 0.0)

		dpx_p2 = max(dpx_n_max2, dpx_p_min2)
		dpx_n2 = max(dpx_n_min2, dpx_p_max2)
		dpy_p2 = max(dpy_n_max2, dpy_p_min2)
		dpy_n2 = max(dpy_n_min2, dpy_p_max2)
		dpz_p2 = max(dpz_n_max2, dpz_p_min2)
		dpz_n2 = max(dpz_n_min2, dpz_p_max2)

		dp_p = sqrt(dpx_p2 + dpy_p2 + dpz_p2)
		dp_n = sqrt(dpx_n2 + dpy_n2 + dpz_n2)

		if(phi0_(i, j, k) > 0.0) then
			dpl = dp_p
			phis = 1.0
		else if(phi0_(i, j, k) < 0.0) then
			dpl = dp_n
			phis = -1.0
		else
			dpl = 1.0
			phis = 0.0
		endif

		phi1_(i, j, k) = phi0_(i, j, k) - phis*(dpl - 1.0)*dtau
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine lsm_smooth_core

subroutine lsm_extend_core(data1_, data0_, dtau, phi, nx, ny, nz, psi, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: data1_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: data0_
	real										:: dtau
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: nx, ny, nz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: psi
	real										:: s
	real										:: wx, wy, wz
	real										:: ddx_p1, ddy_p1, ddz_p1
	real										:: ddx_n1, ddy_n1, ddz_n1
	real										:: ddx_p2, ddy_p2, ddz_p2
	real										:: ddx_n2, ddy_n2, ddz_n2
	real										:: ddxx_n, ddxx_c, ddxx_p
	real										:: ddyy_n, ddyy_c, ddyy_p
	real										:: ddzz_n, ddzz_c, ddzz_p
	real										:: ddx_p, ddy_p, ddz_p
	real										:: ddx_n, ddy_n, ddz_n
	real										:: lsm_getsign
	real										:: lsm_getminmod
	real										:: lsm_getweno3
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(s) &
!$omp					,private(wx, wy, wz) &
!$omp					,private(ddx_p1, ddy_p1, ddz_p1) &
!$omp					,private(ddx_n1, ddy_n1, ddz_n1) &
!$omp					,private(ddx_p2, ddy_p2, ddz_p2) &
!$omp					,private(ddx_n2, ddy_n2, ddz_n2) &
!$omp					,private(ddxx_n, ddxx_c, ddxx_p) &
!$omp					,private(ddyy_n, ddyy_c, ddyy_p) &
!$omp					,private(ddzz_n, ddzz_c, ddzz_p) &
!$omp					,private(ddx_p, ddy_p, ddz_p) &
!$omp					,private(ddx_n, ddy_n, ddz_n) 
!$omp do schedule(dynamic, 1)
	do k=1, kx
	do j=1, jx
	do i=1, ix
		if( psi(i, j, k) >= 0 ) then
			data1_(i, j, k) = data0_(i, j, k)
			cycle
		endif

		if( abs(phi(i, j, k)) > 4.0 ) then
!			cycle
		endif

		s  = lsm_getsign( phi(i, j, k) )
		wx = nx(i, j, k)*s
		wy = ny(i, j, k)*s
		wz = nz(i, j, k)*s

		ddx_p1 = data0_(i+1, j, k) - data0_(i  , j, k)
		ddx_n1 = data0_(i  , j, k) - data0_(i-1, j, k)
		ddy_p1 = data0_(i, j+1, k) - data0_(i, j  , k)
		ddy_n1 = data0_(i, j  , k) - data0_(i, j-1, k)
		ddz_p1 = data0_(i, j, k+1) - data0_(i, j, k  )
		ddz_n1 = data0_(i, j, k  ) - data0_(i, j, k-1)
		ddxx_p = data0_(i+2, j, k) + data0_(i  , j, k) - 2.0*data0_(i+1, j, k)
		ddxx_c = data0_(i+1, j, k) + data0_(i-1, j, k) - 2.0*data0_(i  , j, k)
		ddxx_n = data0_(i  , j, k) + data0_(i-2, j, k) - 2.0*data0_(i-1, j, k)
		ddyy_p = data0_(i, j+2, k) + data0_(i, j  , k) - 2.0*data0_(i, j+1, k)
		ddyy_c = data0_(i, j+1, k) + data0_(i, j-1, k) - 2.0*data0_(i, j  , k)
		ddyy_n = data0_(i, j  , k) + data0_(i, j-2, k) - 2.0*data0_(i, j-1, k)
		ddzz_p = data0_(i, j, k+2) + data0_(i, j, k  ) - 2.0*data0_(i, j, k+1)
		ddzz_c = data0_(i, j, k+1) + data0_(i, j, k-1) - 2.0*data0_(i, j, k  )
		ddzz_n = data0_(i, j, k  ) + data0_(i, j, k-2) - 2.0*data0_(i, j, k-1)
		ddx_p2 = ddx_p1 - 0.5*lsm_getminmod(ddxx_c, ddxx_p)
		ddx_n2 = ddx_n1 + 0.5*lsm_getminmod(ddxx_c, ddxx_n)
		ddy_p2 = ddy_p1 - 0.5*lsm_getminmod(ddyy_c, ddyy_p)
		ddy_n2 = ddy_n1 + 0.5*lsm_getminmod(ddyy_c, ddyy_n)
		ddz_p2 = ddz_p1 - 0.5*lsm_getminmod(ddzz_c, ddzz_p)
		ddz_n2 = ddz_n1 + 0.5*lsm_getminmod(ddzz_c, ddzz_n)
		ddx_p = ddx_p1
		ddx_n = ddx_n1
		ddy_p = ddy_p1
		ddy_n = ddy_n1
		ddz_p = ddz_p1
		ddz_n = ddz_n1

		if( psi(i+1, j, k) >= 0 ) then
			ddx_p = 0.0
		endif
		if( psi(i-1, j, k) >= 0 ) then
			ddx_n = 0.0
		endif
		if( psi(i, j+1, k) >= 0 ) then
			ddy_p = 0.0
		endif
		if( psi(i, j-1, k) >= 0 ) then
			ddy_n = 0.0
		endif
		if( psi(i, j, k+1) >= 0 ) then
			ddz_p = 0.0
		endif
		if( psi(i, j, k-1) >= 0 ) then
			ddz_n = 0.0
		endif

		data1_(i, j, k) = data0_(i, j, k) -( &
																			+ (max(wx, 0.0)*ddx_n + min(wx, 0.0)*ddx_p) &
																			+ (max(wy, 0.0)*ddy_n + min(wy, 0.0)*ddy_p) &
																			+ (max(wz, 0.0)*ddz_n + min(wz, 0.0)*ddz_p) )*dtau
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine lsm_extend_core

subroutine lsm_extend_s_core(data1_, data0_, dtau, psi, psinx, psiny, psinz, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: data1_
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: data0_
	real										:: dtau
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: psi
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: psinx, psiny, psinz
	real										:: s
	real										:: wx, wy, wz
	real										:: ddx_p1, ddy_p1, ddz_p1
	real										:: ddx_n1, ddy_n1, ddz_n1
	real										:: ddx_p2, ddy_p2, ddz_p2
	real										:: ddx_n2, ddy_n2, ddz_n2
	real										:: ddxx_n, ddxx_c, ddxx_p
	real										:: ddyy_n, ddyy_c, ddyy_p
	real										:: ddzz_n, ddzz_c, ddzz_p
	real										:: ddx_p, ddy_p, ddz_p
	real										:: ddx_n, ddy_n, ddz_n
	real										:: lsm_getsign
	real										:: lsm_getminmod
	real										:: lsm_getweno3
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(s) &
!$omp					,private(wx, wy, wz) &
!$omp					,private(ddx_p1, ddy_p1, ddz_p1) &
!$omp					,private(ddx_n1, ddy_n1, ddz_n1) &
!$omp					,private(ddx_p2, ddy_p2, ddz_p2) &
!$omp					,private(ddx_n2, ddy_n2, ddz_n2) &
!$omp					,private(ddxx_n, ddxx_c, ddxx_p) &
!$omp					,private(ddyy_n, ddyy_c, ddyy_p) &
!$omp					,private(ddzz_n, ddzz_c, ddzz_p) &
!$omp					,private(ddx_p, ddy_p, ddz_p) &
!$omp					,private(ddx_n, ddy_n, ddz_n) 
!$omp do schedule(dynamic, 1)
	do k=1, kx
	do j=1, jx
	do i=1, ix
		if( psi(i, j, k) < 0 ) then
			data1_(i, j, k) = data0_(i, j, k)
			cycle
		endif

		if( abs(psi(i, j, k)) > 8.0 ) then
!			cycle
		endif

		wx = psinx(i, j, k)
		wy = psiny(i, j, k)
		wz = psinz(i, j, k)

		ddx_p1 = data0_(i+1, j, k) - data0_(i  , j, k)
		ddx_n1 = data0_(i  , j, k) - data0_(i-1, j, k)
		ddy_p1 = data0_(i, j+1, k) - data0_(i, j  , k)
		ddy_n1 = data0_(i, j  , k) - data0_(i, j-1, k)
		ddz_p1 = data0_(i, j, k+1) - data0_(i, j, k  )
		ddz_n1 = data0_(i, j, k  ) - data0_(i, j, k-1)
		ddxx_p = data0_(i+2, j, k) + data0_(i  , j, k) - 2.0*data0_(i+1, j, k)
		ddxx_c = data0_(i+1, j, k) + data0_(i-1, j, k) - 2.0*data0_(i  , j, k)
		ddxx_n = data0_(i  , j, k) + data0_(i-2, j, k) - 2.0*data0_(i-1, j, k)
		ddyy_p = data0_(i, j+2, k) + data0_(i, j  , k) - 2.0*data0_(i, j+1, k)
		ddyy_c = data0_(i, j+1, k) + data0_(i, j-1, k) - 2.0*data0_(i, j  , k)
		ddyy_n = data0_(i, j  , k) + data0_(i, j-2, k) - 2.0*data0_(i, j-1, k)
		ddzz_p = data0_(i, j, k+2) + data0_(i, j, k  ) - 2.0*data0_(i, j, k+1)
		ddzz_c = data0_(i, j, k+1) + data0_(i, j, k-1) - 2.0*data0_(i, j, k  )
		ddzz_n = data0_(i, j, k  ) + data0_(i, j, k-2) - 2.0*data0_(i, j, k-1)
		ddx_p2 = ddx_p1 - 0.5*lsm_getminmod(ddxx_c, ddxx_p)
		ddx_n2 = ddx_n1 + 0.5*lsm_getminmod(ddxx_c, ddxx_n)
		ddy_p2 = ddy_p1 - 0.5*lsm_getminmod(ddyy_c, ddyy_p)
		ddy_n2 = ddy_n1 + 0.5*lsm_getminmod(ddyy_c, ddyy_n)
		ddz_p2 = ddz_p1 - 0.5*lsm_getminmod(ddzz_c, ddzz_p)
		ddz_n2 = ddz_n1 + 0.5*lsm_getminmod(ddzz_c, ddzz_n)
		ddx_p = ddx_p1
		ddx_n = ddx_n1
		ddy_p = ddy_p1
		ddy_n = ddy_n1
		ddz_p = ddz_p1
		ddz_n = ddz_n1

		data1_(i, j, k) = data0_(i, j, k) -( &
																			+ (max(wx, 0.0)*ddx_n + min(wx, 0.0)*ddx_p) &
																			+ (max(wy, 0.0)*ddy_n + min(wy, 0.0)*ddy_p) &
																			+ (max(wz, 0.0)*ddz_n + min(wz, 0.0)*ddz_p) )*dtau
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine lsm_extend_s_core

