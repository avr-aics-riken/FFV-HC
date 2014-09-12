!>  @file  bsdf.f90
!!  @brief Basic subprograms for initializing the signed distance functions for Cartesian grid data structure
!<

subroutine bsdf_liquid(phi, psi, xc, yc, zc, rc, lc, wc, hc, sz, g, b, gsz)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	integer									:: i1, j1, k1
	integer, dimension(3)		:: b
	integer									:: gix, gjx, gkx
	integer, dimension(3)		:: gsz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi, psi
	real										:: xc, yc, zc, rc, lc, wc, hc
	real										:: lx, ly, lz
	real										:: gx, gy, gz
	real										:: x0, y0, z0
	real										:: r0
	real										:: rx, ry, rz, r2, rl
	real										:: phi0, phi1, phi2, phi3
	real										:: psi0, psi1, psi2, psi3
	gix = gsz(1)
	gjx = gsz(2)
	gkx = gsz(3)
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i1 = b(1)*sz(1)
	j1 = b(2)*sz(2)
	k1 = b(3)*sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(lx, ly, lz) &
!$omp					,private(gx, gy, gz) &
!$omp					,private(x0, y0, z0) &
!$omp					,private(r0) &
!$omp					,private(rx, ry, rz, r2, rl) &
!$omp					,private(phi0, phi1, phi2, phi3) &
!$omp					,private(psi0, psi1, psi2, psi3) 
!$omp do schedule(static, 1)
	do k=1-g, kx+g
	do j=1-g, jx+g
	do i=1-g, ix+g
		lx = i - 0.5
		ly = j - 0.5
		lz = k - 0.5
		gx = i + i1 - 0.5
		gy = j + j1 - 0.5
		gz = k + k1 - 0.5

		rx = gx - xc
		ry = gy - yc
		rz = gz - zc
		r2 = rx*rx + ry*ry + rz*rz
		rl = sqrt(r2)
		psi(i, j, k) = -(rl + 100.0)
		phi(i, j, k) =   rl + 100.0
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bsdf_liquid

subroutine bsdf_pipe(phi, psi, xc, yc, zc, rc, lc, wc, hc, sz, g, b, gsz)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	integer									:: i1, j1, k1
	integer, dimension(3)		:: b
	integer									:: gix, gjx, gkx
	integer, dimension(3)		:: gsz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi, psi
	real										:: xc, yc, zc, rc, lc, wc, hc
	real										:: lx, ly, lz
	real										:: gx, gy, gz
	real										:: x0, y0, z0
	real										:: r0
	real										:: rx, ry, rz, r2, rl
	real										:: phi0, phi1, phi2, phi3
	real										:: psi0, psi1, psi2, psi3
	gix = gsz(1)
	gjx = gsz(2)
	gkx = gsz(3)
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i1 = b(1)*sz(1)
	j1 = b(2)*sz(2)
	k1 = b(3)*sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(lx, ly, lz) &
!$omp					,private(gx, gy, gz) &
!$omp					,private(x0, y0, z0) &
!$omp					,private(r0) &
!$omp					,private(rx, ry, rz, r2, rl) &
!$omp					,private(phi0, phi1, phi2, phi3) &
!$omp					,private(psi0, psi1, psi2, psi3) 
!$omp do schedule(static, 1)
	do k=1-g, kx+g
	do j=1-g, jx+g
	do i=1-g, ix+g
		lx = i - 0.5
		ly = j - 0.5
		lz = k - 0.5
		gx = i + i1 - 0.5
		gy = j + j1 - 0.5
		gz = k + k1 - 0.5

		rx = gx - xc
		ry = gy - yc
		rz = gz - zc
		r2 = rx*rx + ry*ry + rz*rz
		r2 = rx*rx + ry*ry
		r2 = ry*ry + rz*rz
		rl = sqrt(r2)
		psi(i, j, k) = rl - rc

		phi0 =-(gx - (lc - 0.5*wc))
		phi1 = (gx - (lc + 0.5*wc))
		phi(i, j, k) = phi0
		if( abs(phi1) < abs(phi0) ) then
			phi(i, j, k) = phi1
		endif
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bsdf_pipe

subroutine bsdf_pipe2(phi, psi, xc, yc, zc, rc, lc, wc, hc, sz, g, b, gsz)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	integer									:: i1, j1, k1
	integer, dimension(3)		:: b
	integer									:: gix, gjx, gkx
	integer, dimension(3)		:: gsz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi, psi
	real										:: xc, yc, zc, rc, lc, wc, hc
	real										:: lx, ly, lz
	real										:: gx, gy, gz
	real										:: x0, y0, z0
	real										:: r0
	real										:: rx, ry, rz, r2, rl
	real										:: phi0, phi1, phi2, phi3
	real										:: psi0, psi1, psi2, psi3
	gix = gsz(1)
	gjx = gsz(2)
	gkx = gsz(3)
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i1 = b(1)*sz(1)
	j1 = b(2)*sz(2)
	k1 = b(3)*sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(lx, ly, lz) &
!$omp					,private(gx, gy, gz) &
!$omp					,private(x0, y0, z0) &
!$omp					,private(r0) &
!$omp					,private(rx, ry, rz, r2, rl) &
!$omp					,private(phi0, phi1, phi2, phi3) &
!$omp					,private(psi0, psi1, psi2, psi3) 
!$omp do schedule(static, 1)
	do k=1-g, kx+g
	do j=1-g, jx+g
	do i=1-g, ix+g
		lx = i - 0.5
		ly = j - 0.5
		lz = k - 0.5
		gx = i + i1 - 0.5
		gy = j + j1 - 0.5
		gz = k + k1 - 0.5

		rx = gx - xc
		ry = gy - yc
		rz = gz - zc
		r2 = rx*rx + ry*ry + rz*rz
		r2 = rx*rx + ry*ry
		r2 = ry*ry + rz*rz
		rl = sqrt(r2)
		psi(i, j, k) = rl - rc

!		phi(i, j, k) = -(gx - wc)

		rx = gx - xc
		ry = gy - yc
		rz = gz - zc
		r2 = rx*rx + ry*ry + rz*rz
		rl = sqrt(r2)
		phi0 = rl - rc
		phi1 = -(gx - wc)

		phi(i, j, k) = phi1
		if( gx < wc + rc ) then
			phi(i, j, k) = phi0
		end if

!		phi(i, j, k) = phi0
!		if( abs(phi1) < abs(phi0) ) then
!			phi(i, j, k) = phi1
!		endif
!		phi(i, j, k) = phi1
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bsdf_pipe2

subroutine bsdf_bubble2d(phi, psi, xc, yc, zc, rc, lc, wc, hc, sz, g, b, gsz)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	integer									:: i1, j1, k1
	integer, dimension(3)		:: b
	integer									:: gix, gjx, gkx
	integer, dimension(3)		:: gsz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi, psi
	real										:: xc, yc, zc, rc, lc, wc, hc
	real										:: lx, ly, lz
	real										:: gx, gy, gz
	real										:: x0, y0, z0
	real										:: r0
	real										:: rx, ry, rz, r2, rl
	real										:: phi0, phi1, phi2, phi3
	real										:: psi0, psi1, psi2, psi3
	gix = gsz(1)
	gjx = gsz(2)
	gkx = gsz(3)
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i1 = b(1)*sz(1)
	j1 = b(2)*sz(2)
	k1 = b(3)*sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(lx, ly, lz) &
!$omp					,private(gx, gy, gz) &
!$omp					,private(x0, y0, z0) &
!$omp					,private(r0) &
!$omp					,private(rx, ry, rz, r2, rl) &
!$omp					,private(phi0, phi1, phi2, phi3) &
!$omp					,private(psi0, psi1, psi2, psi3) 
!$omp do schedule(static, 1)
	do k=1-g, kx+g
	do j=1-g, jx+g
	do i=1-g, ix+g
		lx = i - 0.5
		ly = j - 0.5
		lz = k - 0.5
		gx = i + i1 - 0.5
		gy = j + j1 - 0.5
		gz = k + k1 - 0.5

		rx = gx - xc
		ry = gy - yc
		rz = gz - zc
		r2 = rx*rx + ry*ry + rz*rz
		r2 = rx*rx + ry*ry
		r2 = rx*rx + rz*rz
		rl = sqrt(r2)

		phi(i, j, k) =  (rl - rc)
		psi(i, j, k) = -(rl + 100.0)
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bsdf_bubble2d

subroutine bsdf_drop(phi, psi, xc, yc, zc, rc, lc, wc, hc, sz, g, b, gsz)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	integer									:: i1, j1, k1
	integer, dimension(3)		:: b
	integer									:: gix, gjx, gkx
	integer, dimension(3)		:: gsz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi, psi
	real										:: xc, yc, zc, rc, lc, wc, hc
	real										:: lx, ly, lz
	real										:: gx, gy, gz
	real										:: x0, y0, z0
	real										:: r0
	real										:: rx, ry, rz, r2, rl
	real										:: phi0, phi1, phi2, phi3
	real										:: psi0, psi1, psi2, psi3
	gix = gsz(1)
	gjx = gsz(2)
	gkx = gsz(3)
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i1 = b(1)*sz(1)
	j1 = b(2)*sz(2)
	k1 = b(3)*sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(lx, ly, lz) &
!$omp					,private(gx, gy, gz) &
!$omp					,private(x0, y0, z0) &
!$omp					,private(r0) &
!$omp					,private(rx, ry, rz, r2, rl) &
!$omp					,private(phi0, phi1, phi2, phi3) &
!$omp					,private(psi0, psi1, psi2, psi3) 
!$omp do schedule(static, 1)
	do k=1-g, kx+g
	do j=1-g, jx+g
	do i=1-g, ix+g
		lx = i - 0.5
		ly = j - 0.5
		lz = k - 0.5
		gx = i + i1 - 0.5
		gy = j + j1 - 0.5
		gz = k + k1 - 0.5

		rx = gx - xc
		ry = gy - yc
		rz = gz - zc
		r2 = rx*rx + ry*ry + rz*rz
		rl = sqrt(r2)

		phi(i, j, k) = -(rl - rc)
		psi(i, j, k) = -(rl + 100.0)
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bsdf_drop

subroutine bsdf_bubble(phi, psi, xc, yc, zc, rc, lc, wc, hc, sz, g, b, gsz)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	integer									:: i1, j1, k1
	integer, dimension(3)		:: b
	integer									:: gix, gjx, gkx
	integer, dimension(3)		:: gsz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi, psi
	real										:: xc, yc, zc, rc, lc, wc, hc
	real										:: lx, ly, lz
	real										:: gx, gy, gz
	real										:: x0, y0, z0
	real										:: r0
	real										:: rx, ry, rz, r2, rl
	real										:: phi0, phi1, phi2, phi3
	real										:: psi0, psi1, psi2, psi3
	gix = gsz(1)
	gjx = gsz(2)
	gkx = gsz(3)
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i1 = b(1)*sz(1)
	j1 = b(2)*sz(2)
	k1 = b(3)*sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(lx, ly, lz) &
!$omp					,private(gx, gy, gz) &
!$omp					,private(x0, y0, z0) &
!$omp					,private(r0) &
!$omp					,private(rx, ry, rz, r2, rl) &
!$omp					,private(phi0, phi1, phi2, phi3) &
!$omp					,private(psi0, psi1, psi2, psi3) 
!$omp do schedule(static, 1)
	do k=1-g, kx+g
	do j=1-g, jx+g
	do i=1-g, ix+g
		lx = i - 0.5
		ly = j - 0.5
		lz = k - 0.5
		gx = i + i1 - 0.5
		gy = j + j1 - 0.5
		gz = k + k1 - 0.5

		rx = gx - xc
		ry = gy - yc
		rz = gz - zc
		r2 = rx*rx + ry*ry + rz*rz
		rl = sqrt(r2)

		phi(i, j, k) =  (rl - rc)
		psi(i, j, k) = -(rl + 100.0)
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bsdf_bubble

subroutine bsdf_beaker(phi, psi, xc, yc, zc, rc, lc, wc, hc, sz, g, b, gsz)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	integer									:: i1, j1, k1
	integer, dimension(3)		:: b
	integer									:: gix, gjx, gkx
	integer, dimension(3)		:: gsz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi, psi
	real										:: xc, yc, zc, rc, lc, wc, hc
	real										:: lx, ly, lz
	real										:: gx, gy, gz
	real										:: x0, y0, z0
	real										:: r0
	real										:: rx, ry, rz, r2, rl
	real										:: phi0, phi1, phi2, phi3
	real										:: psi0, psi1, psi2, psi3
	gix = gsz(1)
	gjx = gsz(2)
	gkx = gsz(3)
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i1 = b(1)*sz(1)
	j1 = b(2)*sz(2)
	k1 = b(3)*sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(lx, ly, lz) &
!$omp					,private(gx, gy, gz) &
!$omp					,private(x0, y0, z0) &
!$omp					,private(r0) &
!$omp					,private(rx, ry, rz, r2, rl) &
!$omp					,private(phi0, phi1, phi2, phi3) &
!$omp					,private(psi0, psi1, psi2, psi3) 
!$omp do schedule(static, 1)
	do k=1-g, kx+g
	do j=1-g, jx+g
	do i=1-g, ix+g
		lx = i - 0.5
		ly = j - 0.5
		lz = k - 0.5
		gx = i + i1 - 0.5
		gy = j + j1 - 0.5
		gz = k + k1 - 0.5

		rx = gx - xc
		ry = gy - yc
		rz = gz - zc
		r2 = rx*rx + ry*ry + rz*rz
		r2 = ry*ry + rz*rz
		r2 = rx*rx + ry*ry
		rl = sqrt(r2)
		psi(i, j, k) = -rl

		phi0 =-(gz - (lc - 0.5*wc))
		phi1 = (gz - (lc + 0.5*wc))
		phi(i, j, k) = phi0
		if( abs(phi1) < abs(phi0) ) then
			phi(i, j, k) = phi1
		endif
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bsdf_beaker

subroutine bsdf_ducky(phi, psi, xc, yc, zc, rc, lc, wc, hc, sz, g, b, gsz)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	integer									:: i1, j1, k1
	integer, dimension(3)		:: b
	integer									:: gix, gjx, gkx
	integer, dimension(3)		:: gsz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi, psi
	real										:: xc, yc, zc, rc, lc, wc, hc
	real										:: lx, ly, lz
	real										:: gx, gy, gz
	real										:: x0, y0, z0
	real										:: r0
	real										:: rx, ry, rz, r2, rl
	real										:: phi0, phi1, phi2, phi3
	real										:: psi0, psi1, psi2, psi3
	gix = gsz(1)
	gjx = gsz(2)
	gkx = gsz(3)
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i1 = b(1)*sz(1)
	j1 = b(2)*sz(2)
	k1 = b(3)*sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(lx, ly, lz) &
!$omp					,private(gx, gy, gz) &
!$omp					,private(x0, y0, z0) &
!$omp					,private(r0) &
!$omp					,private(rx, ry, rz, r2, rl) &
!$omp					,private(phi0, phi1, phi2, phi3) &
!$omp					,private(psi0, psi1, psi2, psi3) 
!$omp do schedule(static, 1)
	do k=1-g, kx+g
	do j=1-g, jx+g
	do i=1-g, ix+g
		lx = i - 0.5
		ly = j - 0.5
		lz = k - 0.5
		gx = i + i1 - 0.5
		gy = j + j1 - 0.5
		gz = k + k1 - 0.5

		rx = gx - xc
		ry = gy - yc
		rz = gz - zc
		r2 = rx*rx + ry*ry + rz*rz
		r2 = ry*ry + rz*rz
		r2 = rx*rx + ry*ry
		rl = sqrt(r2)
		psi(i, j, k) = -rl

		phi0 = (gz - (lc - 0.5*wc))
		phi1 =-(gz - (lc + 0.5*wc))
		phi(i, j, k) = phi0
		if( abs(phi1) < abs(phi0) ) then
			phi(i, j, k) = phi1
		endif
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bsdf_ducky

subroutine bsdf_dambreak2d(phi, psi, xc, yc, zc, rc, lc, wc, hc, sz, g, b, gsz)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	integer									:: i1, j1, k1
	integer, dimension(3)		:: b
	integer									:: gix, gjx, gkx
	integer, dimension(3)		:: gsz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: phi, psi
	real										:: xc, yc, zc, rc, lc, wc, hc
	real										:: lx, ly, lz
	real										:: gx, gy, gz
	real										:: x0, y0, z0
	real										:: r0
	real										:: rx, ry, rz, r2, rl
	real										:: phi0, phi1, phi2, phi3
	real										:: psi0, psi1, psi2, psi3
	gix = gsz(1)
	gjx = gsz(2)
	gkx = gsz(3)
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i1 = b(1)*sz(1)
	j1 = b(2)*sz(2)
	k1 = b(3)*sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(lx, ly, lz) &
!$omp					,private(gx, gy, gz) &
!$omp					,private(x0, y0, z0) &
!$omp					,private(r0) &
!$omp					,private(rx, ry, rz, r2, rl) &
!$omp					,private(phi0, phi1, phi2, phi3) &
!$omp					,private(psi0, psi1, psi2, psi3) 
!$omp do schedule(static, 1)
	do k=1-g, kx+g
	do j=1-g, jx+g
	do i=1-g, ix+g
		lx = i - 0.5
		ly = j - 0.5
		lz = k - 0.5
		gx = i + i1 - 0.5
		gy = j + j1 - 0.5
		gz = k + k1 - 0.5

		rx = gx - xc
		ry = gy - yc
		rz = gz - zc
		r2 = rx*rx + rz*rz
		rl = sqrt(r2)

		phi0 = -(gx - wc)
		phi1 = -(gz - hc)

		if( phi0 >= 0 .and. phi1 >= 0 ) then
			phi2 = phi0
			if( abs(phi1) < abs(phi0) ) then
				phi2 = phi1
			endif
		end if
		if( phi0 >= 0 .and. phi1 < 0 ) then
			phi2 = phi1
		end if
		if( phi0 < 0 .and. phi1 >= 0 ) then
			phi2 = phi0
		end if
		if( phi0 < 0 .and. phi1 < 0 ) then
			rx = gx - wc
			rz = gz - hc
			r2 = rx*rx + rz*rz
			rl = sqrt(r2)
			phi2 = -r2
		end if

		phi(i, j, k) = phi2
		psi(i, j, k) = -(rl + 100.0)
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bsdf_dambreak2d

