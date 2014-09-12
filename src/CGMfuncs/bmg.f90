!>  @file  bmg.f90
!!  @brief Basic subprograms for the multigrid method for Cartesian grid data structure
!<

function mg_geth(phi)
	implicit none
	real					:: mg_geth
	real					:: phi
	mg_geth = 1.0
	if( phi < 0 ) then
		mg_geth = 0.0
	endif
	return
end function mg_geth

subroutine mg_restriction_1(xc, x, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)/2+g, 1-g:sz(2)/2+g, 1-g:sz(3)/2+g)	:: xc
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)				:: x
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j, k)
!$omp do
	do k=1,kx/2
	do j=1,jx/2
	do i=1,ix/2
		xc(i, j, k) = 0.125*( x(2*i    , 2*j    , 2*k    ) &
												+ x(2*i - 1, 2*j    , 2*k    ) &
												+ x(2*i    , 2*j - 1, 2*k    ) &
												+ x(2*i    , 2*j    , 2*k - 1) &
												+ x(2*i - 1, 2*j - 1, 2*k    ) &
												+ x(2*i    , 2*j - 1, 2*k - 1) &
												+ x(2*i - 1, 2*j    , 2*k - 1) &
												+ x(2*i - 1, 2*j - 1, 2*k - 1) )
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine mg_restriction_1

subroutine mg_prolongation_1(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)/2+g, 1-g:sz(2)/2+g, 1-g:sz(3)/2+g)	:: xc
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)				:: x
	real										:: a
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	a = 4.0
!$omp parallel private(i, j, k)
!$omp do
	do k=1,kx/2
	do j=1,jx/2
	do i=1,ix/2
		x(2*i    , 2*j    , 2*k    ) = x(2*i    , 2*j    , 2*k    ) + xc(i, j, k)*a
		x(2*i - 1, 2*j    , 2*k    ) = x(2*i - 1, 2*j    , 2*k    ) + xc(i, j, k)*a
		x(2*i    , 2*j - 1, 2*k    ) = x(2*i    , 2*j - 1, 2*k    ) + xc(i, j, k)*a
		x(2*i    , 2*j    , 2*k - 1) = x(2*i    , 2*j    , 2*k - 1) + xc(i, j, k)*a
		x(2*i - 1, 2*j - 1, 2*k    ) = x(2*i - 1, 2*j - 1, 2*k    ) + xc(i, j, k)*a
		x(2*i    , 2*j - 1, 2*k - 1) = x(2*i    , 2*j - 1, 2*k - 1) + xc(i, j, k)*a
		x(2*i - 1, 2*j    , 2*k - 1) = x(2*i - 1, 2*j    , 2*k - 1) + xc(i, j, k)*a
		x(2*i - 1, 2*j - 1, 2*k - 1) = x(2*i - 1, 2*j - 1, 2*k - 1) + xc(i, j, k)*a
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine mg_prolongation_1

subroutine mg_restriction_2(bc, r, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)/2+g, 1-g:sz(2)/2+g, 1-g:sz(3)/2+g)	:: bc
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)				:: r
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j, k)
!$omp do 
	do k=1,kx/2
	do j=1,jx/2
	do i=1,ix/2
		bc(i, j, k) = 0.125*(1.0/64.0)*( &
										27.0*( &
													r(2*i    , 2*j    , 2*k    ) &
												+ r(2*i - 1, 2*j    , 2*k    ) &
												+ r(2*i    , 2*j - 1, 2*k    ) &
												+ r(2*i    , 2*j    , 2*k - 1) &
												+ r(2*i - 1, 2*j - 1, 2*k    ) &
												+ r(2*i    , 2*j - 1, 2*k - 1) &
												+ r(2*i - 1, 2*j    , 2*k - 1) &
												+ r(2*i - 1, 2*j - 1, 2*k - 1) &
										) &
									+  9.0*( &
													r(2*i + 1, 2*j    , 2*k    ) + r(2*i    , 2*j + 1, 2*k    ) + r(2*i    , 2*j    , 2*k + 1) &
												+ r(2*i - 2, 2*j    , 2*k    ) + r(2*i - 1, 2*j + 1, 2*k    ) + r(2*i - 1, 2*j    , 2*k + 1) &
												+ r(2*i + 1, 2*j - 1, 2*k    ) + r(2*i    , 2*j - 2, 2*k    ) + r(2*i    , 2*j - 1, 2*k + 1) &
												+ r(2*i + 1, 2*j    , 2*k - 1) + r(2*i    , 2*j + 1, 2*k - 1) + r(2*i    , 2*j    , 2*k - 2) &
												+ r(2*i - 2, 2*j - 1, 2*k    ) + r(2*i - 1, 2*j - 2, 2*k    ) + r(2*i - 1, 2*j - 1, 2*k + 1) &
												+ r(2*i + 1, 2*j - 1, 2*k - 1) + r(2*i    , 2*j - 2, 2*k - 1) + r(2*i    , 2*j - 1, 2*k - 2) &
												+ r(2*i - 2, 2*j    , 2*k - 1) + r(2*i - 1, 2*j + 1, 2*k - 1) + r(2*i - 1, 2*j    , 2*k - 2) &
												+ r(2*i - 2, 2*j - 1, 2*k - 1) + r(2*i - 1, 2*j - 2, 2*k - 1) + r(2*i - 1, 2*j - 1, 2*k - 2) &
										) &
									+  3.0*( &
													r(2*i + 1, 2*j + 1, 2*k    ) + r(2*i    , 2*j + 1, 2*k + 1) + r(2*i + 1, 2*j    , 2*k + 1) &
												+ r(2*i - 2, 2*j + 1, 2*k    ) + r(2*i - 1, 2*j + 1, 2*k + 1) + r(2*i - 2, 2*j    , 2*k + 1) &
												+ r(2*i + 1, 2*j - 2, 2*k    ) + r(2*i    , 2*j - 2, 2*k + 1) + r(2*i + 1, 2*j - 1, 2*k + 1) &
												+ r(2*i + 1, 2*j + 1, 2*k - 1) + r(2*i    , 2*j + 1, 2*k - 2) + r(2*i + 1, 2*j    , 2*k - 2) &
												+ r(2*i - 2, 2*j - 2, 2*k    ) + r(2*i - 1, 2*j - 2, 2*k + 1) + r(2*i - 2, 2*j - 1, 2*k + 1) &
												+ r(2*i + 1, 2*j - 2, 2*k - 1) + r(2*i    , 2*j - 2, 2*k - 2) + r(2*i + 1, 2*j - 1, 2*k - 2) &
												+ r(2*i - 2, 2*j + 1, 2*k - 1) + r(2*i - 1, 2*j + 1, 2*k - 2) + r(2*i - 2, 2*j    , 2*k - 2) &
												+ r(2*i - 2, 2*j - 2, 2*k - 1) + r(2*i - 1, 2*j - 2, 2*k - 2) + r(2*i - 2, 2*j - 1, 2*k - 2) &
										) &
									+  1.0*( &
													r(2*i + 1, 2*j + 1, 2*k + 1) &
												+ r(2*i - 2, 2*j + 1, 2*k + 1) &
												+ r(2*i + 1, 2*j - 2, 2*k + 1) &
												+ r(2*i + 1, 2*j + 1, 2*k - 2) &
												+ r(2*i - 2, 2*j - 2, 2*k + 1) &
												+ r(2*i + 1, 2*j - 2, 2*k - 2) &
												+ r(2*i - 2, 2*j + 1, 2*k - 2) &
												+ r(2*i - 2, 2*j - 2, 2*k - 2) &
										) &
									)
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine mg_restriction_2

subroutine mg_prolongation_2(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)/2+g, 1-g:sz(2)/2+g, 1-g:sz(3)/2+g)	:: xc
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)				:: x
	real										:: a
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	a = 4.0/64.0
!$omp parallel private(i, j, k)
!$omp do 
	do k=1,kx/2
	do j=1,jx/2
	do i=1,ix/2
		x(2*i    , 2*j    , 2*k    ) = x(2*i    , 2*j    , 2*k    ) &
																		+ a*27.0*xc(i, j, k) &
																		+ a*9.0*(xc(i+1, j, k) + xc(i, j+1, k) + xc(i, j, k+1)) &
																		+ a*3.0*(xc(i+1, j+1, k) + xc(i, j+1, k+1) + xc(i+1, j, k+1)) &
																		+ a*xc(i+1, j+1, k+1)
		x(2*i - 1, 2*j    , 2*k    ) = x(2*i - 1, 2*j    , 2*k    ) &
																		+ a*27.0*xc(i, j, k) &
																		+ a*9.0*(xc(i-1, j, k) + xc(i, j+1, k) + xc(i, j, k+1)) &
																		+ a*3.0*(xc(i-1, j+1, k) + xc(i, j+1, k+1) + xc(i-1, j, k+1)) &
																		+ a*xc(i-1, j+1, k+1)
		x(2*i    , 2*j - 1, 2*k    ) = x(2*i    , 2*j - 1, 2*k    ) &
																		+ a*27.0*xc(i, j, k) &
																		+ a*9.0*(xc(i+1, j, k) + xc(i, j-1, k) + xc(i, j, k+1)) &
																		+ a*3.0*(xc(i+1, j-1, k) + xc(i, j-1, k+1) + xc(i+1, j, k+1)) &
																		+ a*xc(i+1, j-1, k+1)
		x(2*i    , 2*j    , 2*k - 1) = x(2*i    , 2*j    , 2*k - 1) &
																		+ a*27.0*xc(i, j, k) &
																		+ a*9.0*(xc(i+1, j, k) + xc(i, j+1, k) + xc(i, j, k-1)) &
																		+ a*3.0*(xc(i+1, j+1, k) + xc(i, j+1, k-1) + xc(i+1, j, k-1)) &
																		+ a*xc(i+1, j+1, k-1)
		x(2*i - 1, 2*j - 1, 2*k    ) = x(2*i - 1, 2*j - 1, 2*k    ) &
																		+ a*27.0*xc(i, j, k) &
																		+ a*9.0*(xc(i-1, j, k) + xc(i, j-1, k) + xc(i, j, k+1)) &
																		+ a*3.0*(xc(i-1, j-1, k) + xc(i, j-1, k+1) + xc(i-1, j, k+1)) &
																		+ a*xc(i-1, j-1, k+1)
		x(2*i    , 2*j - 1, 2*k - 1) = x(2*i    , 2*j - 1, 2*k - 1) &
																		+ a*27.0*xc(i, j, k) &
																		+ a*9.0*(xc(i+1, j, k) + xc(i, j-1, k) + xc(i, j, k-1)) &
																		+ a*3.0*(xc(i+1, j-1, k) + xc(i, j-1, k-1) + xc(i+1, j, k-1)) &
																		+ a*xc(i+1, j-1, k-1)
		x(2*i - 1, 2*j    , 2*k - 1) = x(2*i - 1, 2*j    , 2*k - 1) &
																		+ a*27.0*xc(i, j, k) &
																		+ a*9.0*(xc(i-1, j, k) + xc(i, j+1, k) + xc(i, j, k-1)) &
																		+ a*3.0*(xc(i-1, j+1, k) + xc(i, j+1, k-1) + xc(i-1, j, k-1)) &
																		+ a*xc(i-1, j+1, k-1)
		x(2*i - 1, 2*j - 1, 2*k - 1) = x(2*i - 1, 2*j - 1, 2*k - 1) &
																		+ a*27.0*xc(i, j, k) &
																		+ a*9.0*(xc(i-1, j, k) + xc(i, j-1, k) + xc(i, j, k-1)) &
																		+ a*3.0*(xc(i-1, j-1, k) + xc(i, j-1, k-1) + xc(i-1, j, k-1)) &
																		+ a*xc(i-1, j-1, k-1)
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine mg_prolongation_2

subroutine mg_restriction_1_mp_0(xc, x, phi, psi, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)/2+g, 1-g:sz(2)/2+g, 1-g:sz(3)/2+g)	:: xc
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)				:: x
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)				:: phi, psi
	real										:: psic0
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(psic0)
!$omp do
	do k=1,kx/2
	do j=1,jx/2
	do i=1,ix/2
		xc(i, j, k) = 0.125*( x(2*i    , 2*j    , 2*k    ) &
												+ x(2*i - 1, 2*j    , 2*k    ) &
												+ x(2*i    , 2*j - 1, 2*k    ) &
												+ x(2*i    , 2*j    , 2*k - 1) &
												+ x(2*i - 1, 2*j - 1, 2*k    ) &
												+ x(2*i    , 2*j - 1, 2*k - 1) &
												+ x(2*i - 1, 2*j    , 2*k - 1) &
												+ x(2*i - 1, 2*j - 1, 2*k - 1) )
		psic0 = 0.125*( psi(2*i    , 2*j    , 2*k    ) &
									+ psi(2*i - 1, 2*j    , 2*k    ) &
									+ psi(2*i    , 2*j - 1, 2*k    ) &
									+ psi(2*i    , 2*j    , 2*k - 1) &
									+ psi(2*i - 1, 2*j - 1, 2*k    ) &
									+ psi(2*i    , 2*j - 1, 2*k - 1) &
									+ psi(2*i - 1, 2*j    , 2*k - 1) &
									+ psi(2*i - 1, 2*j - 1, 2*k - 1) )*0.5
		if( psic0 >= 0 ) then
			xc(i, j, k) = 0.0
		endif
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine mg_restriction_1_mp_0

subroutine mg_prolongation_1_mp_0(x, xc, phi, psi, rhol, rhog, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)/2+g, 1-g:sz(2)/2+g, 1-g:sz(3)/2+g)	:: xc
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)				:: x
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)				:: phi, psi
	real										:: rhol, rhog
	real										:: a
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	a = 4.0
!$omp parallel private(i, j, k)
!$omp do
	do k=1,kx/2
	do j=1,jx/2
	do i=1,ix/2
		x(2*i    , 2*j    , 2*k    ) = x(2*i    , 2*j    , 2*k    ) + xc(i, j, k)*a
		x(2*i - 1, 2*j    , 2*k    ) = x(2*i - 1, 2*j    , 2*k    ) + xc(i, j, k)*a
		x(2*i    , 2*j - 1, 2*k    ) = x(2*i    , 2*j - 1, 2*k    ) + xc(i, j, k)*a
		x(2*i    , 2*j    , 2*k - 1) = x(2*i    , 2*j    , 2*k - 1) + xc(i, j, k)*a
		x(2*i - 1, 2*j - 1, 2*k    ) = x(2*i - 1, 2*j - 1, 2*k    ) + xc(i, j, k)*a
		x(2*i    , 2*j - 1, 2*k - 1) = x(2*i    , 2*j - 1, 2*k - 1) + xc(i, j, k)*a
		x(2*i - 1, 2*j    , 2*k - 1) = x(2*i - 1, 2*j    , 2*k - 1) + xc(i, j, k)*a
		x(2*i - 1, 2*j - 1, 2*k - 1) = x(2*i - 1, 2*j - 1, 2*k - 1) + xc(i, j, k)*a
		if( psi(2*i    , 2*j    , 2*k    ) >= 0 ) then
			x(2*i    , 2*j    , 2*k    ) = 0.0
		endif
		if( psi(2*i - 1, 2*j    , 2*k    ) >= 0 ) then
			x(2*i - 1, 2*j    , 2*k    ) = 0.0
		endif
		if( psi(2*i    , 2*j - 1, 2*k    ) >= 0 ) then
			x(2*i    , 2*j - 1, 2*k    ) = 0.0
		endif
		if( psi(2*i    , 2*j    , 2*k - 1) >= 0 ) then
			x(2*i    , 2*j    , 2*k - 1) = 0.0
		endif
		if( psi(2*i - 1, 2*j - 1, 2*k    ) >= 0 ) then
			x(2*i - 1, 2*j - 1, 2*k    ) = 0.0
		endif
		if( psi(2*i    , 2*j - 1, 2*k - 1) >= 0 ) then
			x(2*i    , 2*j - 1, 2*k - 1) = 0.0
		endif
		if( psi(2*i - 1, 2*j    , 2*k - 1) >= 0 ) then
			x(2*i - 1, 2*j    , 2*k - 1) = 0.0
		endif
		if( psi(2*i - 1, 2*j - 1, 2*k - 1) >= 0 ) then
			x(2*i - 1, 2*j - 1, 2*k - 1) = 0.0
		endif
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine mg_prolongation_1_mp_0

subroutine mg_restriction_1_mp(xc, x, phi, psi, rhol, rhog, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)/2+g, 1-g:sz(2)/2+g, 1-g:sz(3)/2+g)	:: xc
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)				:: x
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)				:: phi, psi
	real										:: rhol, rhog
	real										:: x_ppp, x_ppn, x_pnp, x_npp, x_pnn, x_npn, x_nnp, x_nnn
	real										:: phi_ppp, phi_ppn, phi_pnp, phi_npp, phi_pnn, phi_npn, phi_nnp, phi_nnn
	real										:: psi_ppp, psi_ppn, psi_pnp, psi_npp, psi_pnn, psi_npn, psi_nnp, psi_nnn
	real										:: Hphi_ppp, Hphi_ppn, Hphi_pnp, Hphi_npp, Hphi_pnn, Hphi_npn, Hphi_nnp, Hphi_nnn
	real										:: rho_ppp, rho_ppn, rho_pnp, rho_npp, rho_pnn, rho_npn, rho_nnp, rho_nnn
	real										:: xc_a, xc_b
	real										:: phic0, psic0
	real										:: mg_geth
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(x_ppp, x_ppn, x_pnp, x_npp, x_pnn, x_npn, x_nnp, x_nnn) &
!$omp					,private(phi_ppp, phi_ppn, phi_pnp, phi_npp, phi_pnn, phi_npn, phi_nnp, phi_nnn) &
!$omp					,private(psi_ppp, psi_ppn, psi_pnp, psi_npp, psi_pnn, psi_npn, psi_nnp, psi_nnn) &
!$omp					,private(Hphi_ppp, Hphi_ppn, Hphi_pnp, Hphi_npp, Hphi_pnn, Hphi_npn, Hphi_nnp, Hphi_nnn) &
!$omp					,private(rho_ppp, rho_ppn, rho_pnp, rho_npp, rho_pnn, rho_npn, rho_nnp, rho_nnn) &
!$omp					,private(xc_a, xc_b) &
!$omp					,private(phic0, psic0)
!$omp do
	do k=1,kx/2
	do j=1,jx/2
	do i=1,ix/2
		phi_ppp = phi(2*i, 2*j, 2*k)
		phi_ppn = phi(2*i, 2*j, 2*k-1)
		phi_pnp = phi(2*i, 2*j-1, 2*k)
		phi_npp = phi(2*i-1, 2*j, 2*k)
		phi_pnn = phi(2*i, 2*j-1, 2*k-1)
		phi_npn = phi(2*i-1, 2*j, 2*k-1)
		phi_nnp = phi(2*i-1, 2*j-1, 2*k)
		phi_nnn = phi(2*i-1, 2*j-1, 2*k-1)

		psi_ppp = psi(2*i, 2*j, 2*k)
		psi_ppn = psi(2*i, 2*j, 2*k-1)
		psi_pnp = psi(2*i, 2*j-1, 2*k)
		psi_npp = psi(2*i-1, 2*j, 2*k)
		psi_pnn = psi(2*i, 2*j-1, 2*k-1)
		psi_npn = psi(2*i-1, 2*j, 2*k-1)
		psi_nnp = psi(2*i-1, 2*j-1, 2*k)
		psi_nnn = psi(2*i-1, 2*j-1, 2*k-1)

		Hphi_ppp = mg_geth(phi_ppp)
		Hphi_ppn = mg_geth(phi_ppn)
		Hphi_pnp = mg_geth(phi_pnp)
		Hphi_npp = mg_geth(phi_npp)
		Hphi_pnn = mg_geth(phi_pnn)
		Hphi_npn = mg_geth(phi_npn)
		Hphi_nnp = mg_geth(phi_nnp)
		Hphi_nnn = mg_geth(phi_nnn)

		rho_ppp = rhog + (rhol - rhog)*Hphi_ppp
		rho_ppn = rhog + (rhol - rhog)*Hphi_ppn
		rho_pnp = rhog + (rhol - rhog)*Hphi_pnp
		rho_npp = rhog + (rhol - rhog)*Hphi_npp
		rho_pnn = rhog + (rhol - rhog)*Hphi_pnn
		rho_npn = rhog + (rhol - rhog)*Hphi_npn
		rho_nnp = rhog + (rhol - rhog)*Hphi_nnp
		rho_nnn = rhog + (rhol - rhog)*Hphi_nnn

		phic0 = 0.125*( phi_ppp + phi_nnn &
									+ phi_ppn + phi_nnp &
									+ phi_pnp + phi_npn &
									+ phi_npp + phi_pnn )*0.5
		psic0 = 0.125*( psi_ppp + psi_nnn &
									+ psi_ppn + psi_nnp &
									+ psi_pnp + psi_npn &
									+ psi_npp + psi_pnn )*0.5

		x_ppp = x(2*i  , 2*j  , 2*k  )
		x_ppn = x(2*i  , 2*j  , 2*k-1)
		x_pnp = x(2*i  , 2*j-1, 2*k  )
		x_npp = x(2*i-1, 2*j  , 2*k  )
		x_pnn = x(2*i  , 2*j-1, 2*k-1)
		x_npn = x(2*i-1, 2*j  , 2*k-1)
		x_nnp = x(2*i-1, 2*j-1, 2*k  )
		x_nnn = x(2*i-1, 2*j-1, 2*k-1)

		xc_a = (x_ppp + x_nnn &
					+ x_ppn + x_nnp &
					+ x_pnp + x_npn &
					+ x_npp + x_pnn )
		xc_b = 8.0

		xc_a = (x_ppp/rho_ppp + x_nnn/rho_nnn &
					+ x_ppn/rho_ppn + x_nnp/rho_nnp &
					+ x_pnp/rho_pnp + x_npn/rho_npn &
					+ x_npp/rho_npp + x_pnn/rho_pnn )
		xc_b = (  1.0/rho_ppp +   1.0/rho_nnn &
					+   1.0/rho_ppn +   1.0/rho_nnp &
					+   1.0/rho_pnp +   1.0/rho_npn &
					+   1.0/rho_npp +   1.0/rho_pnn )

		xc_a = (x_ppp*rho_ppp + x_nnn*rho_nnn &
					+ x_ppn*rho_ppn + x_nnp*rho_nnp &
					+ x_pnp*rho_pnp + x_npn*rho_npn &
					+ x_npp*rho_npp + x_pnn*rho_pnn )
		xc_b = (  1.0*rho_ppp +   1.0*rho_nnn &
					+   1.0*rho_ppn +   1.0*rho_nnp &
					+   1.0*rho_pnp +   1.0*rho_npn &
					+   1.0*rho_npp +   1.0*rho_pnn )

		xc_a = (x_ppp + x_nnn &
					+ x_ppn + x_nnp &
					+ x_pnp + x_npn &
					+ x_npp + x_pnn )
		xc_b = 8.0

		xc(i, j, k) = xc_a/xc_b

		if( psic0 >= 0 ) then
			xc(i, j, k) = 0.0
		endif

	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine mg_restriction_1_mp

subroutine mg_prolongation_1_mp(x, xc, phi, psi, rhol, rhog, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)/2+g, 1-g:sz(2)/2+g, 1-g:sz(3)/2+g)	:: xc
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)				:: x
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)				:: phi, psi
	real										:: rhol, rhog
	real										:: a
	real										:: x_ppp, x_ppn, x_pnp, x_npp, x_pnn, x_npn, x_nnp, x_nnn
	real										:: phi_ppp, phi_ppn, phi_pnp, phi_npp, phi_pnn, phi_npn, phi_nnp, phi_nnn
	real										:: psi_ppp, psi_ppn, psi_pnp, psi_npp, psi_pnn, psi_npn, psi_nnp, psi_nnn
	real										:: Hphi_ppp, Hphi_ppn, Hphi_pnp, Hphi_npp, Hphi_pnn, Hphi_npn, Hphi_nnp, Hphi_nnn
	real										:: rho_ppp, rho_ppn, rho_pnp, rho_npp, rho_pnn, rho_npn, rho_nnp, rho_nnn
	real										:: dx_ppp, dx_ppn, dx_pnp, dx_npp, dx_pnn, dx_npn, dx_nnp, dx_nnn
	real										:: xc_a, xc_b
	real										:: phic0, psic0
	real										:: mg_geth
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	a = 4.0
!$omp parallel private(i, j, k) &
!$omp					,private(x_ppp, x_ppn, x_pnp, x_npp, x_pnn, x_npn, x_nnp, x_nnn) &
!$omp					,private(phi_ppp, phi_ppn, phi_pnp, phi_npp, phi_pnn, phi_npn, phi_nnp, phi_nnn) &
!$omp					,private(psi_ppp, psi_ppn, psi_pnp, psi_npp, psi_pnn, psi_npn, psi_nnp, psi_nnn) &
!$omp					,private(Hphi_ppp, Hphi_ppn, Hphi_pnp, Hphi_npp, Hphi_pnn, Hphi_npn, Hphi_nnp, Hphi_nnn) &
!$omp					,private(rho_ppp, rho_ppn, rho_pnp, rho_npp, rho_pnn, rho_npn, rho_nnp, rho_nnn) &
!$omp					,private(dx_ppp, dx_ppn, dx_pnp, dx_npp, dx_pnn, dx_npn, dx_nnp, dx_nnn) &
!$omp					,private(xc_a, xc_b) &
!$omp					,private(phic0, psic0)
!$omp do
	do k=1,kx/2
	do j=1,jx/2
	do i=1,ix/2
		phi_ppp = phi(2*i, 2*j, 2*k)
		phi_ppn = phi(2*i, 2*j, 2*k-1)
		phi_pnp = phi(2*i, 2*j-1, 2*k)
		phi_npp = phi(2*i-1, 2*j, 2*k)
		phi_pnn = phi(2*i, 2*j-1, 2*k-1)
		phi_npn = phi(2*i-1, 2*j, 2*k-1)
		phi_nnp = phi(2*i-1, 2*j-1, 2*k)
		phi_nnn = phi(2*i-1, 2*j-1, 2*k-1)

		psi_ppp = psi(2*i, 2*j, 2*k)
		psi_ppn = psi(2*i, 2*j, 2*k-1)
		psi_pnp = psi(2*i, 2*j-1, 2*k)
		psi_npp = psi(2*i-1, 2*j, 2*k)
		psi_pnn = psi(2*i, 2*j-1, 2*k-1)
		psi_npn = psi(2*i-1, 2*j, 2*k-1)
		psi_nnp = psi(2*i-1, 2*j-1, 2*k)
		psi_nnn = psi(2*i-1, 2*j-1, 2*k-1)

		Hphi_ppp = mg_geth(phi_ppp)
		Hphi_ppn = mg_geth(phi_ppn)
		Hphi_pnp = mg_geth(phi_pnp)
		Hphi_npp = mg_geth(phi_npp)
		Hphi_pnn = mg_geth(phi_pnn)
		Hphi_npn = mg_geth(phi_npn)
		Hphi_nnp = mg_geth(phi_nnp)
		Hphi_nnn = mg_geth(phi_nnn)

		rho_ppp = rhog + (rhol - rhog)*Hphi_ppp
		rho_ppn = rhog + (rhol - rhog)*Hphi_ppn
		rho_pnp = rhog + (rhol - rhog)*Hphi_pnp
		rho_npp = rhog + (rhol - rhog)*Hphi_npp
		rho_pnn = rhog + (rhol - rhog)*Hphi_pnn
		rho_npn = rhog + (rhol - rhog)*Hphi_npn
		rho_nnp = rhog + (rhol - rhog)*Hphi_nnp
		rho_nnn = rhog + (rhol - rhog)*Hphi_nnn

		phic0 = 0.125*( phi_ppp + phi_nnn &
									+ phi_ppn + phi_nnp &
									+ phi_pnp + phi_npn &
									+ phi_npp + phi_pnn )*0.5
		psic0 = 0.125*( psi_ppp + psi_nnn &
									+ psi_ppn + psi_nnp &
									+ psi_pnp + psi_npn &
									+ psi_npp + psi_pnn )*0.5

		x_ppp = x(2*i  , 2*j  , 2*k  )
		x_ppn = x(2*i  , 2*j  , 2*k-1)
		x_pnp = x(2*i  , 2*j-1, 2*k  )
		x_npp = x(2*i-1, 2*j  , 2*k  )
		x_pnn = x(2*i  , 2*j-1, 2*k-1)
		x_npn = x(2*i-1, 2*j  , 2*k-1)
		x_nnp = x(2*i-1, 2*j-1, 2*k  )
		x_nnn = x(2*i-1, 2*j-1, 2*k-1)

		xc_a = (x_ppp/rho_ppp + x_nnn/rho_nnn &
					+ x_ppn/rho_ppn + x_nnp/rho_nnp &
					+ x_pnp/rho_pnp + x_npn/rho_npn &
					+ x_npp/rho_npp + x_pnn/rho_pnn )
		xc_b = (  1.0/rho_ppp +   1.0/rho_nnn &
					+   1.0/rho_ppn +   1.0/rho_nnp &
					+   1.0/rho_pnp +   1.0/rho_npn &
					+   1.0/rho_npp +   1.0/rho_pnn )

		dx_ppp = xc(i, j, k)
		dx_ppn = xc(i, j, k)
		dx_pnp = xc(i, j, k)
		dx_npp = xc(i, j, k)
		dx_pnn = xc(i, j, k)
		dx_npn = xc(i, j, k)
		dx_nnp = xc(i, j, k)
		dx_nnn = xc(i, j, k)

		dx_ppp = xc(i, j, k)/rho_ppp/xc_b*8.0
		dx_ppn = xc(i, j, k)/rho_ppn/xc_b*8.0
		dx_pnp = xc(i, j, k)/rho_pnp/xc_b*8.0
		dx_npp = xc(i, j, k)/rho_npp/xc_b*8.0
		dx_pnn = xc(i, j, k)/rho_pnn/xc_b*8.0
		dx_npn = xc(i, j, k)/rho_npn/xc_b*8.0
		dx_nnp = xc(i, j, k)/rho_nnp/xc_b*8.0
		dx_nnn = xc(i, j, k)/rho_nnn/xc_b*8.0

		dx_ppp = xc(i, j, k)
		dx_ppn = xc(i, j, k)
		dx_pnp = xc(i, j, k)
		dx_npp = xc(i, j, k)
		dx_pnn = xc(i, j, k)
		dx_npn = xc(i, j, k)
		dx_nnp = xc(i, j, k)
		dx_nnn = xc(i, j, k)

		if( psi_ppp >= 0 ) then
			dx_ppp = 0.0
		endif
		if( psi_ppn >= 0 ) then
			dx_ppn = 0.0
		endif
		if( psi_pnp >= 0 ) then
			dx_pnp = 0.0
		endif
		if( psi_npp >= 0 ) then
			dx_npp = 0.0
		endif
		if( psi_pnn >= 0 ) then
			dx_pnn = 0.0
		endif
		if( psi_npn >= 0 ) then
			dx_npn = 0.0
		endif
		if( psi_nnp >= 0 ) then
			dx_nnp = 0.0
		endif
		if( psi_nnn >= 0 ) then
			dx_nnn = 0.0
		endif

		x_ppp = x_ppp + dx_ppp*a
		x_ppn = x_ppn + dx_ppn*a
		x_pnp = x_pnp + dx_pnp*a
		x_npp = x_npp + dx_npp*a
		x_pnn = x_pnn + dx_npn*a
		x_npn = x_npn + dx_npn*a
		x_nnp = x_nnp + dx_nnp*a
		x_nnn = x_nnn + dx_nnn*a

		x(2*i  , 2*j  , 2*k  ) = x_ppp 
		x(2*i  , 2*j  , 2*k-1) = x_ppn
		x(2*i  , 2*j-1, 2*k  ) = x_pnp
		x(2*i-1, 2*j  , 2*k  ) = x_npp
		x(2*i  , 2*j-1, 2*k-1) = x_pnn
		x(2*i-1, 2*j  , 2*k-1) = x_npn
		x(2*i-1, 2*j-1, 2*k  ) = x_nnp
		x(2*i-1, 2*j-1, 2*k-1) = x_nnn
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine mg_prolongation_1_mp






subroutine mg_restriction_1_1(bc, r, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)/2+g, 1-g:sz(2)/2+g, 1-g:sz(3)/2+g)	:: bc
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)				:: r
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j, k)
!$omp do
	do k=2,kx,2
	do j=2,jx,2
	do i=2,ix,2
		bc(i/2, j/2, k/2) = 0.125*( r(i    , j    , k    ) &
															+ r(i - 1, j    , k    ) &
															+ r(i    , j - 1, k    ) &
															+ r(i    , j    , k - 1) &
															+ r(i - 1, j - 1, k    ) &
															+ r(i    , j - 1, k - 1) &
															+ r(i - 1, j    , k - 1) &
															+ r(i - 1, j - 1, k - 1) )
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine mg_restriction_1_1

subroutine mg_prolongation_1_1(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)/2+g, 1-g:sz(2)/2+g, 1-g:sz(3)/2+g)	:: xc
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)				:: x
	real										:: a
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	a = 4.0
!$omp parallel private(i, j, k)
!$omp do
	do k=2,kx,2
	do j=2,jx,2
	do i=2,ix,2
		x(i    , j    , k    ) = x(i    , j    , k    ) + xc(i/2, j/2, k/2)*a
		x(i - 1, j    , k    ) = x(i - 1, j    , k    ) + xc(i/2, j/2, k/2)*a
		x(i    , j - 1, k    ) = x(i    , j - 1, k    ) + xc(i/2, j/2, k/2)*a
		x(i    , j    , k - 1) = x(i    , j    , k - 1) + xc(i/2, j/2, k/2)*a
		x(i - 1, j - 1, k    ) = x(i - 1, j - 1, k    ) + xc(i/2, j/2, k/2)*a
		x(i    , j - 1, k - 1) = x(i    , j - 1, k - 1) + xc(i/2, j/2, k/2)*a
		x(i - 1, j    , k - 1) = x(i - 1, j    , k - 1) + xc(i/2, j/2, k/2)*a
		x(i - 1, j - 1, k - 1) = x(i - 1, j - 1, k - 1) + xc(i/2, j/2, k/2)*a
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine mg_prolongation_1_1

subroutine mg_restriction_1_2(bc, r, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)/2+g, 1-g:sz(2)/2+g, 1-g:sz(3)/2+g)	:: bc
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)				:: r
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j, k)
!$omp do schedule(static, 1)
	do k=1,kx
	do j=1,jx
	do i=1,ix
		bc((i+1)/2, (j+1)/2, (k+1)/2) = bc((i+1)/2, (j+1)/2, (k+1)/2) + 0.125*r(i, j, k)
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine mg_restriction_1_2

subroutine mg_prolongation_1_2(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)/2+g, 1-g:sz(2)/2+g, 1-g:sz(3)/2+g)	:: xc
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)				:: x
	real										:: a
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	a = 4.0
!$omp parallel private(i, j, k)
!$omp do schedule(static, 1)
	do k=1,kx
	do j=1,jx
	do i=1,ix
		x(i, j, k) = x(i, j, k) + xc((i+1)/2, (j+1)/2, (k+1)/2)*a
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine mg_prolongation_1_2


subroutine mg_restriction2(Ac, bc, A, r, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)/2+g, 1-g:sz(2)/2+g, 1-g:sz(3)/2+g, 0:3)	:: Ac
	real, dimension(1-g:sz(1)/2+g, 1-g:sz(2)/2+g, 1-g:sz(3)/2+g)			:: bc
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:3)				:: A
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)						:: r
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j, k)
!$omp do 
	do k=1,kx/2
	do j=1,jx/2
	do i=1,ix/2
		bc(i, j, k) = 0.125*( r(2*i    , 2*j    , 2*k    ) &
												+ r(2*i + 1, 2*j    , 2*k    ) &
												+ r(2*i    , 2*j + 1, 2*k    ) &
												+ r(2*i    , 2*j    , 2*k + 1) &
												+ r(2*i + 1, 2*j + 1, 2*k    ) &
												+ r(2*i    , 2*j + 1, 2*k + 1) &
												+ r(2*i + 1, 2*j    , 2*k + 1) &
												+ r(2*i + 1, 2*j + 1, 2*k + 1) )
		Ac(i, j, k, 1) = 0.25*( &
													A(2*i + 1, 2*j    , 2*k    , 1) &
												+ A(2*i + 1, 2*j + 1, 2*k    , 1) &
												+ A(2*i + 1, 2*j    , 2*k + 1, 1) &
												+ A(2*i + 1, 2*j + 1, 2*k + 1, 1) )
		Ac(i, j, k, 2) = 0.25*( &
													A(2*i    , 2*j + 1, 2*k    , 2) &
												+ A(2*i + 1, 2*j + 1, 2*k    , 2) &
												+ A(2*i    , 2*j + 1, 2*k + 1, 2) &
												+ A(2*i + 1, 2*j + 1, 2*k + 1, 2) )
		Ac(i, j, k, 3) = 0.25*( &
													A(2*i    , 2*j    , 2*k + 1, 3) &
												+ A(2*i    , 2*j + 1, 2*k + 1, 3) &
												+ A(2*i + 1, 2*j    , 2*k + 1, 3) &
												+ A(2*i + 1, 2*j + 1, 2*k + 1, 3) )
		Ac(i, j, k, 0) = 0.125*( &
													A(2*i    , 2*j    , 2*k    , 0) &
												+ A(2*i + 1, 2*j    , 2*k    , 0) &
												+ A(2*i    , 2*j + 1, 2*k    , 0) &
												+ A(2*i    , 2*j    , 2*k + 1, 0) &
												+ A(2*i + 1, 2*j + 1, 2*k    , 0) &
												+ A(2*i    , 2*j + 1, 2*k + 1, 0) &
												+ A(2*i + 1, 2*j    , 2*k + 1, 0) &
												+ A(2*i + 1, 2*j + 1, 2*k + 1, 0) )
		Ac(i, j, k, 0) = - Ac(i, j, k, 1) - Ac(i-1, j, k, 1) &
										 - Ac(i, j, k, 2) - Ac(i, j-1, k, 2) &
										 - Ac(i, j, k, 3) - Ac(i, j, k-1, 3)
!		Ac(1, i, j, k) = 0.125*( &
!													A(1, 2*i    , 2*j    , 2*k    ) &
!												+ A(1, 2*i + 1, 2*j    , 2*k    ) &
!												+ A(1, 2*i    , 2*j + 1, 2*k    ) &
!												+ A(1, 2*i    , 2*j    , 2*k + 1) &
!												+ A(1, 2*i + 1, 2*j + 1, 2*k    ) &
!												+ A(1, 2*i    , 2*j + 1, 2*k + 1) &
!												+ A(1, 2*i + 1, 2*j    , 2*k + 1) &
!												+ A(1, 2*i + 1, 2*j + 1, 2*k + 1) )
!		Ac(2, i, j, k) = 0.125*( &
!													A(2, 2*i    , 2*j    , 2*k    ) &
!												+ A(2, 2*i + 1, 2*j    , 2*k    ) &
!												+ A(2, 2*i    , 2*j + 1, 2*k    ) &
!												+ A(2, 2*i    , 2*j    , 2*k + 1) &
!												+ A(2, 2*i + 1, 2*j + 1, 2*k    ) &
!												+ A(2, 2*i    , 2*j + 1, 2*k + 1) &
!												+ A(2, 2*i + 1, 2*j    , 2*k + 1) &
!												+ A(2, 2*i + 1, 2*j + 1, 2*k + 1) )
!		Ac(3, i, j, k) = 0.125*( &
!													A(3, 2*i    , 2*j    , 2*k    ) &
!												+ A(3, 2*i + 1, 2*j    , 2*k    ) &
!												+ A(3, 2*i    , 2*j + 1, 2*k    ) &
!												+ A(3, 2*i    , 2*j    , 2*k + 1) &
!												+ A(3, 2*i + 1, 2*j + 1, 2*k    ) &
!												+ A(3, 2*i    , 2*j + 1, 2*k + 1) &
!												+ A(3, 2*i + 1, 2*j    , 2*k + 1) &
!												+ A(3, 2*i + 1, 2*j + 1, 2*k + 1) )
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine mg_restriction2

subroutine mg_restriction_1_ls(xc, x, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)/2+g, 1-g:sz(2)/2+g, 1-g:sz(3)/2+g)	:: xc
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)				:: x
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j, k)
!$omp do
	do k=1,kx/2
	do j=1,jx/2
	do i=1,ix/2
		xc(i, j, k) = 0.125*( x(2*i    , 2*j    , 2*k    ) &
												+ x(2*i - 1, 2*j    , 2*k    ) &
												+ x(2*i    , 2*j - 1, 2*k    ) &
												+ x(2*i    , 2*j    , 2*k - 1) &
												+ x(2*i - 1, 2*j - 1, 2*k    ) &
												+ x(2*i    , 2*j - 1, 2*k - 1) &
												+ x(2*i - 1, 2*j    , 2*k - 1) &
												+ x(2*i - 1, 2*j - 1, 2*k - 1) )*0.5
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine mg_restriction_1_ls


