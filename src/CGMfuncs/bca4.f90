!>  @file  bca4.f90
!!  @brief Functions to compute the elements of the coefficient matrix at the outer boundaries of the simulation box for Cartesian grid data structure (for symmetric 7-band matrices)
!<

subroutine bc_A1_d_4(A, b, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:3)	:: A
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)			:: b
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i = sz(1)
!$omp parallel private(j, k)
!$omp do
	do k=1, kx
	do j=1, jx
		A(i, j, k, 0) = A(i, j, k, 0) - A(i, j, k, 1)
		b(i, j, k) = b(i, j, k) - A(i, j, k, 1)*xc*2.0
		A(i, j, k, 1) = 0.0
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_A1_d_4

subroutine bc_A3_d_4(A, b, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:3)	:: A
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)			:: b
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i = 1
!$omp parallel private(j, k)
!$omp do
	do k=1, kx
	do j=1, jx
		A(i, j, k, 0) = A(i, j, k, 0) - A(i-1, j, k, 1)
		b(i, j, k) = b(i, j, k) - A(i-1, j, k, 1)*xc*2.0
		A(i-1, j, k, 1) = 0.0
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_A3_d_4

subroutine bc_A2_d_4(A, b, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:3)	:: A
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)			:: b
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	j = sz(2)
!$omp parallel private(i, k)
!$omp do
	do k=1, kx
	do i=1, ix
		A(i, j, k, 0) = A(i, j, k, 0) - A(i, j, k, 2)
		b(i, j, k) = b(i, j, k) - A(i, j, k, 2)*xc*2.0
		A(i, j, k, 2) = 0.0
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_A2_d_4

subroutine bc_A4_d_4(A, b, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:3)	:: A
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)			:: b
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	j = 1
!$omp parallel private(i, k)
!$omp do
	do k=1, kx
	do i=1, ix
		A(i, j, k, 0) = A(i, j, k, 0) - A(i, j-1, k, 2)
		b(i, j, k) = b(i, j, k) - A(i, j-1, k, 2)*xc*2.0
		A(i, j-1, k, 2) = 0.0
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_A4_d_4

subroutine bc_A5_d_4(A, b, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:3)	:: A
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)			:: b
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	k = sz(3)
!$omp parallel private(i, j)
!$omp do
	do j=1, jx
	do i=1, ix
		A(i, j, k, 0) = A(i, j, k, 0) - A(i, j, k, 3)
		b(i, j, k) = b(i, j, k) - A(i, j, k, 3)*xc*2.0
		A(i, j, k, 3) = 0.0
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_A5_d_4

subroutine bc_A6_d_4(A, b, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:3)	:: A
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)			:: b
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	k = 1
!$omp parallel private(i, j)
!$omp do
	do j=1, jx
	do i=1, ix
		A(i, j, k, 0) = A(i, j, k, 0) - A(i, j, k-1, 3)
		b(i, j, k) = b(i, j, k) - A(i, j, k-1, 3)*xc*2.0
		A(i, j, k-1, 3) = 0.0
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_A6_d_4

subroutine bc_A1_n_4(A, b, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:3)	:: A
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)			:: b
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i = sz(1)
!$omp parallel private(j, k)
!$omp do
	do k=1, kx
	do j=1, jx
		A(i, j, k, 0) = A(i, j, k, 0) + A(i, j, k, 1)
		b(i, j, k) = b(i, j, k) - A(i, j, k, 1)*xc
		A(i, j, k, 1) = 0.0
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_A1_n_4

subroutine bc_A3_n_4(A, b, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:3)	:: A
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)			:: b
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i = 1
!$omp parallel private(j, k)
!$omp do
	do k=1, kx
	do j=1, jx
		A(i, j, k, 0) = A(i, j, k, 0) + A(i-1, j, k, 1)
		b(i, j, k) = b(i, j, k) + A(i-1, j, k, 1)*xc
		A(i-1, j, k, 1) = 0.0
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_A3_n_4

subroutine bc_A2_n_4(A, b, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:3)	:: A
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)			:: b
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	j = sz(2)
!$omp parallel private(i, k)
!$omp do
	do k=1, kx
	do i=1, ix
		A(i, j, k, 0) = A(i, j, k, 0) + A(i, j, k, 2)
		b(i, j, k) = b(i, j, k) - A(i, j, k, 2)*xc
		A(i, j, k, 2) = 0.0
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_A2_n_4

subroutine bc_A4_n_4(A, b, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:3)	:: A
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)			:: b
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	j = 1
!$omp parallel private(i, k)
!$omp do
	do k=1, kx
	do i=1, ix
		A(i, j, k, 0) = A(i, j, k, 0) + A(i, j-1, k, 2)
		b(i, j, k) = b(i, j, k) + A(i, j-1, k, 2)*xc
		A(i, j-1, k, 2) = 0.0
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_A4_n_4

subroutine bc_A5_n_4(A, b, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:3)	:: A
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)			:: b
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	k = sz(3)
!$omp parallel private(i, j)
!$omp do
	do j=1, jx
	do i=1, ix
		A(i, j, k, 0) = A(i, j, k, 0) + A(i, j, k, 3)
		b(i, j, k) = b(i, j, k) - A(i, j, k, 3)*xc
		A(i, j, k, 3) = 0.0
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_A5_n_4

subroutine bc_A6_n_4(A, b, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:3)	:: A
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)			:: b
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	k = 1
!$omp parallel private(i, j)
!$omp do
	do j=1, jx
	do i=1, ix
		A(i, j, k, 0) = A(i, j, k, 0) + A(i, j, k-1, 3)
		b(i, j, k) = b(i, j, k) + A(i, j, k-1, 3)*xc
		A(i, j, k-1, 3) = 0.0
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_A6_n_4

