!>  @file  bcx.f90
!!  @brief Functions to compute the values at the guide cells (BCX) at the outer boundaries of the simulation box for Cartesian grid data structure
!< 

subroutine bc_x1_d(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i = sz(1)
!$omp parallel private(j, k)
!$omp do
	do k=1-g, kx+g
	do j=1-g, jx+g
		x(i+1, j, k) = 2.0*xc - x(i, j, k)
		x(i+2, j, k) = xc
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x1_d

subroutine bc_x3_d(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i = 1
!$omp parallel private(j, k)
!$omp do
	do k=1-g, kx+g
	do j=1-g, jx+g
		x(i-1, j, k) = 2.0*xc - x(i, j, k)
		x(i-2, j, k) = xc
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x3_d

subroutine bc_x2_d(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	j = sz(2)
!$omp parallel private(i, k)
!$omp do
	do k=1-g, kx+g
	do i=1-g, ix+g
		x(i, j+1, k) = 2.0*xc - x(i, j, k)
		x(i, j+2, k) = xc
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x2_d

subroutine bc_x4_d(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	j = 1
!$omp parallel private(i, k)
!$omp do
	do k=1-g, kx+g
	do i=1-g, ix+g
		x(i, j-1, k) = 2.0*xc - x(i, j, k)
		x(i, j-2, k) = xc
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x4_d

subroutine bc_x5_d(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	k = kx
!$omp parallel private(i, j)
!$omp do
	do j=1-g, jx+g
	do i=1-g, ix+g
		x(i, j, k+1) = 2.0*xc - x(i, j, k)
		x(i, j, k+2) = xc
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x5_d

subroutine bc_x6_d(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	k = 1
!$omp parallel private(i, j)
!$omp do
	do j=1-g, jx+g
	do i=1-g, ix+g
		x(i, j, k-1) = 2.0*xc - x(i, j, k)
		x(i, j, k-2) = xc
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x6_d


subroutine bc_x1_d_2nd(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i = sz(1)
!$omp parallel private(j, k)
!$omp do
	do k=1-g, kx+g
	do j=1-g, jx+g
		x(i+1, j, k) = (8.0*xc - 6.0*x(i, j, k) + x(i-1, j, k))/3.0
		x(i+2, j, k) = xc
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x1_d_2nd

subroutine bc_x3_d_2nd(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i = 1
!$omp parallel private(j, k)
!$omp do
	do k=1-g, kx+g
	do j=1-g, jx+g
		x(i-1, j, k) = (8.0*xc - 6.0*x(i, j, k) + x(i+1, j, k))/3.0
		x(i-2, j, k) = xc
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x3_d_2nd

subroutine bc_x2_d_2nd(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	j = sz(2)
!$omp parallel private(i, k)
!$omp do
	do k=1-g, kx+g
	do i=1-g, ix+g
		x(i, j+1, k) = (8.0*xc - 6.0*x(i, j, k) + x(i, j-1, k))/3.0
		x(i, j+2, k) = xc
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x2_d_2nd

subroutine bc_x4_d_2nd(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	j = 1
!$omp parallel private(i, k)
!$omp do
	do k=1-g, kx+g
	do i=1-g, ix+g
		x(i, j-1, k) = (8.0*xc - 6.0*x(i, j, k) + x(i, j+1, k))/3.0
		x(i, j-2, k) = xc
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x4_d_2nd

subroutine bc_x5_d_2nd(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	k = kx
!$omp parallel private(i, j)
!$omp do
	do j=1-g, jx+g
	do i=1-g, ix+g
		x(i, j, k+1) = (8.0*xc - 6.0*x(i, j, k) + x(i, j, k-1))/3.0
		x(i, j, k+2) = xc
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x5_d_2nd

subroutine bc_x6_d_2nd(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	k = 1
!$omp parallel private(i, j)
!$omp do
	do j=1-g, jx+g
	do i=1-g, ix+g
		x(i, j, k-1) = (8.0*xc - 6.0*x(i, j, k) + x(i, j, k+1))/3.0
		x(i, j, k-2) = xc
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x6_d_2nd

subroutine bc_x1_n(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i = sz(1)
!$omp parallel private(j, k)
!$omp do
	do k=1-g, kx+g
	do j=1-g, jx+g
		x(i+1, j, k) = x(i, j, k) + xc
		x(i+2, j, k) = x(i, j, k) + xc*2.0
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x1_n

subroutine bc_x3_n(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i = 1
!$omp parallel private(j, k)
!$omp do
	do k=1-g, kx+g
	do j=1-g, jx+g
		x(i-1, j, k) = x(i, j, k) - xc
		x(i-2, j, k) = x(i, j, k) - xc*2.0
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x3_n

subroutine bc_x2_n(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	j = sz(2)
!$omp parallel private(i, k)
!$omp do
	do k=1-g, kx+g
	do i=1-g, ix+g
		x(i, j+1, k) = x(i, j, k) + xc
		x(i, j+2, k) = x(i, j, k) + xc*2.0
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x2_n

subroutine bc_x4_n(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	j = 1
!$omp parallel private(i, k)
!$omp do
	do k=1-g, kx+g
	do i=1-g, ix+g
		x(i, j-1, k) = x(i, j, k) - xc
		x(i, j-2, k) = x(i, j, k) - xc*2.0
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x4_n

subroutine bc_x5_n(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	k = sz(3)
!$omp parallel private(i, j)
!$omp do
	do j=1-g, jx+g
	do i=1-g, ix+g
		x(i, j, k+1) = x(i, j, k) + xc
		x(i, j, k+2) = x(i, j, k) + xc*2.0
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x5_n

subroutine bc_x6_n(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	k = 1
!$omp parallel private(i, j)
!$omp do
	do j=1-g, jx+g
	do i=1-g, ix+g
		x(i, j, k-1) = x(i, j, k) - xc
		x(i, j, k-2) = x(i, j, k) - xc*2.0
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x6_n

subroutine bc_x1_d_f(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i = sz(1)
!$omp parallel private(j, k)
!$omp do
	do k=1-g, kx+g
	do j=1-g, jx+g
		x(i+1, j, k) = xc
		x(i+2, j, k) = xc
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x1_d_f

subroutine bc_x3_d_f(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i = 1
!$omp parallel private(j, k)
!$omp do
	do k=1-g, kx+g
	do j=1-g, jx+g
		x(i  , j, k) = xc
		x(i-1, j, k) = xc
		x(i-2, j, k) = xc
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x3_d_f

subroutine bc_x2_d_f(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	j = sz(2)
!$omp parallel private(i, k)
!$omp do
	do k=1-g, kx+g
	do i=1-g, ix+g
		x(i, j+1, k) = xc
		x(i, j+2, k) = xc
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x2_d_f

subroutine bc_x4_d_f(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	j = 1
!$omp parallel private(i, k)
!$omp do
	do k=1-g, kx+g
	do i=1-g, ix+g
		x(i, j  , k) = xc
		x(i, j-1, k) = xc
		x(i, j-2, k) = xc
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x4_d_f

subroutine bc_x5_d_f(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	k = sz(3)
!$omp parallel private(i, j)
!$omp do
	do j=1-g, jx+g
	do i=1-g, ix+g
		x(i, j, k+1) = xc
		x(i, j, k+2) = xc
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x5_d_f

subroutine bc_x6_d_f(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	k = 1
!$omp parallel private(i, j)
!$omp do
	do j=1-g, jx+g
	do i=1-g, ix+g
		x(i, j, k  ) = xc
		x(i, j, k-1) = xc
		x(i, j, k-2) = xc
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x6_d_f

subroutine bc_x1_p(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(j, k)
!$omp do
	do k=1-g, kx+g
	do j=1-g, jx+g
		x(ix+1, j, k) = x(1, j, k)
		x(ix+2, j, k) = x(2, j, k)
		x( 0, j, k) = x(ix  , j, k)
		x(-1, j, k) = x(ix-1, j, k)
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x1_p

subroutine bc_x2_p(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, k)
!$omp do
	do k=1-g, kx+g
	do i=1-g, ix+g
		x(i, jx+1, k) = x(i, 1, k)
		x(i, jx+2, k) = x(i, 2, k)
		x(i,  0, k) = x(i, jx  , k)
		x(i, -1, k) = x(i, jx-1, k)
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x2_p

subroutine bc_x5_p(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j)
!$omp do
	do j=1-g, jx+g
	do i=1-g, ix+g
		x(i, j, kx+1) = x(i, j, 1)
		x(i, j, kx+2) = x(i, j, 2)
		x(i, j,  0) = x(i, j, kx  )
		x(i, j, -1) = x(i, j, kx-1)
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x5_p

subroutine bc_x1_ls(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i = sz(1)
!$omp parallel private(j, k)
!$omp do
	do k=1, kx
	do j=1, jx
		x(i+1, j, k) = x(i, j, k) + x(i, j, k) - x(i-1, j, k)
		x(i+2, j, k) = x(i, j, k) 
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x1_ls

subroutine bc_x3_ls(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i = 1
!$omp parallel private(j, k)
!$omp do
	do k=1, kx
	do j=1, jx
		x(i-1, j, k) = x(i, j, k) + x(i, j, k) - x(i+1, j, k)
		x(i-2, j, k) = x(i, j, k)
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x3_ls

subroutine bc_x2_ls(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	j = sz(2)
!$omp parallel private(i, k)
!$omp do
	do k=1, kx
	do i=1, ix
		x(i, j+1, k) = x(i, j, k)*2.0 - x(i, j-1, k)
		x(i, j+2, k) = x(i, j+1, k)
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x2_ls

subroutine bc_x4_ls(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	j = 1
!$omp parallel private(i, k)
!$omp do
	do k=1, kx
	do i=1, ix
		x(i, j-1, k) = x(i, j, k)*2.0 - x(i, j+1, k)
		x(i, j-2, k) = x(i, j-1, k)
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x4_ls

subroutine bc_x5_ls(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	k = sz(3)
!$omp parallel private(i, j)
!$omp do
	do j=1, jx
	do i=1, ix
		x(i, j, k+1) = x(i, j, k) + x(i, j, k) - x(i, j, k-1)
		x(i, j, k+2) = x(i, j, k) 
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x5_ls

subroutine bc_x6_ls(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	real										:: pi, thetas
	real										:: dpx, dpy
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	k = 1
	pi = atan(1.0)*4.0
	thetas = 163.0/180.0*pi
	xc = -cos(thetas)
!$omp parallel private(i, j) &
!$omp          private(dpx, dpy)
!$omp do
	do j=1, jx
	do i=1, ix
		dpx = 0.5*(x(i+1, j, k) - x(i-1, j, k))
		dpy = 0.5*(x(i, j+1, k) - x(i, j-1, k))
		xc = -cos(thetas)/sin(thetas)*sqrt(dpx*dpx + dpy*dpy)
		x(i, j, k-1) = x(i, j, k) - xc
		x(i, j, k-2) = x(i, j, k) - xc*2.0

		xc = -cos(thetas)
		x(i, j, k-1) = x(i, j, k+1) - xc*2.0
		x(i, j, k-2) = x(i, j, k+1) - xc*3.0
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x6_ls

subroutine bc_x6_s(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	k = 1
!$omp parallel private(i, j)
!$omp do
	do j=1, jx
	do i=1, ix
!		if( (i-ix/2) + (j-jx/2) <  ix/8.0*sqrt(2.0) .and. &
!		    (i-ix/2) + (j-jx/2) > -ix/8.0*sqrt(2.0) .and. &
!			 -(i-ix/2) + (j-jx/2) <  ix/8.0*sqrt(2.0) .and. &
!			 -(i-ix/2) + (j-jx/2) > -ix/8.0*sqrt(2.0) ) then
		if( abs(i-ix/2.0) < ix/8.0 .and. abs(j-jx/2.0) < jx/8.0 ) then
!		if( abs(i-0.5-ix/2.0) < ix/20.0 ) then
  		x(i, j, k-1) = 2.0*xc - x(i, j, k)
  		x(i, j, k-2) = xc
		else
  		x(i, j, k-1) = x(i, j, k) - xc
  		x(i, j, k-2) = x(i, j, k) - xc*2.0
		end if
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine bc_x6_s

