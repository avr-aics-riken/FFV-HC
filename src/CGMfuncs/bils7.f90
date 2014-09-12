!>  @file  bils7.f90
!!  @brief Basic subprograms for iterative linear solvers (BILS) for Cartesian grid data structure (for symmetric/asymmetric 7-band matrices) 
!<

subroutine rbgs_smoother_7_b(x, A, b, param, color, offset, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:6)	:: A
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
	real										:: param
	integer									:: color, offset
	real										:: r
	real										:: x_tmp
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(r, x_tmp)	
!$omp do schedule(static, 1)
	do k=1, kx
	do j=1, jx
	do i=1+mod(k+j+color+offset,2), ix, 2
   r = b(i, j, k) &
          - A(i, j, k, 1)*x(i+1, j, k) &
				  - A(i, j, k, 3)*x(i-1, j, k) &
				  - A(i, j, k, 2)*x(i, j+1, k) &
				  - A(i, j, k, 4)*x(i, j-1, k) &
				  - A(i, j, k, 5)*x(i, j, k+1) &
				  - A(i, j, k, 6)*x(i, j, k-1)
		x_tmp = r/A(i, j, k, 0)
		x(i, j, k) = x(i, j, k) - param*(x(i, j, k) - x_tmp)
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine rbgs_smoother_7_b

subroutine calc_ax_7_b(y, A, x, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: y 
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:6)	:: A
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: ax
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(ax)	
!$omp do schedule(static, 1)
	do k=1, kx
	do j=1, jx
	do i=1, ix
		ax = A(i, j, k, 0)*x(i, j, k) &
				+A(i, j, k, 1)*x(i+1, j, k) &
				+A(i, j, k, 3)*x(i-1, j, k) &
				+A(i, j, k, 2)*x(i, j+1, k) &
				+A(i, j, k, 4)*x(i, j-1, k) &
				+A(i, j, k, 5)*x(i, j, k+1) &
				+A(i, j, k, 6)*x(i, j, k-1)
		y(i, j, k) = ax
	end do
	end do
	end do
!$omp end do
!$omp end parallel
! 13 flop
end subroutine calc_ax_7_b

subroutine calc_r_7_b(r, A, b, x, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: r
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:6)	:: A
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j, k) 
!$omp do schedule(static, 1)
	do k=1, kx
	do j=1, jx
	do i=1, ix
		r(i, j, k) = b(i, j, k) &
								- A(i, j, k, 0)*x(i, j, k) &
								- A(i, j, k, 1)*x(i+1, j, k) &
								- A(i, j, k, 3)*x(i-1, j, k) &
								- A(i, j, k, 2)*x(i, j+1, k) &
								- A(i, j, k, 4)*x(i, j-1, k) &
								- A(i, j, k, 5)*x(i, j, k+1) &
								- A(i, j, k, 6)*x(i, j, k-1)
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine calc_r_7_b

subroutine calc_rr_7_b(rr, A, b, x, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:6)	:: A
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: rr
	real										:: r
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	rr = 0.0
!$omp parallel private(i, j, k) &
!$omp					,private(r)
!$omp do schedule(static, 1), reduction(+:rr)
	do k=1, kx
	do j=1, jx
	do i=1, ix
		r = b(i, j, k) &
			- A(i, j, k, 0)*x(i, j, k) &
			- A(i, j, k, 1)*x(i+1, j, k) &
			- A(i, j, k, 3)*x(i-1, j, k) &
			- A(i, j, k, 2)*x(i, j+1, k) &
			- A(i, j, k, 4)*x(i, j-1, k) &
			- A(i, j, k, 5)*x(i, j, k+1) &
			- A(i, j, k, 6)*x(i, j, k-1)
		rr = rr + r*r
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine calc_rr_7_b


subroutine rbgs_smoother_7_s(x, A0, A1, A2, A3, A4, A5, A6, b, param, color, offset, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
  integer                 :: ii, jj, kk
  integer                 :: it, jt, kt
  integer, dimension(3)   :: szt
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: A0, A1, A2, A3, A4, A5, A6
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
	real										:: param
	integer									:: color, offset
	real										:: ax
	real										:: x_old, x_new, x_tmp
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
  it = 1
  jt = 16
	kt = 1
!$omp parallel private(i, j, k) &
!$omp					,private(ii, jj, kk) &
!$omp					,private(ax) &
!$omp					,private(x_old, x_tmp, x_new)
!$omp do schedule(static, 1)
  do jj=1, jx, jt
  do k=1, kx
	do j=jj, min(jj + jt - 1, jx)
	do i=1+mod(k+j+color+offset,2), ix, 2
		ax = A1(i, j, k)*x(i+1, j, k) &
				+A3(i, j, k)*x(i-1, j, k) &
				+A2(i, j, k)*x(i, j+1, k) &
				+A4(i, j, k)*x(i, j-1, k) &
				+A5(i, j, k)*x(i, j, k+1) &
				+A6(i, j, k)*x(i, j, k-1)
		x_old = x(i, j, k)
		x_tmp = (b(i, j, k) - ax)/A0(i, j, k)
		x_new = x_old - param*(x_old - x_tmp)
		x(i, j, k) = x_new
	end do
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine rbgs_smoother_7_s

subroutine calc_ax_7_s(y, A0, A1, A2, A3, A4, A5, A6, x, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: y 
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: A0, A1, A2, A3, A4, A5, A6
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: ax
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j, k) &
!$omp					,private(ax)	
!$omp do schedule(static, 1)
	do k=1, kx
	do j=1, jx
	do i=1, ix
		ax = A0(i, j, k)*x(i, j, k) &
				+A1(i, j, k)*x(i+1, j, k) &
				+A3(i, j, k)*x(i-1, j, k) &
				+A2(i, j, k)*x(i, j+1, k) &
				+A4(i, j, k)*x(i, j-1, k) &
				+A5(i, j, k)*x(i, j, k+1) &
				+A6(i, j, k)*x(i, j, k-1)
		y(i, j, k) = ax
	end do
	end do
	end do
!$omp end do
!$omp end parallel
! 13 flop
end subroutine calc_ax_7_s

subroutine calc_r_7_s(r, A0, A1, A2, A3, A4, A5, A6, b, x, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: r
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: A0, A1, A2, A3, A4, A5, A6
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j, k) 
!$omp do schedule(static, 1)
	do k=1, kx
	do j=1, jx
	do i=1, ix
		r(i, j, k) &
			= b(i, j, k) &
			- A0(i, j, k)*x(i, j, k) &
			- A1(i, j, k)*x(i+1, j, k) &
			- A3(i, j, k)*x(i-1, j, k) &
			- A2(i, j, k)*x(i, j+1, k) &
			- A4(i, j, k)*x(i, j-1, k) &
			- A5(i, j, k)*x(i, j, k+1) &
			- A6(i, j, k)*x(i, j, k-1)
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine calc_r_7_s

subroutine calc_rr_7_s(rr, A0, A1, A2, A3, A4, A5, A6, b, x, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: A0, A1, A2, A3, A4, A5, A6
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: rr
	real										:: r
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	rr = 0.0
!$omp parallel private(i, j, k) &
!$omp					,private(r)
!$omp do schedule(static, 1), reduction(+:rr)
	do k=1, kx
	do j=1, jx
	do i=1, ix
		r = b(i, j, k) &
			- A0(i, j, k)*x(i, j, k) &
			- A1(i, j, k)*x(i+1, j, k) &
			- A3(i, j, k)*x(i-1, j, k) &
			- A2(i, j, k)*x(i, j+1, k) &
			- A4(i, j, k)*x(i, j-1, k) &
			- A5(i, j, k)*x(i, j, k+1) &
			- A6(i, j, k)*x(i, j, k-1)
		rr = rr + r*r
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine calc_rr_7_s


