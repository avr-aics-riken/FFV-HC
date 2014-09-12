!>  @file  bils4.f90
!!  @brief Basic subprograms for iterative linear solvers (BILS) for Cartesian grid data structure (for symmetric 7-band matrices) 
!<

subroutine rbgs_smoother_4_d(x, A, b, param, color, offset, sz, g, &
															mx, &
															sx1, sx3, sx2, sx4, sx5, sx6, &
															rx1, rx3, rx2, rx4, rx5, rx6, &
															node)
	implicit none
	include 'mpif.h'
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:3)	:: A
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
	real										:: param
	integer									:: color, offset
	real										:: r
	real										:: x_tmp
	integer									:: m
	integer									:: mx
	real, dimension(1:mx, 1-g:sz(2)+g, 1-g:sz(3)+g) :: sx1, sx3, rx1, rx3
	real, dimension(1:mx, 1-g:sz(1)+g, 1-g:sz(3)+g) :: sx2, sx4, rx2, rx4
	real, dimension(1:mx, 1-g:sz(1)+g, 1-g:sz(2)+g) :: sx5, sx6, rx5, rx6
	integer, dimension(0:6)	:: node
	integer									:: stats(mpi_status_size, 2)
	integer									:: req(12)
	integer									:: ierror
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
#ifdef __K_FAPP
call fapp_start("rbgs_smoother_4_d", 0, 0)
#endif
!$omp parallel private(i, j, k) &
!$omp					,private(r) &
!$omp					,private(color) &
!$omp					,private(x_tmp)
!$omp do schedule(static, 1)
	do k=0, kx+1
	do j=0, jx+1
	do i=0+mod(k+j+color+offset+1,2), ix+1, 2
   r = b(i, j, k) &
          - A(i  , j, k, 1)*x(i+1, j, k) &
				  - A(i-1, j, k, 1)*x(i-1, j, k) &
				  - A(i, j  , k, 2)*x(i, j+1, k) &
				  - A(i, j-1, k, 2)*x(i, j-1, k) &
				  - A(i, j, k  , 3)*x(i, j, k+1) &
				  - A(i, j, k-1, 3)*x(i, j, k-1)
		x_tmp = r/A(i, j, k, 0)
		x(i, j, k) = (1.0 - param)*x(i, j, k) + param*x_tmp
	end do
	end do
	end do
!$omp end do
	color = 1 - color
!$omp do schedule(static, 1)
	do k=1, kx
	do j=1, jx
!ocl norecurrence(A)
!ocl norecurrence(x)
	do i=1+mod(k+j+color+offset+1,2), ix, 2
   r = b(i, j, k) &
          - A(i  , j, k, 1)*x(i+1, j, k) &
				  - A(i-1, j, k, 1)*x(i-1, j, k) &
				  - A(i, j  , k, 2)*x(i, j+1, k) &
				  - A(i, j-1, k, 2)*x(i, j-1, k) &
				  - A(i, j, k  , 3)*x(i, j, k+1) &
				  - A(i, j, k-1, 3)*x(i, j, k-1)
		x_tmp = r/A(i, j, k, 0)
		x(i, j, k) = (1.0 - param)*x(i, j, k) + param*x_tmp
	end do
	end do
	end do
!$omp end do
!$omp end parallel
#ifdef __K_FAPP
call fapp_stop("rbgs_smoother_4_d", 0, 0)
#endif

#ifdef __K_FAPP
call fapp_start("rbgs_smoother_4_d", 1, 0)
#endif
!$omp parallel private(i, j, k) &
!$omp					,private(m)
!$omp do
	do k=1-g, kx+g
	do j=1-g, jx+g
	do m=1, mx
		sx1(m, j, k) = x(ix-mx+m, j, k)
		sx3(m, j, k) = x(m      , j, k)
	end do
	end do
	end do
!$omp end do
!$omp do
	do k=1-g, kx+g
	do i=1-g, ix+g
	do m=1, mx
		sx2(m, i, k) = x(i, jx-mx+m, k)
		sx4(m, i, k) = x(i, m      , k)
	end do
	end do
	end do
!$omp end do
!$omp do
	do j=1-g, jx+g
	do i=1-g, ix+g
	do m=1, mx
		sx5(m, i, j) = x(i, j, kx-mx+m)
		sx6(m, i, j) = x(i, j, m      )
	end do
	end do
	end do
!$omp end do
!$omp end parallel
#ifdef __K_FAPP
call fapp_stop("rbgs_smoother_4_d", 1, 0)
#endif

#ifdef __K_FAPP
call fapp_start("rbgs_smoother_4_d", 2, 0)
#endif

#ifdef _REAL_IS_DOUBLE_
		call mpi_sendrecv(sx1, (jx+2*g)*(kx+2*g)*mx, mpi_double_precision, node(1), 1, &
											rx3, (jx+2*g)*(kx+2*g)*mx, mpi_double_precision, node(3), 1, &
											mpi_comm_world, mpi_statuses_ignore, ierror)

		call mpi_sendrecv(sx3, (jx+2*g)*(kx+2*g)*mx, mpi_double_precision, node(3), 3, &
											rx1, (jx+2*g)*(kx+2*g)*mx, mpi_double_precision, node(1), 3, &
											mpi_comm_world, mpi_statuses_ignore, ierror)

		call mpi_sendrecv(sx2, (ix+2*g)*(kx+2*g)*mx, mpi_double_precision, node(2), 2, &
											rx4, (ix+2*g)*(kx+2*g)*mx, mpi_double_precision, node(4), 2, &
											mpi_comm_world, mpi_statuses_ignore, ierror)

		call mpi_sendrecv(sx4, (ix+2*g)*(kx+2*g)*mx, mpi_double_precision, node(4), 4, &
											rx2, (ix+2*g)*(kx+2*g)*mx, mpi_double_precision, node(2), 4, &
											mpi_comm_world, mpi_statuses_ignore, ierror)

		call mpi_sendrecv(sx5, (ix+2*g)*(jx+2*g)*mx, mpi_double_precision, node(5), 5, &
											rx6, (ix+2*g)*(jx+2*g)*mx, mpi_double_precision, node(6), 5, &
											mpi_comm_world, mpi_statuses_ignore, ierror)

		call mpi_sendrecv(sx6, (ix+2*g)*(jx+2*g)*mx, mpi_double_precision, node(6), 6, &
											rx5, (ix+2*g)*(jx+2*g)*mx, mpi_double_precision, node(5), 6, &
											mpi_comm_world, mpi_statuses_ignore, ierror)
#else
!		call mpi_sendrecv(sx1, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(1), 1, &
!											rx3, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(3), 1, &
!											mpi_comm_world, mpi_statuses_ignore, ierror)
!		call mpi_sendrecv(sx3, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(3), 3, &
!											rx1, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(1), 3, &
!											mpi_comm_world, mpi_statuses_ignore, ierror)
!		call mpi_sendrecv(sx2, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(2), 2, &
!											rx4, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(4), 2, &
!											mpi_comm_world, mpi_statuses_ignore, ierror)
!		call mpi_sendrecv(sx4, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(4), 4, &
!											rx2, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(2), 4, &
!											mpi_comm_world, mpi_statuses_ignore, ierror)
!		call mpi_sendrecv(sx5, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(5), 5, &
!											rx6, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(6), 5, &
!											mpi_comm_world, mpi_statuses_ignore, ierror)
!		call mpi_sendrecv(sx6, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(6), 6, &
!											rx5, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(5), 6, &
!											mpi_comm_world, mpi_statuses_ignore, ierror)

		call mpi_irecv(rx1, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(1), 3, &
									 mpi_comm_world, req(4), ierror)
		call mpi_irecv(rx3, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(3), 1, &
									 mpi_comm_world, req(2), ierror)
		call mpi_isend(sx1, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(1), 1, &
									 mpi_comm_world, req(1), ierror)
		call mpi_isend(sx3, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(3), 3, &
									 mpi_comm_world, req(3), ierror)
		call mpi_waitall(4, req, mpi_statuses_ignore, ierror)

		call mpi_irecv(rx2, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(2), 4, &
									 mpi_comm_world, req(4), ierror)
		call mpi_irecv(rx4, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(4), 2, &
									 mpi_comm_world, req(2), ierror)
		call mpi_isend(sx2, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(2), 2, &
									 mpi_comm_world, req(1), ierror)
		call mpi_isend(sx4, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(4), 4, &
									 mpi_comm_world, req(3), ierror)
		call mpi_waitall(4, req, mpi_statuses_ignore, ierror)

		call mpi_irecv(rx5, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(5), 6, &
									 mpi_comm_world, req(4), ierror)
		call mpi_irecv(rx6, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(6), 5, &
									 mpi_comm_world, req(2), ierror)
		call mpi_isend(sx5, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(5), 5, &
									 mpi_comm_world, req(1), ierror)
		call mpi_isend(sx6, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(6), 6, &
									 mpi_comm_world, req(3), ierror)
		call mpi_waitall(4, req, mpi_statuses_ignore, ierror)
#endif

#ifdef __K_FAPP
call fapp_stop("rbgs_smoother_4_d", 2, 0)
#endif

#ifdef __K_FAPP
call fapp_start("rbgs_smoother_4_d", 3, 0)
#endif
!$omp parallel private(i, j, k) &
!$omp					,private(m)
!$omp do
	do k=1-g, kx+g
	do j=1-g, jx+g
	do m=1, mx
		x( ix+m, j, k) = rx1(m, j, k)
		x(-mx+m, j, k) = rx3(m, j, k)
	end do
	end do
	end do
!$omp end do
!$omp do
	do k=1-g, kx+g
	do i=1-g, ix+g
	do m=1, mx
		x(i,  jx+m, k) = rx2(m, i, k)
		x(i, -mx+m, k) = rx4(m, i, k)
	end do
	end do
	end do
!$omp end do
!$omp do
	do j=1-g, jx+g
	do i=1-g, ix+g
	do m=1, mx
		x(i, j,  kx+m) = rx5(m, i, j)
		x(i, j, -mx+m) = rx6(m, i, j)
	end do
	end do
	end do
!$omp end do
!$omp end parallel
#ifdef __K_FAPP
call fapp_stop("rbgs_smoother_4_d", 3, 0)
#endif
end subroutine rbgs_smoother_4_d

subroutine rbgs_smoother_4_n(x, A, b, param, color, offset, sz, g, &
															mx, &
															sx1, sx3, sx2, sx4, sx5, sx6, &
															rx1, rx3, rx2, rx4, rx5, rx6, &
															node)
	implicit none
	include 'mpif.h'
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:3)	:: A
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
	real										:: param
	integer									:: color, offset
	real										:: r
	real										:: x_tmp
	integer									:: m
	integer									:: mx
	real, dimension(1:mx, 1-g:sz(2)+g, 1-g:sz(3)+g) :: sx1, sx3, rx1, rx3
	real, dimension(1:mx, 1-g:sz(1)+g, 1-g:sz(3)+g) :: sx2, sx4, rx2, rx4
	real, dimension(1:mx, 1-g:sz(1)+g, 1-g:sz(2)+g) :: sx5, sx6, rx5, rx6
	integer, dimension(0:6)	:: node
	integer									:: stats(mpi_status_size, 2)
	integer									:: req(12)
	integer									:: ierror
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
#ifdef __K_FAPP
call fapp_start("rbgs_smoother_4_n", 0, 0)
#endif
!$omp parallel private(i, j, k) &
!$omp					,private(r) &
!$omp					,private(x_tmp)
!$omp do schedule(static, 1)
	do k=1, kx
	do j=1, jx
	do i=1+mod(k+j+color+offset+1,2), ix, 2
   r = b(i, j, k) &
          - A(i  , j, k, 1)*x(i+1, j, k) &
				  - A(i-1, j, k, 1)*x(i-1, j, k) &
				  - A(i, j  , k, 2)*x(i, j+1, k) &
				  - A(i, j-1, k, 2)*x(i, j-1, k) &
				  - A(i, j, k  , 3)*x(i, j, k+1) &
				  - A(i, j, k-1, 3)*x(i, j, k-1)
		x_tmp = r/A(i, j, k, 0)
		x(i, j, k) = (1.0 - param)*x(i, j, k) + param*x_tmp
	end do
	end do
	end do
!$omp end do
!$omp end parallel
#ifdef __K_FAPP
call fapp_stop("rbgs_smoother_4_n", 0, 0)
#endif

#ifdef __K_FAPP
call fapp_start("rbgs_smoother_4_n", 1, 0)
#endif
!$omp parallel private(i, j, k) &
!$omp					,private(m)
!$omp do
	do k=1-g, kx+g
	do j=1-g, jx+g
	do m=1, mx
		sx1(m, j, k) = x(ix-mx+m, j, k)
		sx3(m, j, k) = x(m      , j, k)
	end do
	end do
	end do
!$omp end do
!$omp do
	do k=1-g, kx+g
	do i=1-g, ix+g
	do m=1, mx
		sx2(m, i, k) = x(i, jx-mx+m, k)
		sx4(m, i, k) = x(i, m      , k)
	end do
	end do
	end do
!$omp end do
!$omp do
	do j=1-g, jx+g
	do i=1-g, ix+g
	do m=1, mx
		sx5(m, i, j) = x(i, j, kx-mx+m)
		sx6(m, i, j) = x(i, j, m      )
	end do
	end do
	end do
!$omp end do
!$omp end parallel
#ifdef __K_FAPP
call fapp_stop("rbgs_smoother_4_n", 1, 0)
#endif

#ifdef __K_FAPP
call fapp_start("rbgs_smoother_4_n", 2, 0)
#endif

#ifdef _REAL_IS_DOUBLE_
		call mpi_sendrecv(sx1, (jx+2*g)*(kx+2*g)*mx, mpi_double_precision, node(1), 1, &
											rx3, (jx+2*g)*(kx+2*g)*mx, mpi_double_precision, node(3), 1, &
											mpi_comm_world, mpi_statuses_ignore, ierror)

		call mpi_sendrecv(sx3, (jx+2*g)*(kx+2*g)*mx, mpi_double_precision, node(3), 3, &
											rx1, (jx+2*g)*(kx+2*g)*mx, mpi_double_precision, node(1), 3, &
											mpi_comm_world, mpi_statuses_ignore, ierror)

		call mpi_sendrecv(sx2, (ix+2*g)*(kx+2*g)*mx, mpi_double_precision, node(2), 2, &
											rx4, (ix+2*g)*(kx+2*g)*mx, mpi_double_precision, node(4), 2, &
											mpi_comm_world, mpi_statuses_ignore, ierror)

		call mpi_sendrecv(sx4, (ix+2*g)*(kx+2*g)*mx, mpi_double_precision, node(4), 4, &
											rx2, (ix+2*g)*(kx+2*g)*mx, mpi_double_precision, node(2), 4, &
											mpi_comm_world, mpi_statuses_ignore, ierror)

		call mpi_sendrecv(sx5, (ix+2*g)*(jx+2*g)*mx, mpi_double_precision, node(5), 5, &
											rx6, (ix+2*g)*(jx+2*g)*mx, mpi_double_precision, node(6), 5, &
											mpi_comm_world, mpi_statuses_ignore, ierror)

		call mpi_sendrecv(sx6, (ix+2*g)*(jx+2*g)*mx, mpi_double_precision, node(6), 6, &
											rx5, (ix+2*g)*(jx+2*g)*mx, mpi_double_precision, node(5), 6, &
											mpi_comm_world, mpi_statuses_ignore, ierror)
#else
!		call mpi_sendrecv(sx1, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(1), 1, &
!											rx3, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(3), 1, &
!											mpi_comm_world, mpi_statuses_ignore, ierror)
!		call mpi_sendrecv(sx3, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(3), 3, &
!											rx1, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(1), 3, &
!											mpi_comm_world, mpi_statuses_ignore, ierror)
!		call mpi_sendrecv(sx2, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(2), 2, &
!											rx4, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(4), 2, &
!											mpi_comm_world, mpi_statuses_ignore, ierror)
!		call mpi_sendrecv(sx4, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(4), 4, &
!											rx2, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(2), 4, &
!											mpi_comm_world, mpi_statuses_ignore, ierror)
!		call mpi_sendrecv(sx5, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(5), 5, &
!											rx6, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(6), 5, &
!											mpi_comm_world, mpi_statuses_ignore, ierror)
!		call mpi_sendrecv(sx6, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(6), 6, &
!											rx5, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(5), 6, &
!											mpi_comm_world, mpi_statuses_ignore, ierror)

		call mpi_irecv(rx1, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(1), 3, &
									 mpi_comm_world, req(4), ierror)
		call mpi_irecv(rx3, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(3), 1, &
									 mpi_comm_world, req(2), ierror)
		call mpi_isend(sx1, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(1), 1, &
									 mpi_comm_world, req(1), ierror)
		call mpi_isend(sx3, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(3), 3, &
									 mpi_comm_world, req(3), ierror)
		call mpi_waitall(4, req, mpi_statuses_ignore, ierror)

		call mpi_irecv(rx2, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(2), 4, &
									 mpi_comm_world, req(4), ierror)
		call mpi_irecv(rx4, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(4), 2, &
									 mpi_comm_world, req(2), ierror)
		call mpi_isend(sx2, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(2), 2, &
									 mpi_comm_world, req(1), ierror)
		call mpi_isend(sx4, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(4), 4, &
									 mpi_comm_world, req(3), ierror)
		call mpi_waitall(4, req, mpi_statuses_ignore, ierror)

		call mpi_irecv(rx5, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(5), 6, &
									 mpi_comm_world, req(4), ierror)
		call mpi_irecv(rx6, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(6), 5, &
									 mpi_comm_world, req(2), ierror)
		call mpi_isend(sx5, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(5), 5, &
									 mpi_comm_world, req(1), ierror)
		call mpi_isend(sx6, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(6), 6, &
									 mpi_comm_world, req(3), ierror)
		call mpi_waitall(4, req, mpi_statuses_ignore, ierror)
#endif

#ifdef __K_FAPP
call fapp_stop("rbgs_smoother_4_n", 2, 0)
#endif

#ifdef __K_FAPP
call fapp_start("rbgs_smoother_4_n", 3, 0)
#endif
!$omp parallel private(i, j, k) &
!$omp					,private(m)
!$omp do
	do k=1-g, kx+g
	do j=1-g, jx+g
	do m=1, mx
		x( ix+m, j, k) = rx1(m, j, k)
		x(-mx+m, j, k) = rx3(m, j, k)
	end do
	end do
	end do
!$omp end do
!$omp do
	do k=1-g, kx+g
	do i=1-g, ix+g
	do m=1, mx
		x(i,  jx+m, k) = rx2(m, i, k)
		x(i, -mx+m, k) = rx4(m, i, k)
	end do
	end do
	end do
!$omp end do
!$omp do
	do j=1-g, jx+g
	do i=1-g, ix+g
	do m=1, mx
		x(i, j,  kx+m) = rx5(m, i, j)
		x(i, j, -mx+m) = rx6(m, i, j)
	end do
	end do
	end do
!$omp end do
!$omp end parallel
#ifdef __K_FAPP
call fapp_stop("rbgs_smoother_4_n", 3, 0)
#endif
end subroutine rbgs_smoother_4_n

subroutine rbgs_smoother_4_p(x, A, b, param, color, offset, sz, g, &
															mx, &
															sx1, sx3, sx2, sx4, sx5, sx6, &
															rx1, rx3, rx2, rx4, rx5, rx6, &
															node)
	implicit none
	include 'mpif.h'
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:3)	:: A
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
	real										:: param
	integer									:: color, offset
	real										:: r
	real										:: x_tmp
	integer									:: m
	integer									:: mx
	real, dimension(1:mx, 1-g:sz(2)+g, 1-g:sz(3)+g) :: sx1, sx3, rx1, rx3
	real, dimension(1:mx, 1-g:sz(1)+g, 1-g:sz(3)+g) :: sx2, sx4, rx2, rx4
	real, dimension(1:mx, 1-g:sz(1)+g, 1-g:sz(2)+g) :: sx5, sx6, rx5, rx6
	integer, dimension(0:6)	:: node
	integer									:: stats(mpi_status_size, 2)
	integer									:: req(12)
	integer									:: ierror
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
#ifdef __K_FAPP
call fapp_start("rbgs_smoother_4_p", 0, 0)
#endif
!$omp parallel private(i, j, k) &
!$omp					,private(r) &
!$omp					,private(x_tmp)
	i=ix
!$omp do schedule(static, 1) 
	do k=1, kx
	do j=1+mod(k-1+ix-1+color+offset,2), jx, 2
   r = b(i, j, k) &
          - A(i  , j, k, 1)*x(i+1, j, k) &
				  - A(i-1, j, k, 1)*x(i-1, j, k) &
				  - A(i, j  , k, 2)*x(i, j+1, k) &
				  - A(i, j-1, k, 2)*x(i, j-1, k) &
				  - A(i, j, k  , 3)*x(i, j, k+1) &
				  - A(i, j, k-1, 3)*x(i, j, k-1)
		x_tmp = r/A(i, j, k, 0)
		x(i, j, k) = (1.0 - param)*x(i, j, k) + param*x_tmp
	end do
	end do
!$omp end do
	i=1
!$omp do schedule(static, 1)
	do k=1, kx
	do j=1+mod(k-1+color+offset,2), jx, 2
   r = b(i, j, k) &
          - A(i  , j, k, 1)*x(i+1, j, k) &
				  - A(i-1, j, k, 1)*x(i-1, j, k) &
				  - A(i, j  , k, 2)*x(i, j+1, k) &
				  - A(i, j-1, k, 2)*x(i, j-1, k) &
				  - A(i, j, k  , 3)*x(i, j, k+1) &
				  - A(i, j, k-1, 3)*x(i, j, k-1)
		x_tmp = r/A(i, j, k, 0)
		x(i, j, k) = (1.0 - param)*x(i, j, k) + param*x_tmp
	end do
	end do
!$omp end do
	j=jx
!$omp do schedule(static, 1)
	do k=1, kx
	do i=1+mod(k-1+jx-1+color+offset,2), ix, 2
   r = b(i, j, k) &
          - A(i  , j, k, 1)*x(i+1, j, k) &
				  - A(i-1, j, k, 1)*x(i-1, j, k) &
				  - A(i, j  , k, 2)*x(i, j+1, k) &
				  - A(i, j-1, k, 2)*x(i, j-1, k) &
				  - A(i, j, k  , 3)*x(i, j, k+1) &
				  - A(i, j, k-1, 3)*x(i, j, k-1)
		x_tmp = r/A(i, j, k, 0)
		x(i, j, k) = (1.0 - param)*x(i, j, k) + param*x_tmp
	end do
	end do
!$omp end do
	j=1
!$omp do schedule(static, 1)
	do k=1, kx
	do i=1+mod(k-1+color+offset,2), ix, 2
   r = b(i, j, k) &
          - A(i  , j, k, 1)*x(i+1, j, k) &
				  - A(i-1, j, k, 1)*x(i-1, j, k) &
				  - A(i, j  , k, 2)*x(i, j+1, k) &
				  - A(i, j-1, k, 2)*x(i, j-1, k) &
				  - A(i, j, k  , 3)*x(i, j, k+1) &
				  - A(i, j, k-1, 3)*x(i, j, k-1)
		x_tmp = r/A(i, j, k, 0)
		x(i, j, k) = (1.0 - param)*x(i, j, k) + param*x_tmp
	end do
	end do
!$omp end do
	k=kx
!$omp do schedule(static, 1)
	do j=1, jx
	do i=1+mod(j-1+kx-1+color+offset,2), ix, 2
   r = b(i, j, k) &
          - A(i  , j, k, 1)*x(i+1, j, k) &
				  - A(i-1, j, k, 1)*x(i-1, j, k) &
				  - A(i, j  , k, 2)*x(i, j+1, k) &
				  - A(i, j-1, k, 2)*x(i, j-1, k) &
				  - A(i, j, k  , 3)*x(i, j, k+1) &
				  - A(i, j, k-1, 3)*x(i, j, k-1)
		x_tmp = r/A(i, j, k, 0)
		x(i, j, k) = (1.0 - param)*x(i, j, k) + param*x_tmp
	end do
	end do
!$omp end do
	k=1
!$omp do schedule(static, 1) 
	do j=1, jx
	do i=1+mod(j-1+color+offset,2), ix, 2
   r = b(i, j, k) &
          - A(i  , j, k, 1)*x(i+1, j, k) &
				  - A(i-1, j, k, 1)*x(i-1, j, k) &
				  - A(i, j  , k, 2)*x(i, j+1, k) &
				  - A(i, j-1, k, 2)*x(i, j-1, k) &
				  - A(i, j, k  , 3)*x(i, j, k+1) &
				  - A(i, j, k-1, 3)*x(i, j, k-1)
		x_tmp = r/A(i, j, k, 0)
		x(i, j, k) = (1.0 - param)*x(i, j, k) + param*x_tmp
	end do
	end do
!$omp end do
!$omp end parallel
#ifdef __K_FAPP
call fapp_stop("rbgs_smoother_4_p", 0, 0)
#endif

#ifdef __K_FAPP
call fapp_start("rbgs_smoother_4_p", 1, 0)
#endif
!$omp parallel private(i, j, k) &
!$omp					,private(m)
!$omp do
	do k=1-g, kx+g
	do j=1-g, jx+g
	do m=1, mx
		sx1(m, j, k) = x(ix-mx+m, j, k)
		sx3(m, j, k) = x(m      , j, k)
	end do
	end do
	end do
!$omp end do
!$omp do
	do k=1-g, kx+g
	do i=1-g, ix+g
	do m=1, mx
		sx2(m, i, k) = x(i, jx-mx+m, k)
		sx4(m, i, k) = x(i, m      , k)
	end do
	end do
	end do
!$omp end do
!$omp do
	do j=1-g, jx+g
	do i=1-g, ix+g
	do m=1, mx
		sx5(m, i, j) = x(i, j, kx-mx+m)
		sx6(m, i, j) = x(i, j, m      )
	end do
	end do
	end do
!$omp end do
!$omp end parallel
#ifdef __K_FAPP
call fapp_stop("rbgs_smoother_4_p", 1, 0)
#endif

#ifdef __K_FAPP
call fapp_start("rbgs_smoother_4_p", 2, 0)
#endif
!$omp parallel private(i, j, k) &
!$omp					,private(r) &
!$omp					,private(x_tmp) &
!$omp					,private(ierror)

!$omp master
#ifdef _REAL_IS_DOUBLE_
		call mpi_sendrecv(sx1, (jx+2*g)*(kx+2*g)*mx, mpi_double_precision, node(1), 1, &
											rx3, (jx+2*g)*(kx+2*g)*mx, mpi_double_precision, node(3), 1, &
											mpi_comm_world, mpi_statuses_ignore, ierror)

		call mpi_sendrecv(sx3, (jx+2*g)*(kx+2*g)*mx, mpi_double_precision, node(3), 3, &
											rx1, (jx+2*g)*(kx+2*g)*mx, mpi_double_precision, node(1), 3, &
											mpi_comm_world, mpi_statuses_ignore, ierror)

		call mpi_sendrecv(sx2, (ix+2*g)*(kx+2*g)*mx, mpi_double_precision, node(2), 2, &
											rx4, (ix+2*g)*(kx+2*g)*mx, mpi_double_precision, node(4), 2, &
											mpi_comm_world, mpi_statuses_ignore, ierror)

		call mpi_sendrecv(sx4, (ix+2*g)*(kx+2*g)*mx, mpi_double_precision, node(4), 4, &
											rx2, (ix+2*g)*(kx+2*g)*mx, mpi_double_precision, node(2), 4, &
											mpi_comm_world, mpi_statuses_ignore, ierror)

		call mpi_sendrecv(sx5, (ix+2*g)*(jx+2*g)*mx, mpi_double_precision, node(5), 5, &
											rx6, (ix+2*g)*(jx+2*g)*mx, mpi_double_precision, node(6), 5, &
											mpi_comm_world, mpi_statuses_ignore, ierror)

		call mpi_sendrecv(sx6, (ix+2*g)*(jx+2*g)*mx, mpi_double_precision, node(6), 6, &
											rx5, (ix+2*g)*(jx+2*g)*mx, mpi_double_precision, node(5), 6, &
											mpi_comm_world, mpi_statuses_ignore, ierror)
#else
!		call mpi_sendrecv(sx1, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(1), 1, &
!											rx3, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(3), 1, &
!											mpi_comm_world, mpi_statuses_ignore, ierror)
!		call mpi_sendrecv(sx3, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(3), 3, &
!											rx1, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(1), 3, &
!											mpi_comm_world, mpi_statuses_ignore, ierror)
!		call mpi_sendrecv(sx2, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(2), 2, &
!											rx4, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(4), 2, &
!											mpi_comm_world, mpi_statuses_ignore, ierror)
!		call mpi_sendrecv(sx4, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(4), 4, &
!											rx2, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(2), 4, &
!											mpi_comm_world, mpi_statuses_ignore, ierror)
!		call mpi_sendrecv(sx5, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(5), 5, &
!											rx6, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(6), 5, &
!											mpi_comm_world, mpi_statuses_ignore, ierror)
!		call mpi_sendrecv(sx6, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(6), 6, &
!											rx5, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(5), 6, &
!											mpi_comm_world, mpi_statuses_ignore, ierror)

		call mpi_irecv(rx1, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(1), 3, &
									 mpi_comm_world, req(4), ierror)
		call mpi_irecv(rx3, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(3), 1, &
									 mpi_comm_world, req(2), ierror)
		call mpi_isend(sx1, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(1), 1, &
									 mpi_comm_world, req(1), ierror)
		call mpi_isend(sx3, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(3), 3, &
									 mpi_comm_world, req(3), ierror)
		call mpi_waitall(4, req, mpi_statuses_ignore, ierror)

		call mpi_irecv(rx2, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(2), 4, &
									 mpi_comm_world, req(4), ierror)
		call mpi_irecv(rx4, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(4), 2, &
									 mpi_comm_world, req(2), ierror)
		call mpi_isend(sx2, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(2), 2, &
									 mpi_comm_world, req(1), ierror)
		call mpi_isend(sx4, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(4), 4, &
									 mpi_comm_world, req(3), ierror)
		call mpi_waitall(4, req, mpi_statuses_ignore, ierror)

		call mpi_irecv(rx5, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(5), 6, &
									 mpi_comm_world, req(4), ierror)
		call mpi_irecv(rx6, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(6), 5, &
									 mpi_comm_world, req(2), ierror)
		call mpi_isend(sx5, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(5), 5, &
									 mpi_comm_world, req(1), ierror)
		call mpi_isend(sx6, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(6), 6, &
									 mpi_comm_world, req(3), ierror)
		call mpi_waitall(4, req, mpi_statuses_ignore, ierror)
#endif
!$omp end master

!$omp do schedule(static, 1)
	do k=2, kx-1
	do j=2, jx-1
	do i=2+mod(k+j+color+offset+1,2), ix-1, 2
   r = b(i, j, k) &
          - A(i  , j, k, 1)*x(i+1, j, k) &
				  - A(i-1, j, k, 1)*x(i-1, j, k) &
				  - A(i, j  , k, 2)*x(i, j+1, k) &
				  - A(i, j-1, k, 2)*x(i, j-1, k) &
				  - A(i, j, k  , 3)*x(i, j, k+1) &
				  - A(i, j, k-1, 3)*x(i, j, k-1)
		x_tmp = r/A(i, j, k, 0)
		x(i, j, k) = (1.0 - param)*x(i, j, k) + param*x_tmp
	end do
	end do
	end do
!$omp end do
!$omp end parallel
#ifdef __K_FAPP
call fapp_stop("rbgs_smoother_4_p", 2, 0)
#endif

#ifdef __K_FAPP
call fapp_start("rbgs_smoother_4_p", 3, 0)
#endif
!$omp parallel private(i, j, k) &
!$omp					,private(m)
!$omp do
	do k=1-g, kx+g
	do j=1-g, jx+g
	do m=1, mx
		x( ix+m, j, k) = rx1(m, j, k)
		x(-mx+m, j, k) = rx3(m, j, k)
	end do
	end do
	end do
!$omp end do
!$omp do
	do k=1-g, kx+g
	do i=1-g, ix+g
	do m=1, mx
		x(i,  jx+m, k) = rx2(m, i, k)
		x(i, -mx+m, k) = rx4(m, i, k)
	end do
	end do
	end do
!$omp end do
!$omp do
	do j=1-g, jx+g
	do i=1-g, ix+g
	do m=1, mx
		x(i, j,  kx+m) = rx5(m, i, j)
		x(i, j, -mx+m) = rx6(m, i, j)
	end do
	end do
	end do
!$omp end do
!$omp end parallel
#ifdef __K_FAPP
call fapp_stop("rbgs_smoother_4_p", 3, 0)
#endif
end subroutine rbgs_smoother_4_p

subroutine rbgs_smoother_4_b(x, A, b, param, color, offset, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
  integer                 :: ii, jj, kk
  integer                 :: it, jt, kt
  integer, dimension(3)   :: szt
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:3)	:: A
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
	real										:: param
	integer									:: color, offset
	real										:: r
	real										:: x_tmp
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
  it = 1
  jt = 16
	kt = 1
#ifdef __K_FAPP
call fapp_start("rbgs_smoother_4_b", 0, 0)
#endif
!$omp parallel private(i, j, k) &
!$omp					,private(r, x_tmp)	
!$omp do schedule(static, 1)
  do k=1, kx
  do j=1, jx
  do i=1+mod(k+j+color+offset,2), ix, 2
   r = b(i, j, k) &
          - A(i  , j, k, 1)*x(i+1, j, k) &
				  - A(i-1, j, k, 1)*x(i-1, j, k) &
				  - A(i, j  , k, 2)*x(i, j+1, k) &
				  - A(i, j-1, k, 2)*x(i, j-1, k) &
				  - A(i, j, k  , 3)*x(i, j, k+1) &
				  - A(i, j, k-1, 3)*x(i, j, k-1)
		x_tmp = r/A(i, j, k, 0)
		x(i, j, k) = (1.0 - param)*x(i, j, k) + param*x_tmp
	end do
	end do
	end do
!$omp end do
!$omp end parallel
#ifdef __K_FAPP
call fapp_stop("rbgs_smoother_4_b", 0, 0)
#endif
end subroutine rbgs_smoother_4_b

subroutine calc_ax_4_b(y, A, x, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: y 
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:3)	:: A
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
				+A(i  , j, k, 1)*x(i+1, j, k) &
				+A(i-1, j, k, 1)*x(i-1, j, k) &
				+A(i, j  , k, 2)*x(i, j+1, k) &
				+A(i, j-1, k, 2)*x(i, j-1, k) &
				+A(i, j, k  , 3)*x(i, j, k+1) &
				+A(i, j, k-1, 3)*x(i, j, k-1)
		y(i, j, k) = ax
	end do
	end do
	end do
!$omp end do
!$omp end parallel
! 13 flop
end subroutine calc_ax_4_b

subroutine calc_r_4_b(r, A, b, x, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: r
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:3)	:: A
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
		    -A(i, j, k, 0)*x(i, j, k) &
				-A(i  , j, k, 1)*x(i+1, j, k) &
				-A(i-1, j, k, 1)*x(i-1, j, k) &
				-A(i, j  , k, 2)*x(i, j+1, k) &
				-A(i, j-1, k, 2)*x(i, j-1, k) &
				-A(i, j, k  , 3)*x(i, j, k+1) &
				-A(i, j, k-1, 3)*x(i, j, k-1)
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine calc_r_4_b

subroutine calc_rr_4_b(rr, A, b, x, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:3)	:: A
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
        -A(i, j, k, 0)*x(i, j, k) &
				-A(i  , j, k, 1)*x(i+1, j, k) &
				-A(i-1, j, k, 1)*x(i-1, j, k) &
				-A(i, j  , k, 2)*x(i, j+1, k) &
				-A(i, j-1, k, 2)*x(i, j-1, k) &
				-A(i, j, k  , 3)*x(i, j, k+1) &
				-A(i, j, k-1, 3)*x(i, j, k-1)
		rr = rr + r*r
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine calc_rr_4_b


subroutine rbgs_smoother_4_s(x, A0, A1, A2, A3, b, param, color, offset, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
  integer                 :: ii, jj, kk
  integer                 :: it, jt, kt
  integer, dimension(3)   :: szt
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: A0, A1, A2, A3
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
	real										:: param
	integer									:: color, offset
	real										:: b_ax
	real										:: x_tmp
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
  it = 1
  jt = 8
	kt = 1
!$omp parallel private(i, j, k) &
!$omp					,private(ii, jj, kk) &
!$omp					,private(b_ax) &
!$omp					,private(x_tmp)
!$omp do 
  do jj=1, jx, jt
  do k=1, kx
	do j=jj, min(jj + jt - 1, jx)
	do i=1+mod(k+j+color+offset,2), ix, 2
    b_ax = b(i, j, k) &
          - A1(i  , j, k)*x(i+1, j, k) &
          - A1(i-1, j, k)*x(i-1, j, k) &
          - A2(i, j  , k)*x(i, j+1, k) &
          - A2(i, j-1, k)*x(i, j-1, k) &
          - A3(i, j, k  )*x(i, j, k+1) &
          - A3(i, j, k-1)*x(i, j, k-1)
		x_tmp = b_ax/A0(i, j, k)
		x(i, j, k) = (1.0 - param)*x(i, j, k) + param*x_tmp
	end do
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine rbgs_smoother_4_s

subroutine calc_ax_4_s(y, A0, A1, A2, A3, x, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: y 
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: A0, A1, A2, A3
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
				+A1(i  , j, k)*x(i+1, j, k) &
				+A1(i-1, j, k)*x(i-1, j, k) &
				+A2(i, j  , k)*x(i, j+1, k) &
				+A2(i, j-1, k)*x(i, j-1, k) &
				+A3(i, j, k  )*x(i, j, k+1) &
				+A3(i, j, k-1)*x(i, j, k-1)
		y(i, j, k) = ax
	end do
	end do
	end do
!$omp end do
!$omp end parallel
! 13 flop
end subroutine calc_ax_4_s

subroutine calc_r_4_s(r, A0, A1, A2, A3, b, x, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: r
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: A0, A1, A2, A3
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
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
				+A1(i  , j, k)*x(i+1, j, k) &
				+A1(i-1, j, k)*x(i-1, j, k) &
				+A2(i, j  , k)*x(i, j+1, k) &
				+A2(i, j-1, k)*x(i, j-1, k) &
				+A3(i, j, k  )*x(i, j, k+1) &
				+A3(i, j, k-1)*x(i, j, k-1)
		r(i, j, k) = b(i, j, k) - ax
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine calc_r_4_s

subroutine calc_rr_4_s(rr, A0, A1, A2, A3, b, x, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: A0, A1, A2, A3
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: ax
	real										:: rr
	real										:: r
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	rr = 0.0
!$omp parallel private(i, j, k) &
!$omp					,private(ax) &
!$omp					,private(r)
!$omp do schedule(static, 1), reduction(+:rr)
	do k=1, kx
	do j=1, jx
	do i=1, ix
		ax = A0(i, j, k)*x(i, j, k) &
				+A1(i  , j, k)*x(i+1, j, k) &
				+A1(i-1, j, k)*x(i-1, j, k) &
				+A2(i, j  , k)*x(i, j+1, k) &
				+A2(i, j-1, k)*x(i, j-1, k) &
				+A3(i, j, k  )*x(i, j, k+1) &
				+A3(i, j, k-1)*x(i, j, k-1)
		r = b(i, j, k) - ax
		rr = rr + r*r
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine calc_rr_4_s


