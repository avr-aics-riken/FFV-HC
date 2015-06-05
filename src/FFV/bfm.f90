
subroutine bfm_fan( &
                fx, fy, fz, &
                ux0, uy0, uz0, &
                rid, &
                rid_target, &
								b, &
								nx, ny, nz, &
								dpmax, umax, &
                dx, dt, &
                sz, g)
  implicit none
  integer                  :: i, j, k
  integer                  :: ix, jx, kx
  integer                  :: g
  integer, dimension(3)    :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: fx, fy, fz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: ux0, uy0, uz0
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: rid
	integer									:: rid_target
	real										:: b
	real										:: nx, ny, nz
	real										:: dpmax, umax
  real                    :: dx, dt
	real										:: ux, uy, uz, un
	real										:: dp
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp          private(ux, uy, uz, un) &
!$omp          private(dp)
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
		if( rid(i, j, k) == rid_target ) then
			ux = ux0(i ,j, k)
			uy = uy0(i ,j, k)
			uz = uz0(i ,j, k)
			un = ux*nx + uy*ny + uz*nz
			dp = dpmax*(1.0 - (un/umax)**2)
			dp = dpmax*(1.0 - (un/umax))
			if( un < 0.0 ) then
				dp = dpmax
			else if( un > umax ) then
				dp = 0.0
			end if
      write(*,*) dp
			fx(i, j, k) = (dp/b)*nx
			fy(i, j, k) = (dp/b)*ny
			fz(i, j, k) = (dp/b)*nz
		end if
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bfm_fan

subroutine bfm_hex( &
                fx, fy, fz, &
                ux0, uy0, uz0, &
                rid, &
                rid_target, &
								b, &
								nx, ny, nz, &
								dpmax, umax, &
                dx, dt, &
                sz, g)
  implicit none
  integer                  :: i, j, k
  integer                  :: ix, jx, kx
  integer                  :: g
  integer, dimension(3)    :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: fx, fy, fz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: ux0, uy0, uz0
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: rid
	integer									:: rid_target
	real										:: b
	real										:: nx, ny, nz
	real										:: dpmax, umax
  real                    :: dx, dt
	real										:: ux, uy, uz, un, u2, ul
	real										:: dp, dp_
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp          private(ux, uy, uz, un, u2, ul) &
!$omp          private(dp, dp_)
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
		if( rid(i, j, k) == rid_target ) then
			ux = ux0(i ,j, k)
			uy = uy0(i ,j, k)
			uz = uz0(i ,j, k)
			u2 = ux*ux + uy*uy + uz*uz
			ul = sqrt(u2)
			dp = dpmax*(ul/umax)**2
			dp_= dpmax*(ul/umax**2)
      write(*,*) dp
			fx(i, j, k) = -(dp_/b)*ux
			fy(i, j, k) = -(dp_/b)*uy
			fz(i, j, k) = -(dp_/b)*uz
		end if
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bfm_hex

