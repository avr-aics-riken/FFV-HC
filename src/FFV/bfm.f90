
subroutine bfm_fan( &
                fx, fy, fz, &
                ux0, uy0, uz0, &
                rid, &
                rid_target, &
								b, &
								nx, ny, nz, &
								c0, c1 ,c2, &
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
	real										:: c0, c1, c2
  real                    :: dx, dt
	real										:: ux, uy, uz, un, u2, ul
	real										:: dp
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp          private(ux, uy, uz, un, u2, ul) &
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
			u2 = ux*ux + uy*uy + uz*uz
			ul = sqrt(u2)
			dp = c0 + c1*ul + c2*ul*ul
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
								c0, c1 ,c2, &
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
	real										:: c0, c1, c2
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
			dp = c0 + c1*ul + c2*ul*ul
			dp_= c1 + c2*ul
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

