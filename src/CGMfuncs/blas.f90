!>  @file  blas.f90
!!  @brief Basic subprograms for linear algebra for Cartesian grid data structure (for symmetric/asymmetric 7-band matrices)
!<

subroutine fill(x, a, sz, g)
  implicit none
  integer, dimension(3)   :: sz
  integer                 :: g
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: x
  real                    :: a 
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k)
!$omp do
#else
#endif
  do k=1-g, kx+g
  do j=1-g, jx+g
  do i=1-g, ix+g
    x(i, j, k) = a
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine fill

subroutine fill_vf3d(x, a, sz, g, ne)
  implicit none
  integer, dimension(3)   :: sz
  integer                 :: g
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: l
  integer                 :: ne
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:ne-1)  :: x
  real                    :: a
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(l, i, j, k)
!$omp do
#else
#endif
  do l=0, ne-1
  do k=1-g, kx+g
  do j=1-g, jx+g
  do i=1-g, ix+g
    x(i, j, k, l) = a
  end do
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine fill_vf3d

subroutine copy(y, x, sz, g)
  implicit none
  integer, dimension(3)    :: sz
  integer                  :: g
  integer                  :: i, j, k
  integer                  :: ix, jx, kx
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: y
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: x
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k)
!$omp do
#else
#endif
  do k=1-g, kx+g
  do j=1-g, jx+g
  do i=1-g, ix+g
    y(i, j, k) = x(i, j, k)
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine copy

subroutine copy_integer(y, x, sz, g)
  implicit none
  integer, dimension(3)    :: sz
  integer                  :: g
  integer                  :: i, j, k
  integer                  :: ix, jx, kx
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: y
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: x
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k)
!$omp do
#else
#endif
  do k=1-g, kx+g
  do j=1-g, jx+g
  do i=1-g, ix+g
    y(i, j, k) = x(i, j, k)
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine copy_integer

subroutine add(C, A, B, sz, g)
  implicit none
  integer, dimension(3)    :: sz
  integer                  :: g
  integer                  :: i, j, k
  integer                  :: ix, jx, kx
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: C
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: A, B
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k)
!$omp do
#else
!ocl nouxsimd
!ocl serial
!dir$ noparallel
#endif
  do k=1, kx
  do j=1, jx
  do i=1, ix
    C(i, j, k) = A(i, j, k) + B(i, j, k)
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine add

subroutine adda(y, a, sz, g)
  implicit none
  integer, dimension(3)    :: sz
  integer                  :: g
  integer                  :: i, j, k
  integer                  :: ix, jx, kx
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: y
  real                    :: a
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k)
!$omp do
#else
!ocl nouxsimd
!ocl serial
!dir$ noparallel
#endif
  do k=1, kx
  do j=1, jx
  do i=1, ix
    y(i, j, k) = y(i, j, k) + a
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine adda

subroutine ave(C, A, B, sz, g)
  implicit none
  integer, dimension(3)    :: sz
  integer                  :: g
  integer                  :: i, j, k
  integer                  :: ix, jx, kx
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: C
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: A, B
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k)
!$omp do
#else
!ocl nouxsimd
!ocl serial
!dir$ noparallel
#endif
  do k=1, kx
  do j=1, jx
  do i=1, ix
    C(i, j, k) = 0.5*(A(i, j, k) + B(i, j, k))
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine ave

subroutine triad(C, A, B, d, sz, g)
  implicit none
  integer, dimension(3)    :: sz
  integer                  :: g
  integer                  :: i, j, k
  integer                  :: ix, jx, kx
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: C
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: A, B
  real                    :: d
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k)
!$omp do
#else
!ocl nouxsimd
!ocl serial
!dir$ noparallel
#endif
  do k=1, kx
  do j=1, jx
  do i=1, ix
    C(i, j, k) = A(i, j, k) + d*B(i, j, k)
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine triad

subroutine scal(y, a, sz, g)
  implicit none
  integer, dimension(3)    :: sz
  integer                  :: g
  integer                  :: i, j, k
  integer                  :: ix, jx, kx
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: y
  real                    :: a
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k)
!$omp do
#else
!ocl nouxsimd
!ocl serial
!dir$ noparallel
#endif
  do k=1, kx
  do j=1, jx
  do i=1, ix
    y(i, j, k) = a*y(i, j, k)
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine scal

subroutine avew(y, x, a, sz, g)
  implicit none
  integer, dimension(3)    :: sz
  integer                  :: g
  integer                  :: i, j, k
  integer                  :: ix, jx, kx
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: y
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: x
  real                    :: a
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k)
!$omp do
#else
!ocl nouxsimd
!ocl serial
!dir$ noparallel
#endif
  do k=1, kx
  do j=1, jx
  do i=1, ix
    y(i, j, k) = (1.0 - a)*y(i, j, k) + a*x(i, j, k)
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine avew

subroutine avew_2(y, x0, x1, a, sz, g)
  implicit none
  integer, dimension(3)    :: sz
  integer                  :: g
  integer                  :: i, j, k
  integer                  :: ix, jx, kx
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: y
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: x0, x1
  real                    :: a
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k)
!$omp do
#else
!ocl nouxsimd
!ocl serial
!dir$ noparallel
#endif
  do k=1, kx
  do j=1, jx
  do i=1, ix
    y(i, j, k) = (1.0 - a)*y(i, j, k) + a*x0(i, j, k)*x1(i, j, k)
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine avew_2

subroutine axpy(y, x, a, sz, g)
  implicit none
  integer, dimension(3)    :: sz
  integer                  :: g
  integer                  :: i, j, k
  integer                  :: ix, jx, kx
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: y
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: x
  real                    :: a
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k)
!$omp do
#else
!ocl nouxsimd
!ocl serial
!dir$ noparallel
#endif
  do k=1, kx
  do j=1, jx
  do i=1, ix
    y(i, j, k) = y(i, j, k) + a*x(i, j, k)
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine axpy

subroutine xpay(y, x, a, sz, g)
  implicit none
  integer, dimension(3)    :: sz
  integer                  :: g
  integer                  :: i, j, k
  integer                  :: ix, jx, kx
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: y
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: x
  real                    :: a
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k)
!$omp do
#else
!ocl nouxsimd
!ocl serial
!dir$ noparallel
#endif
  do k=1, kx
  do j=1, jx
  do i=1, ix
    y(i, j, k) = x(i, j, k) + a*y(i, j, k)
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine xpay

subroutine axpyz(z, x, y, a, sz, g)
  implicit none
  integer, dimension(3)    :: sz
  integer                  :: g
  integer                  :: i, j, k
  integer                  :: ix, jx, kx
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: z
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: x
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: y
  real                    :: a
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k)
!$omp do
#else
!ocl nouxsimd
!ocl serial
!dir$ noparallel
#endif
  do k=1, kx
  do j=1, jx
  do i=1, ix
    z(i, j, k) = a*x(i, j, k) + y(i, j, k)
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine axpyz

subroutine axpbypz(z, x, y, a, b, sz, g)
  implicit none
  integer                  :: i, j, k
  integer                  :: ix, jx, kx
  integer                  :: g
  integer, dimension(3)    :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: z
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: x, y
  real                    :: a, b
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k)
!$omp do
#else
!ocl nouxsimd
!ocl serial
!dir$ noparallel
#endif
  do k=1, kx
  do j=1, jx
  do i=1, ix
    z(i, j, k) = a*x(i, j, k) + b*y(i, j, k) + z(i, j, k)
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine axpbypz

subroutine dot(xy, x, y, sz, g)
  implicit none
  integer, dimension(3)    :: sz
  integer                  :: g
  integer                  :: i, j, k
  integer                  :: ix, jx, kx
  real                    :: xy
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: x
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: y
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  xy = 0.0
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k)
!$omp do reduction(+:xy)
#else
!ocl nouxsimd
!ocl serial
!dir$ noparallel
#endif
  do k=1, kx
  do j=1, jx
  do i=1, ix
    xy = xy + x(i, j, k)*y(i, j, k)
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine dot

subroutine dotx2(xy, xz, x, y, z, sz, g)
  implicit none
  integer, dimension(3)    :: sz
  integer                  :: g
  integer                  :: i, j, k
  integer                  :: ix, jx, kx
  real                    :: xy, xz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: x
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: y
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: z
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  xy = 0.0
  xz = 0.0
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k)
!$omp do reduction(+:xy,xz)
#else
!ocl nouxsimd
!ocl serial
!dir$ noparallel
#endif
  do k=1, kx
  do j=1, jx
  do i=1, ix
    xy = xy + x(i, j, k)*y(i, j, k)
    xz = xz + x(i, j, k)*z(i, j, k)
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine dotx2

