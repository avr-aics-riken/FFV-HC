function bcut_getminmod(a, b)
  implicit none
  real          :: bcut_getminmod
  real          :: a, b
!  bcut_getminmod = 0.5*(sign(a, a) + sign(b, b))*min(abs(a), abs(b))
  bcut_getminmod = 0.0d0
  if( abs(a) .lt. abs(b) ) then
    bcut_getminmod = a
  else if( abs(a) .gt. abs(b) ) then
    bcut_getminmod = b
  endif
  return
end function bcut_getminmod

function bcut_getupwind(u, n, p)
  implicit none
  real          :: bcut_getupwind
  real          :: u
  real          :: n, p
  bcut_getupwind = 0.0d0
  if( u .gt. 0 ) then
    bcut_getupwind = n
  else if( u .lt. 0 ) then
    bcut_getupwind = p
  endif
  return
end function bcut_getupwind

function bcut_getweno3(v1, v2, v3)
  implicit none
  real          :: bcut_getweno3
  real          :: v1, v2, v3
  real          :: r
  real          :: w
  real          :: eps
  eps = 0.001
  r = (eps + (v2-v1)**2)/(eps + (v3-v2)**2)
  w = 1.0/(1.0 + 2*r*r)
  bcut_getweno3 = 0.5*(v2 + v3) - 0.5*w*(v1 - 2.0*v2 + v3)
  return
end function bcut_getweno3

subroutine bcut_calc_c_f_u1( &
                fc, &
                f, &
                vw, ve, vs, vn, vb, vt, &
                c0, c1, c2, c3, c4, c5, &
                cid0, cid1, cid2, cid3, cid4, cid5, &
                pid, &
                dx, dt, &
                fi, &
                sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: fc
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: f
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: vw, ve, vs, vn, vb, vt
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
  integer                  :: cidp
  integer                  :: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
  integer                  :: cidw, cide, cids, cidn, cidb, cidt
  integer                  :: cidw0, cide1, cids2, cidn3, cidb4, cidt5
  integer                  :: pidp, pidw, pide, pids, pidn, pidb, pidt
  real                    :: dx, dt
  real                    :: fi
  real                    :: d0, d1, d2, d3, d4, d5
  real                    :: m0, m1, m2, m3, m4, m5
  real                    :: l0, l1, l2, l3, l4, l5
  real                    :: vx0, vx1, vy2, vy3, vz4, vz5
  real                    :: fp, fw, fe, fs, fn, fb, ft
  real                    :: fww, fee, fss, fnn, fbb, ftt
  real                    :: f0, f1, f2, f3, f4, f5
  real                    :: dfx_p, dfx_c, dfx_n
  real                    :: dfy_p, dfy_c, dfy_n
  real                    :: dfz_p, dfz_c, dfz_n
  real                    :: bcut_getupwind
  real                    :: bcut_getminmod
  real                    :: q0, q1, q2, q3, q4, q5
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp           private(cidp) &
!$omp           private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp           private(cidw, cide, cids, cidn, cidb, cidt) &
!$omp           private(cidw0, cide1, cids2, cidn3, cidb4, cidt5) &
!$omp           private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp           private(d0, d1, d2, d3, d4, d5) &
!$omp           private(m0, m1, m2, m3, m4, m5) &
!$omp           private(l0, l1, l2, l3, l4, l5) &
!$omp           private(vx0, vx1, vy2, vy3, vz4, vz5) &
!$omp           private(fp, fw, fe, fs, fn, fb, ft) &
!$omp           private(fww, fee, fss, fnn, fbb, ftt) &
!$omp           private(f0, f1, f2, f3, f4, f5) &
!$omp           private(dfx_p, dfx_c, dfx_n) &
!$omp           private(dfy_p, dfy_c, dfy_n) &
!$omp           private(dfz_p, dfz_c, dfz_n) &
!$omp           private(q0, q1, q2, q3, q4, q5)
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
    vx0 = vw(i, j, k)
    vx1 = ve(i, j, k)
    vy2 = vs(i, j, k)
    vy3 = vn(i, j, k)
    vz4 = vb(i, j, k)
    vz5 = vt(i, j, k)

    fp = f(i, j, k)
    fw = f(i-1, j, k)
    fe = f(i+1, j, k)
    fs = f(i, j-1, k)
    fn = f(i, j+1, k)
    fb = f(i, j, k-1)
    ft = f(i, j, k+1)

    d0 = c0(i, j, k)
    d1 = c1(i, j, k)
    d2 = c2(i, j, k)
    d3 = c3(i, j, k)
    d4 = c4(i, j, k)
    d5 = c5(i, j, k)

    cidp0 = cid0(i, j, k)
    cidp1 = cid1(i, j, k)
    cidp2 = cid2(i, j, k)
    cidp3 = cid3(i, j, k)
    cidp4 = cid4(i, j, k)
    cidp5 = cid5(i, j, k)

    pidp = pid(i, j, k)

    f0 = bcut_getupwind(vx0, fw, fp)
    f1 = bcut_getupwind(vx1, fp, fe)
    f2 = bcut_getupwind(vy2, fs, fp)
    f3 = bcut_getupwind(vy3, fp, fn)
    f4 = bcut_getupwind(vz4, fb, fp)
    f5 = bcut_getupwind(vz5, fp, ft)

    q0 = vx0*f0
    q1 = vx1*f1
    q2 = vy2*f2
    q3 = vy3*f3
    q4 = vz4*f4
    q5 = vz5*f5

    m0 = 0.0d0
    m1 = 0.0d0
    m2 = 0.0d0
    m3 = 0.0d0
    m4 = 0.0d0
    m5 = 0.0d0

    if( cidp0 /= 0 ) then
      q0 = (d0 - 0.5d0)/(d0 + 0.5d0)*q1
      m0 = 1.0d0
    endif
    if( cidp1 /= 0 ) then
      q1 = (d1 - 0.5d0)/(d1 + 0.5d0)*q0
      m1 = 1.0d0
    endif
    if( cidp2 /= 0 ) then
      q2 = (d2 - 0.5d0)/(d2 + 0.5d0)*q3
      m2 = 1.0d0
    endif
    if( cidp3 /= 0 ) then
      q3 = (d3 - 0.5d0)/(d3 + 0.5d0)*q2
      m3 = 1.0d0
    endif
    if( cidp4 /= 0 ) then
      q4 = (d4 - 0.5d0)/(d4 + 0.5d0)*q5
      m4 = 1.0d0
    endif
    if( cidp5 /= 0 ) then
      q5 = (d5 - 0.5d0)/(d5 + 0.5d0)*q4
      m5 = 1.0d0
    endif

    if( cidp0 /= 0 .and. cidp1 /= 0 ) then
      q0 = 0.0d0
      q1 = 0.0d0
    endif
    if( cidp2 /= 0 .and. cidp3 /= 0 ) then
      q2 = 0.0d0
      q3 = 0.0d0
    endif
    if( cidp4 /= 0 .and. cidp5 /= 0 ) then
      q4 = 0.0d0
      q5 = 0.0d0
    endif

    fc(i, j, k) = (q1 - q0)/dx &
                + (q3 - q2)/dx &
                + (q5 - q4)/dx

    if( pidp /= 1 ) then
      fc(i, j, k) = 0.0d0
    endif

  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_calc_c_f_u1

subroutine bcut_calc_c_f_c2( &
                fc, &
                f, &
                vw, ve, vs, vn, vb, vt, &
                c0, c1, c2, c3, c4, c5, &
                cid0, cid1, cid2, cid3, cid4, cid5, &
                pid, &
                dx, dt, &
                fi, &
                sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: fc
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: f
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: vw, ve, vs, vn, vb, vt
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
  integer                  :: cidp
  integer                  :: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
  integer                  :: cidw, cide, cids, cidn, cidb, cidt
  integer                  :: cidw0, cide1, cids2, cidn3, cidb4, cidt5
  integer                  :: pidp, pidw, pide, pids, pidn, pidb, pidt
  real                    :: dx, dt
  real                    :: fi
  real                    :: d0, d1, d2, d3, d4, d5
  real                    :: m0, m1, m2, m3, m4, m5
  real                    :: l0, l1, l2, l3, l4, l5
  real                    :: vx0, vx1, vy2, vy3, vz4, vz5
  real                    :: fp, fw, fe, fs, fn, fb, ft
  real                    :: fww, fee, fss, fnn, fbb, ftt
  real                    :: f0, f1, f2, f3, f4, f5
  real                    :: dfx_p, dfx_c, dfx_n
  real                    :: dfy_p, dfy_c, dfy_n
  real                    :: dfz_p, dfz_c, dfz_n
  real                    :: bcut_getupwind
  real                    :: bcut_getminmod
  real                    :: q0, q1, q2, q3, q4, q5
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp           private(cidp) &
!$omp           private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp           private(cidw, cide, cids, cidn, cidb, cidt) &
!$omp           private(cidw0, cide1, cids2, cidn3, cidb4, cidt5) &
!$omp           private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp           private(d0, d1, d2, d3, d4, d5) &
!$omp           private(m0, m1, m2, m3, m4, m5) &
!$omp           private(l0, l1, l2, l3, l4, l5) &
!$omp           private(vx0, vx1, vy2, vy3, vz4, vz5) &
!$omp           private(fp, fw, fe, fs, fn, fb, ft) &
!$omp           private(fww, fee, fss, fnn, fbb, ftt) &
!$omp           private(f0, f1, f2, f3, f4, f5) &
!$omp           private(dfx_p, dfx_c, dfx_n) &
!$omp           private(dfy_p, dfy_c, dfy_n) &
!$omp           private(dfz_p, dfz_c, dfz_n) &
!$omp           private(q0, q1, q2, q3, q4, q5)
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
    vx0 = vw(i, j, k)
    vx1 = ve(i, j, k)
    vy2 = vs(i, j, k)
    vy3 = vn(i, j, k)
    vz4 = vb(i, j, k)
    vz5 = vt(i, j, k)

    fp = f(i, j, k)
    fw = f(i-1, j, k)
    fe = f(i+1, j, k)
    fs = f(i, j-1, k)
    fn = f(i, j+1, k)
    fb = f(i, j, k-1)
    ft = f(i, j, k+1)

    d0 = c0(i, j, k)
    d1 = c1(i, j, k)
    d2 = c2(i, j, k)
    d3 = c3(i, j, k)
    d4 = c4(i, j, k)
    d5 = c5(i, j, k)

    cidp0 = cid0(i, j, k)
    cidp1 = cid1(i, j, k)
    cidp2 = cid2(i, j, k)
    cidp3 = cid3(i, j, k)
    cidp4 = cid4(i, j, k)
    cidp5 = cid5(i, j, k)

    pidp = pid(i, j, k)

    f0 = 0.5d0*(fp + fw)
    f1 = 0.5d0*(fp + fe)
    f2 = 0.5d0*(fp + fs)
    f3 = 0.5d0*(fp + fn)
    f4 = 0.5d0*(fp + fb)
    f5 = 0.5d0*(fp + ft)

    q0 = vx0*f0
    q1 = vx1*f1
    q2 = vy2*f2
    q3 = vy3*f3
    q4 = vz4*f4
    q5 = vz5*f5

    m0 = 0.0d0
    m1 = 0.0d0
    m2 = 0.0d0
    m3 = 0.0d0
    m4 = 0.0d0
    m5 = 0.0d0

    if( cidp0 /= 0 ) then
      q0 = (d0 - 0.5d0)/(d0 + 0.5d0)*q1
      m0 = 1.0d0
    endif
    if( cidp1 /= 0 ) then
      q1 = (d1 - 0.5d0)/(d1 + 0.5d0)*q0
      m1 = 1.0d0
    endif
    if( cidp2 /= 0 ) then
      q2 = (d2 - 0.5d0)/(d2 + 0.5d0)*q3
      m2 = 1.0d0
    endif
    if( cidp3 /= 0 ) then
      q3 = (d3 - 0.5d0)/(d3 + 0.5d0)*q2
      m3 = 1.0d0
    endif
    if( cidp4 /= 0 ) then
      q4 = (d4 - 0.5d0)/(d4 + 0.5d0)*q5
      m4 = 1.0d0
    endif
    if( cidp5 /= 0 ) then
      q5 = (d5 - 0.5d0)/(d5 + 0.5d0)*q4
      m5 = 1.0d0
    endif

    if( cidp0 /= 0 .and. cidp1 /= 0 ) then
      q0 = 0.0d0
      q1 = 0.0d0
    endif
    if( cidp2 /= 0 .and. cidp3 /= 0 ) then
      q2 = 0.0d0
      q3 = 0.0d0
    endif
    if( cidp4 /= 0 .and. cidp5 /= 0 ) then
      q4 = 0.0d0
      q5 = 0.0d0
    endif

    fc(i, j, k) = (q1 - q0)/dx &
                + (q3 - q2)/dx &
                + (q5 - q4)/dx

    if( pidp /= 1 ) then
      fc(i, j, k) = 0.0d0
    endif

  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_calc_c_f_c2

subroutine bcut_calc_c_f_e3( &
                fc, &
                f, &
                vw, ve, vs, vn, vb, vt, &
                c0, c1, c2, c3, c4, c5, &
                cid0, cid1, cid2, cid3, cid4, cid5, &
                pid, &
                dx, dt, &
                fi, &
                sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: fc
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: f
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: vw, ve, vs, vn, vb, vt
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
  integer                  :: cidp
  integer                  :: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
  integer                  :: cidw, cide, cids, cidn, cidb, cidt
  integer                  :: cidw0, cide1, cids2, cidn3, cidb4, cidt5
  integer                  :: pidp, pidw, pide, pids, pidn, pidb, pidt
  real                    :: dx, dt
  real                    :: fi
  real                    :: d0, d1, d2, d3, d4, d5
  real                    :: m0, m1, m2, m3, m4, m5
  real                    :: l0, l1, l2, l3, l4, l5
  real                    :: vx0, vx1, vy2, vy3, vz4, vz5
  real                    :: fp, fw, fe, fs, fn, fb, ft
  real                    :: fww, fee, fss, fnn, fbb, ftt
  real                    :: f0, f1, f2, f3, f4, f5
  real                    :: dfx_p, dfx_c, dfx_n
  real                    :: dfy_p, dfy_c, dfy_n
  real                    :: dfz_p, dfz_c, dfz_n
  real                    :: bcut_getupwind
  real                    :: bcut_getminmod
  real                    :: q0, q1, q2, q3, q4, q5
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp           private(cidp) &
!$omp           private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp           private(cidw, cide, cids, cidn, cidb, cidt) &
!$omp           private(cidw0, cide1, cids2, cidn3, cidb4, cidt5) &
!$omp           private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp           private(d0, d1, d2, d3, d4, d5) &
!$omp           private(m0, m1, m2, m3, m4, m5) &
!$omp           private(l0, l1, l2, l3, l4, l5) &
!$omp           private(vx0, vx1, vy2, vy3, vz4, vz5) &
!$omp           private(fp, fw, fe, fs, fn, fb, ft) &
!$omp           private(fww, fee, fss, fnn, fbb, ftt) &
!$omp           private(f0, f1, f2, f3, f4, f5) &
!$omp           private(dfx_p, dfx_c, dfx_n) &
!$omp           private(dfy_p, dfy_c, dfy_n) &
!$omp           private(dfz_p, dfz_c, dfz_n) &
!$omp           private(q0, q1, q2, q3, q4, q5)
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
    vx0 = vw(i, j, k)
    vx1 = ve(i, j, k)
    vy2 = vs(i, j, k)
    vy3 = vn(i, j, k)
    vz4 = vb(i, j, k)
    vz5 = vt(i, j, k)

    fp = f(i, j, k)
    fw = f(i-1, j, k)
    fe = f(i+1, j, k)
    fs = f(i, j-1, k)
    fn = f(i, j+1, k)
    fb = f(i, j, k-1)
    ft = f(i, j, k+1)
    fww = f(i-2, j, k)
    fee = f(i+2, j, k)
    fss = f(i, j-2, k)
    fnn = f(i, j+2, k)
    fbb = f(i, j, k-2)
    ftt = f(i, j, k+2)

    d0 = c0(i, j, k)
    d1 = c1(i, j, k)
    d2 = c2(i, j, k)
    d3 = c3(i, j, k)
    d4 = c4(i, j, k)
    d5 = c5(i, j, k)

    cidp0 = cid0(i, j, k)
    cidp1 = cid1(i, j, k)
    cidp2 = cid2(i, j, k)
    cidp3 = cid3(i, j, k)
    cidp4 = cid4(i, j, k)
    cidp5 = cid5(i, j, k)

    cidw0 = cid0(i-1, j, k)
    cide1 = cid1(i+1, j, k)
    cids2 = cid2(i, j-1, k)
    cidn3 = cid3(i, j+1, k)
    cidb4 = cid4(i, j, k-1)
    cidt5 = cid5(i, j, k+1)

    pidp = pid(i, j, k)

    if( cidp0 /= 0 ) then
      fw  = (1.0d0 - 1.0d0/d0)*fp + (1.0d0/d0)*fi
      fww = (1.0d0 - 2.0d0/d0)*fp + (2.0d0/d0)*fi
      m0 = 1.0d0
    endif
    if( cidp1 /= 0 ) then
      fe  = (1.0d0 - 1.0d0/d1)*fp + (1.0d0/d1)*fi
      fee = (1.0d0 - 2.0d0/d1)*fp + (2.0d0/d1)*fi
      m1 = 1.0d0
    endif
    if( cidp2 /= 0 ) then
      fs  = (1.0d0 - 1.0d0/d2)*fp + (1.0d0/d2)*fi
      fss = (1.0d0 - 2.0d0/d2)*fp + (2.0d0/d2)*fi
      m2 = 1.0d0
    endif
    if( cidp3 /= 0 ) then
      fn  = (1.0d0 - 1.0d0/d3)*fp + (1.0d0/d3)*fi
      fnn = (1.0d0 - 2.0d0/d3)*fp + (2.0d0/d3)*fi
      m3 = 1.0d0
    endif
    if( cidp4 /= 0 ) then
      fb  = (1.0d0 - 1.0d0/d4)*fp + (1.0d0/d4)*fi
      fbb = (1.0d0 - 2.0d0/d4)*fp + (2.0d0/d4)*fi
      m4 = 1.0d0
    endif
    if( cidp5 /= 0 ) then
      ft  = (1.0d0 - 1.0d0/d5)*fp + (1.0d0/d5)*fi
      ftt = (1.0d0 - 2.0d0/d5)*fp + (2.0d0/d5)*fi
      m5 = 1.0d0
    endif

    if( cidw0 /= 0 ) then
      fww = fi
    endif
    if( cide1 /= 0 ) then
      fee = fi
    endif
    if( cids2 /= 0 ) then
      fss = fi
    endif
    if( cidn3 /= 0 ) then
      fnn = fi
    endif
    if( cidb4 /= 0 ) then
      fbb = fi
    endif
    if( cidt5 /= 0 ) then
      ftt = fi
    endif

    dfx_p = bcut_getminmod(fee-fe, fe-fp )
    dfx_c = bcut_getminmod(fe -fp, fp-fw )
    dfx_n = bcut_getminmod(fp -fw, fw-fww)
    dfy_p = bcut_getminmod(fnn-fn, fn-fp )
    dfy_c = bcut_getminmod(fn -fp, fp-fs )
    dfy_n = bcut_getminmod(fp -fs, fs-fss)
    dfz_p = bcut_getminmod(ftt-ft, ft-fp )
    dfz_c = bcut_getminmod(ft -fp, fp-fb )
    dfz_n = bcut_getminmod(fp -fb, fb-fbb)

    f0 = bcut_getupwind(vx0, fw + 0.5*dfx_n, fp - 0.5*dfx_c)
    f1 = bcut_getupwind(vx1, fp + 0.5*dfx_c, fe - 0.5*dfx_p)
    f2 = bcut_getupwind(vy2, fs + 0.5*dfy_n, fp - 0.5*dfy_c)
    f3 = bcut_getupwind(vy3, fp + 0.5*dfy_c, fn - 0.5*dfy_p)
    f4 = bcut_getupwind(vz4, fb + 0.5*dfz_n, fp - 0.5*dfz_c)
    f5 = bcut_getupwind(vz5, fp + 0.5*dfz_c, ft - 0.5*dfz_p)

    q0 = vx0*f0
    q1 = vx1*f1
    q2 = vy2*f2
    q3 = vy3*f3
    q4 = vz4*f4
    q5 = vz5*f5

    m0 = 0.0d0
    m1 = 0.0d0
    m2 = 0.0d0
    m3 = 0.0d0
    m4 = 0.0d0
    m5 = 0.0d0

    if( cidp0 /= 0 ) then
      q0 = (d0 - 0.5d0)/(d0 + 0.5d0)*q1
      m0 = 1.0d0
    endif
    if( cidp1 /= 0 ) then
      q1 = (d1 - 0.5d0)/(d1 + 0.5d0)*q0
      m1 = 1.0d0
    endif
    if( cidp2 /= 0 ) then
      q2 = (d2 - 0.5d0)/(d2 + 0.5d0)*q3
      m2 = 1.0d0
    endif
    if( cidp3 /= 0 ) then
      q3 = (d3 - 0.5d0)/(d3 + 0.5d0)*q2
      m3 = 1.0d0
    endif
    if( cidp4 /= 0 ) then
      q4 = (d4 - 0.5d0)/(d4 + 0.5d0)*q5
      m4 = 1.0d0
    endif
    if( cidp5 /= 0 ) then
      q5 = (d5 - 0.5d0)/(d5 + 0.5d0)*q4
      m5 = 1.0d0
    endif

    if( cidp0 /= 0 .and. cidp1 /= 0 ) then
      q0 = 0.0d0
      q1 = 0.0d0
    endif
    if( cidp2 /= 0 .and. cidp3 /= 0 ) then
      q2 = 0.0d0
      q3 = 0.0d0
    endif
    if( cidp4 /= 0 .and. cidp5 /= 0 ) then
      q4 = 0.0d0
      q5 = 0.0d0
    endif

    fc(i, j, k) = (q1 - q0)/dx &
                + (q3 - q2)/dx &
                + (q5 - q4)/dx

    if( pidp /= 1 ) then
      fc(i, j, k) = 0.0d0
    endif

!    f0 = bcut_getupwind(vx0, fw, fp)
!    f1 = bcut_getupwind(vx1, fp, fe)
!    f2 = bcut_getupwind(vy2, fs, fp)
!    f3 = bcut_getupwind(vy3, fp, fn)
!    f4 = bcut_getupwind(vz4, fb, fp)
!    f5 = bcut_getupwind(vz5, fp, ft)
!
!    f0 = 0.5d0*(fp + fw)
!    f1 = 0.5d0*(fp + fe)
!    f2 = 0.5d0*(fp + fs)
!    f3 = 0.5d0*(fp + fn)
!    f4 = 0.5d0*(fp + fb)
!    f5 = 0.5d0*(fp + ft)
!
!    f0 = 0.5d0*(fp + fw)*0.95 + bcut_getupwind(vx0, fw, fp)*0.05
!    f1 = 0.5d0*(fp + fe)*0.95 + bcut_getupwind(vx1, fp, fe)*0.05
!    f2 = 0.5d0*(fp + fs)*0.95 + bcut_getupwind(vy2, fs, fp)*0.05
!    f3 = 0.5d0*(fp + fn)*0.95 + bcut_getupwind(vy3, fp, fn)*0.05
!    f4 = 0.5d0*(fp + fb)*0.95 + bcut_getupwind(vz4, fb, fp)*0.05
!    f5 = 0.5d0*(fp + ft)*0.95 + bcut_getupwind(vz5, fp, ft)*0.05
!
!    fc(i, j, k) = (f1*vx1 - f0*vx0)/dx &
!                + (f3*vy3 - f2*vy2)/dx &
!                + (f5*vz5 - f4*vz4)/dx
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_calc_c_f_e3

subroutine bcut_calc_c_f_quick( &
                fc, &
                f, &
                vw, ve, vs, vn, vb, vt, &
                c0, c1, c2, c3, c4, c5, &
                cid0, cid1, cid2, cid3, cid4, cid5, &
                pid, &
                dx, dt, &
                fi, &
                sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: fc
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: f
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: vw, ve, vs, vn, vb, vt
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
  integer                  :: cidp
  integer                  :: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
  integer                  :: cidw, cide, cids, cidn, cidb, cidt
  integer                  :: cidw0, cide1, cids2, cidn3, cidb4, cidt5
  integer                  :: pidp, pidw, pide, pids, pidn, pidb, pidt
  real                    :: dx, dt
  real                    :: fi
  real                    :: d0, d1, d2, d3, d4, d5
  real                    :: m0, m1, m2, m3, m4, m5
  real                    :: l0, l1, l2, l3, l4, l5
  real                    :: vx0, vx1, vy2, vy3, vz4, vz5
  real                    :: fp, fw, fe, fs, fn, fb, ft
  real                    :: fww, fee, fss, fnn, fbb, ftt
  real                    :: f0, f1, f2, f3, f4, f5
  real                    :: dfx_p, dfx_c, dfx_n
  real                    :: dfy_p, dfy_c, dfy_n
  real                    :: dfz_p, dfz_c, dfz_n
  real                    :: bcut_getupwind
  real                    :: bcut_getminmod
  real                    :: q0, q1, q2, q3, q4, q5
  real                    :: vx, vy, vz
  real                    :: bcut_getweno3
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp           private(cidp) &
!$omp           private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp           private(cidw, cide, cids, cidn, cidb, cidt) &
!$omp           private(cidw0, cide1, cids2, cidn3, cidb4, cidt5) &
!$omp           private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp           private(d0, d1, d2, d3, d4, d5) &
!$omp           private(m0, m1, m2, m3, m4, m5) &
!$omp           private(l0, l1, l2, l3, l4, l5) &
!$omp           private(vx0, vx1, vy2, vy3, vz4, vz5) &
!$omp           private(fp, fw, fe, fs, fn, fb, ft) &
!$omp           private(fww, fee, fss, fnn, fbb, ftt) &
!$omp           private(f0, f1, f2, f3, f4, f5) &
!$omp           private(dfx_p, dfx_c, dfx_n) &
!$omp           private(dfy_p, dfy_c, dfy_n) &
!$omp           private(dfz_p, dfz_c, dfz_n) &
!$omp           private(q0, q1, q2, q3, q4, q5) &
!$omp           private(vx, vy, vz)
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
    vx0 = vw(i, j, k)
    vx1 = ve(i, j, k)
    vy2 = vs(i, j, k)
    vy3 = vn(i, j, k)
    vz4 = vb(i, j, k)
    vz5 = vt(i, j, k)

    fp = f(i, j, k)
    fw = f(i-1, j, k)
    fe = f(i+1, j, k)
    fs = f(i, j-1, k)
    fn = f(i, j+1, k)
    fb = f(i, j, k-1)
    ft = f(i, j, k+1)
    fww = f(i-2, j, k)
    fee = f(i+2, j, k)
    fss = f(i, j-2, k)
    fnn = f(i, j+2, k)
    fbb = f(i, j, k-2)
    ftt = f(i, j, k+2)

    d0 = c0(i, j, k)
    d1 = c1(i, j, k)
    d2 = c2(i, j, k)
    d3 = c3(i, j, k)
    d4 = c4(i, j, k)
    d5 = c5(i, j, k)

    cidp0 = cid0(i, j, k)
    cidp1 = cid1(i, j, k)
    cidp2 = cid2(i, j, k)
    cidp3 = cid3(i, j, k)
    cidp4 = cid4(i, j, k)
    cidp5 = cid5(i, j, k)

    cidw0 = cid0(i-1, j, k)
    cide1 = cid1(i+1, j, k)
    cids2 = cid2(i, j-1, k)
    cidn3 = cid3(i, j+1, k)
    cidb4 = cid4(i, j, k-1)
    cidt5 = cid5(i, j, k+1)

    pidp = pid(i, j, k)

    f0 = bcut_getupwind(vx0, 0.125*(3.0*fp + 6.0*fw - fww), 0.125*(6.0*fp + 3.0*fw - fe) )
    f1 = bcut_getupwind(vx1, 0.125*(6.0*fp + 3.0*fe - fw) , 0.125*(3.0*fp + 6.0*fe - fee))
    f2 = bcut_getupwind(vy2, 0.125*(3.0*fp + 6.0*fs - fss), 0.125*(6.0*fp + 3.0*fs - fn) )
    f3 = bcut_getupwind(vy3, 0.125*(6.0*fp + 3.0*fn - fs) , 0.125*(3.0*fp + 6.0*fn - fnn))
    f4 = bcut_getupwind(vz4, 0.125*(3.0*fp + 6.0*fb - fbb), 0.125*(6.0*fp + 3.0*fb - ft) )
    f5 = bcut_getupwind(vz5, 0.125*(6.0*fp + 3.0*ft - fb) , 0.125*(3.0*fp + 6.0*ft - ftt))

    if( cidw0 /= 0 ) then
			f0 = 0.5*(fp + fw)
    endif
    if( cide1 /= 0 ) then
			f1 = 0.5*(fp + fe)
    endif
    if( cids2 /= 0 ) then
			f2 = 0.5*(fp + fs)
    endif
    if( cidn3 /= 0 ) then
			f3 = 0.5*(fp + fn)
    endif
    if( cidb4 /= 0 ) then
			f4 = 0.5*(fp + fb)
    endif
    if( cidt5 /= 0 ) then
			f5 = 0.5*(fp + ft)
    endif

    q0 = vx0*f0
    q1 = vx1*f1
    q2 = vy2*f2
    q3 = vy3*f3
    q4 = vz4*f4
    q5 = vz5*f5

    if( cidp0 /= 0 ) then
      q0 = (d0 - 0.5d0)/(d0 + 0.5d0)*q1
    endif
    if( cidp1 /= 0 ) then
      q1 = (d1 - 0.5d0)/(d1 + 0.5d0)*q0
    endif
    if( cidp2 /= 0 ) then
      q2 = (d2 - 0.5d0)/(d2 + 0.5d0)*q3
    endif
    if( cidp3 /= 0 ) then
      q3 = (d3 - 0.5d0)/(d3 + 0.5d0)*q2
    endif
    if( cidp4 /= 0 ) then
      q4 = (d4 - 0.5d0)/(d4 + 0.5d0)*q5
    endif
    if( cidp5 /= 0 ) then
      q5 = (d5 - 0.5d0)/(d5 + 0.5d0)*q4
    endif

    if( cidp0 /= 0 .and. cidp1 /= 0 ) then
      q0 = 0.0d0
      q1 = 0.0d0
    endif
    if( cidp2 /= 0 .and. cidp3 /= 0 ) then
      q2 = 0.0d0
      q3 = 0.0d0
    endif
    if( cidp4 /= 0 .and. cidp5 /= 0 ) then
      q4 = 0.0d0
      q5 = 0.0d0
    endif

    fc(i, j, k) = (q1 - q0)/dx &
                + (q3 - q2)/dx &
                + (q5 - q4)/dx

    if( pidp /= 1 ) then
      fc(i, j, k) = 0.0d0
    endif

  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_calc_c_f_quick

subroutine bcut_calc_c_f_w3( &
                fc, &
                f, &
                vw, ve, vs, vn, vb, vt, &
                c0, c1, c2, c3, c4, c5, &
                cid0, cid1, cid2, cid3, cid4, cid5, &
                pid, &
                dx, dt, &
                fi, &
                sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: fc
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: f
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: vw, ve, vs, vn, vb, vt
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
  integer                  :: cidp
  integer                  :: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
  integer                  :: cidw, cide, cids, cidn, cidb, cidt
  integer                  :: cidw0, cide1, cids2, cidn3, cidb4, cidt5
  integer                  :: pidp, pidw, pide, pids, pidn, pidb, pidt
  real                    :: dx, dt
  real                    :: fi
  real                    :: d0, d1, d2, d3, d4, d5
  real                    :: m0, m1, m2, m3, m4, m5
  real                    :: l0, l1, l2, l3, l4, l5
  real                    :: vx0, vx1, vy2, vy3, vz4, vz5
  real                    :: fp, fw, fe, fs, fn, fb, ft
  real                    :: fww, fee, fss, fnn, fbb, ftt
  real                    :: f0, f1, f2, f3, f4, f5
  real                    :: dfx_p, dfx_c, dfx_n
  real                    :: dfy_p, dfy_c, dfy_n
  real                    :: dfz_p, dfz_c, dfz_n
  real                    :: bcut_getupwind
  real                    :: bcut_getminmod
  real                    :: q0, q1, q2, q3, q4, q5
  real                    :: vx, vy, vz
  real                    :: bcut_getweno3
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp           private(cidp) &
!$omp           private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp           private(cidw, cide, cids, cidn, cidb, cidt) &
!$omp           private(cidw0, cide1, cids2, cidn3, cidb4, cidt5) &
!$omp           private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp           private(d0, d1, d2, d3, d4, d5) &
!$omp           private(m0, m1, m2, m3, m4, m5) &
!$omp           private(l0, l1, l2, l3, l4, l5) &
!$omp           private(vx0, vx1, vy2, vy3, vz4, vz5) &
!$omp           private(fp, fw, fe, fs, fn, fb, ft) &
!$omp           private(fww, fee, fss, fnn, fbb, ftt) &
!$omp           private(f0, f1, f2, f3, f4, f5) &
!$omp           private(dfx_p, dfx_c, dfx_n) &
!$omp           private(dfy_p, dfy_c, dfy_n) &
!$omp           private(dfz_p, dfz_c, dfz_n) &
!$omp           private(q0, q1, q2, q3, q4, q5) &
!$omp           private(vx, vy, vz)
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
    vx0 = vw(i, j, k)
    vx1 = ve(i, j, k)
    vy2 = vs(i, j, k)
    vy3 = vn(i, j, k)
    vz4 = vb(i, j, k)
    vz5 = vt(i, j, k)

    fp = f(i, j, k)
    fw = f(i-1, j, k)
    fe = f(i+1, j, k)
    fs = f(i, j-1, k)
    fn = f(i, j+1, k)
    fb = f(i, j, k-1)
    ft = f(i, j, k+1)
    fww = f(i-2, j, k)
    fee = f(i+2, j, k)
    fss = f(i, j-2, k)
    fnn = f(i, j+2, k)
    fbb = f(i, j, k-2)
    ftt = f(i, j, k+2)

    m0 = 0.0d0
    m1 = 0.0d0
    m2 = 0.0d0
    m3 = 0.0d0
    m4 = 0.0d0
    m5 = 0.0d0

    d0 = c0(i, j, k)
    d1 = c1(i, j, k)
    d2 = c2(i, j, k)
    d3 = c3(i, j, k)
    d4 = c4(i, j, k)
    d5 = c5(i, j, k)

    cidp0 = cid0(i, j, k)
    cidp1 = cid1(i, j, k)
    cidp2 = cid2(i, j, k)
    cidp3 = cid3(i, j, k)
    cidp4 = cid4(i, j, k)
    cidp5 = cid5(i, j, k)

    cidw0 = cid0(i-1, j, k)
    cide1 = cid1(i+1, j, k)
    cids2 = cid2(i, j-1, k)
    cidn3 = cid3(i, j+1, k)
    cidb4 = cid4(i, j, k-1)
    cidt5 = cid5(i, j, k+1)

    pidp = pid(i, j, k)

    if( cidp0 /= 0 ) then
      fw  = (1.0d0 - 1.0d0/d0)*fp + (1.0d0/d0)*fi
      fww = (1.0d0 - 2.0d0/d0)*fp + (2.0d0/d0)*fi
      m0 = 1.0d0
    endif
    if( cidp1 /= 0 ) then
      fe  = (1.0d0 - 1.0d0/d1)*fp + (1.0d0/d1)*fi
      fee = (1.0d0 - 2.0d0/d1)*fp + (2.0d0/d1)*fi
      m1 = 1.0d0
    endif
    if( cidp2 /= 0 ) then
      fs  = (1.0d0 - 1.0d0/d2)*fp + (1.0d0/d2)*fi
      fss = (1.0d0 - 2.0d0/d2)*fp + (2.0d0/d2)*fi
      m2 = 1.0d0
    endif
    if( cidp3 /= 0 ) then
      fn  = (1.0d0 - 1.0d0/d3)*fp + (1.0d0/d3)*fi
      fnn = (1.0d0 - 2.0d0/d3)*fp + (2.0d0/d3)*fi
      m3 = 1.0d0
    endif
    if( cidp4 /= 0 ) then
      fb  = (1.0d0 - 1.0d0/d4)*fp + (1.0d0/d4)*fi
      fbb = (1.0d0 - 2.0d0/d4)*fp + (2.0d0/d4)*fi
      m4 = 1.0d0
    endif
    if( cidp5 /= 0 ) then
      ft  = (1.0d0 - 1.0d0/d5)*fp + (1.0d0/d5)*fi
      ftt = (1.0d0 - 2.0d0/d5)*fp + (2.0d0/d5)*fi
      m5 = 1.0d0
    endif

    if( cidw0 /= 0 ) then
      fww = fi
    endif
    if( cide1 /= 0 ) then
      fee = fi
    endif
    if( cids2 /= 0 ) then
      fss = fi
    endif
    if( cidn3 /= 0 ) then
      fnn = fi
    endif
    if( cidb4 /= 0 ) then
      fbb = fi
    endif
    if( cidt5 /= 0 ) then
      ftt = fi
    endif

    dfx_p = bcut_getweno3(fee - fe, fe - fp, fp - fw)/dx
    dfx_n = bcut_getweno3(fw - fww, fp - fw, fe - fp)/dx
    dfy_p = bcut_getweno3(fnn - fn, fn - fp, fp - fs)/dx
    dfy_n = bcut_getweno3(fs - fss, fp - fs, fn - fp)/dx
    dfz_p = bcut_getweno3(ftt - ft, ft - fp, fp - fb)/dx
    dfz_n = bcut_getweno3(fb - fbb, fp - fb, ft - fp)/dx

    vx = 0.5*(vx0 + vx1)
    vy = 0.5*(vy2 + vy3)
    vz = 0.5*(vz4 + vz5)

    fc(i, j, k) = (max(vx, 0.0)*dfx_n + min(vx, 0.0)*dfx_p) &
                + (max(vy, 0.0)*dfy_n + min(vy, 0.0)*dfy_p) &
                + (max(vz, 0.0)*dfz_n + min(vz, 0.0)*dfz_p) 


!    f0 = bcut_getupwind(vx0, fw + 0.5*dfx_n, fp - 0.5*dfx_c)
!    f1 = bcut_getupwind(vx1, fp + 0.5*dfx_c, fe - 0.5*dfx_p)
!    f2 = bcut_getupwind(vy2, fs + 0.5*dfy_n, fp - 0.5*dfy_c)
!    f3 = bcut_getupwind(vy3, fp + 0.5*dfy_c, fn - 0.5*dfy_p)
!    f4 = bcut_getupwind(vz4, fb + 0.5*dfz_n, fp - 0.5*dfz_c)
!    f5 = bcut_getupwind(vz5, fp + 0.5*dfz_c, ft - 0.5*dfz_p)
!
!    q0 = vx0*f0
!    q1 = vx1*f1
!    q2 = vy2*f2
!    q3 = vy3*f3
!    q4 = vz4*f4
!    q5 = vz5*f5
!
!    m0 = 0.0d0
!    m1 = 0.0d0
!    m2 = 0.0d0
!    m3 = 0.0d0
!    m4 = 0.0d0
!    m5 = 0.0d0
!
!    if( cidp0 /= 0 ) then
!      q0 = (d0 - 0.5d0)/(d0 + 0.5d0)*q1
!      m0 = 1.0d0
!    endif
!    if( cidp1 /= 0 ) then
!      q1 = (d1 - 0.5d0)/(d1 + 0.5d0)*q0
!      m1 = 1.0d0
!    endif
!    if( cidp2 /= 0 ) then
!      q2 = (d2 - 0.5d0)/(d2 + 0.5d0)*q3
!      m2 = 1.0d0
!    endif
!    if( cidp3 /= 0 ) then
!      q3 = (d3 - 0.5d0)/(d3 + 0.5d0)*q2
!      m3 = 1.0d0
!    endif
!    if( cidp4 /= 0 ) then
!      q4 = (d4 - 0.5d0)/(d4 + 0.5d0)*q5
!      m4 = 1.0d0
!    endif
!    if( cidp5 /= 0 ) then
!      q5 = (d5 - 0.5d0)/(d5 + 0.5d0)*q4
!      m5 = 1.0d0
!    endif
!
!    fc(i, j, k) = (q1 - q0)/dx &
!                + (q3 - q2)/dx &
!                + (q5 - q4)/dx

    if( pidp /= 1 ) then
      fc(i, j, k) = 0.0d0
    endif

  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_calc_c_f_w3

subroutine bcut_calc_c_f_blend( &
                fc, &
                f, &
                vw, ve, vs, vn, vb, vt, &
                c0, c1, c2, c3, c4, c5, &
                cid0, cid1, cid2, cid3, cid4, cid5, &
                pid, &
                dx, dt, &
                fi, &
                alpha, &
                sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: fc
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: f
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: vw, ve, vs, vn, vb, vt
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
  integer                  :: cidp
  integer                  :: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
  integer                  :: cidw, cide, cids, cidn, cidb, cidt
  integer                  :: cidw0, cide1, cids2, cidn3, cidb4, cidt5
  integer                  :: pidp, pidw, pide, pids, pidn, pidb, pidt
  real                    :: dx, dt
  real                    :: fi
  real                    :: alpha
  real                    :: d0, d1, d2, d3, d4, d5
  real                    :: m0, m1, m2, m3, m4, m5
  real                    :: l0, l1, l2, l3, l4, l5
  real                    :: vx0, vx1, vy2, vy3, vz4, vz5
  real                    :: fp, fw, fe, fs, fn, fb, ft
  real                    :: fww, fee, fss, fnn, fbb, ftt
  real                    :: f0, f1, f2, f3, f4, f5
  real                    :: dfx_p, dfx_c, dfx_n
  real                    :: dfy_p, dfy_c, dfy_n
  real                    :: dfz_p, dfz_c, dfz_n
  real                    :: bcut_getupwind
  real                    :: bcut_getminmod
  real                    :: q0, q1, q2, q3, q4, q5
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp           private(cidp) &
!$omp           private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp           private(cidw, cide, cids, cidn, cidb, cidt) &
!$omp           private(cidw0, cide1, cids2, cidn3, cidb4, cidt5) &
!$omp           private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp           private(d0, d1, d2, d3, d4, d5) &
!$omp           private(m0, m1, m2, m3, m4, m5) &
!$omp           private(l0, l1, l2, l3, l4, l5) &
!$omp           private(vx0, vx1, vy2, vy3, vz4, vz5) &
!$omp           private(fp, fw, fe, fs, fn, fb, ft) &
!$omp           private(fww, fee, fss, fnn, fbb, ftt) &
!$omp           private(f0, f1, f2, f3, f4, f5) &
!$omp           private(dfx_p, dfx_c, dfx_n) &
!$omp           private(dfy_p, dfy_c, dfy_n) &
!$omp           private(dfz_p, dfz_c, dfz_n) &
!$omp           private(q0, q1, q2, q3, q4, q5)
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
    vx0 = vw(i, j, k)
    vx1 = ve(i, j, k)
    vy2 = vs(i, j, k)
    vy3 = vn(i, j, k)
    vz4 = vb(i, j, k)
    vz5 = vt(i, j, k)

    fp = f(i, j, k)
    fw = f(i-1, j, k)
    fe = f(i+1, j, k)
    fs = f(i, j-1, k)
    fn = f(i, j+1, k)
    fb = f(i, j, k-1)
    ft = f(i, j, k+1)

    d0 = c0(i, j, k)
    d1 = c1(i, j, k)
    d2 = c2(i, j, k)
    d3 = c3(i, j, k)
    d4 = c4(i, j, k)
    d5 = c5(i, j, k)

    cidp0 = cid0(i, j, k)
    cidp1 = cid1(i, j, k)
    cidp2 = cid2(i, j, k)
    cidp3 = cid3(i, j, k)
    cidp4 = cid4(i, j, k)
    cidp5 = cid5(i, j, k)

    pidp = pid(i, j, k)

    f0 = 0.5d0*(fp + fw)*alpha + bcut_getupwind(vx0, fw, fp)*(1.0d0 - alpha)
    f1 = 0.5d0*(fp + fe)*alpha + bcut_getupwind(vx1, fp, fe)*(1.0d0 - alpha)
    f2 = 0.5d0*(fp + fs)*alpha + bcut_getupwind(vy2, fs, fp)*(1.0d0 - alpha)
    f3 = 0.5d0*(fp + fn)*alpha + bcut_getupwind(vy3, fp, fn)*(1.0d0 - alpha)
    f4 = 0.5d0*(fp + fb)*alpha + bcut_getupwind(vz4, fb, fp)*(1.0d0 - alpha)
    f5 = 0.5d0*(fp + ft)*alpha + bcut_getupwind(vz5, fp, ft)*(1.0d0 - alpha)

    q0 = vx0*f0
    q1 = vx1*f1
    q2 = vy2*f2
    q3 = vy3*f3
    q4 = vz4*f4
    q5 = vz5*f5

    m0 = 0.0d0
    m1 = 0.0d0
    m2 = 0.0d0
    m3 = 0.0d0
    m4 = 0.0d0
    m5 = 0.0d0

    if( cidp0 /= 0 ) then
      q0 = (d0 - 0.5d0)/(d0 + 0.5d0)*q1
      m0 = 1.0d0
    endif
    if( cidp1 /= 0 ) then
      q1 = (d1 - 0.5d0)/(d1 + 0.5d0)*q0
      m1 = 1.0d0
    endif
    if( cidp2 /= 0 ) then
      q2 = (d2 - 0.5d0)/(d2 + 0.5d0)*q3
      m2 = 1.0d0
    endif
    if( cidp3 /= 0 ) then
      q3 = (d3 - 0.5d0)/(d3 + 0.5d0)*q2
      m3 = 1.0d0
    endif
    if( cidp4 /= 0 ) then
      q4 = (d4 - 0.5d0)/(d4 + 0.5d0)*q5
      m4 = 1.0d0
    endif
    if( cidp5 /= 0 ) then
      q5 = (d5 - 0.5d0)/(d5 + 0.5d0)*q4
      m5 = 1.0d0
    endif

    if( cidp0 /= 0 .and. cidp1 /= 0 ) then
      q0 = 0.0d0
      q1 = 0.0d0
    endif
    if( cidp2 /= 0 .and. cidp3 /= 0 ) then
      q2 = 0.0d0
      q3 = 0.0d0
    endif
    if( cidp4 /= 0 .and. cidp5 /= 0 ) then
      q4 = 0.0d0
      q5 = 0.0d0
    endif

    fc(i, j, k) = (q1 - q0)/dx &
                + (q3 - q2)/dx &
                + (q5 - q4)/dx

    if( pidp /= 1 ) then
      fc(i, j, k) = 0.0d0
    endif

  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_calc_c_f_blend

subroutine bcut_calc_c_u_quick( &
                fc, &
                f, &
                vw, ve, vs, vn, vb, vt, &
                c0, c1, c2, c3, c4, c5, &
                cid0, cid1, cid2, cid3, cid4, cid5, &
                pid, &
                dx, dt, &
                fi, &
                sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: fc
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: f
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: vw, ve, vs, vn, vb, vt
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
  integer                  :: cidp
  integer                  :: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
  integer                  :: cidw, cide, cids, cidn, cidb, cidt
  integer                  :: cidw0, cide1, cids2, cidn3, cidb4, cidt5
  integer                  :: pidp, pidw, pide, pids, pidn, pidb, pidt
  real                    :: dx, dt
  real                    :: fi
  real                    :: d0, d1, d2, d3, d4, d5
  real                    :: m0, m1, m2, m3, m4, m5
  real                    :: l0, l1, l2, l3, l4, l5
  real                    :: vx0, vx1, vy2, vy3, vz4, vz5
  real                    :: fp, fw, fe, fs, fn, fb, ft
  real                    :: fww, fee, fss, fnn, fbb, ftt
  real                    :: f0, f1, f2, f3, f4, f5
  real                    :: dfx_p, dfx_c, dfx_n
  real                    :: dfy_p, dfy_c, dfy_n
  real                    :: dfz_p, dfz_c, dfz_n
  real                    :: bcut_getupwind
  real                    :: bcut_getminmod
  real                    :: q0, q1, q2, q3, q4, q5
  real                    :: vx, vy, vz
  real                    :: bcut_getweno3
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp           private(cidp) &
!$omp           private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp           private(cidw, cide, cids, cidn, cidb, cidt) &
!$omp           private(cidw0, cide1, cids2, cidn3, cidb4, cidt5) &
!$omp           private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp           private(d0, d1, d2, d3, d4, d5) &
!$omp           private(m0, m1, m2, m3, m4, m5) &
!$omp           private(l0, l1, l2, l3, l4, l5) &
!$omp           private(vx0, vx1, vy2, vy3, vz4, vz5) &
!$omp           private(fp, fw, fe, fs, fn, fb, ft) &
!$omp           private(fww, fee, fss, fnn, fbb, ftt) &
!$omp           private(f0, f1, f2, f3, f4, f5) &
!$omp           private(dfx_p, dfx_c, dfx_n) &
!$omp           private(dfy_p, dfy_c, dfy_n) &
!$omp           private(dfz_p, dfz_c, dfz_n) &
!$omp           private(q0, q1, q2, q3, q4, q5) &
!$omp           private(vx, vy, vz)
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
    vx0 = vw(i, j, k)
    vx1 = ve(i, j, k)
    vy2 = vs(i, j, k)
    vy3 = vn(i, j, k)
    vz4 = vb(i, j, k)
    vz5 = vt(i, j, k)

    fp = f(i, j, k)
    fw = f(i-1, j, k)
    fe = f(i+1, j, k)
    fs = f(i, j-1, k)
    fn = f(i, j+1, k)
    fb = f(i, j, k-1)
    ft = f(i, j, k+1)
    fww = f(i-2, j, k)
    fee = f(i+2, j, k)
    fss = f(i, j-2, k)
    fnn = f(i, j+2, k)
    fbb = f(i, j, k-2)
    ftt = f(i, j, k+2)

    d0 = c0(i, j, k)
    d1 = c1(i, j, k)
    d2 = c2(i, j, k)
    d3 = c3(i, j, k)
    d4 = c4(i, j, k)
    d5 = c5(i, j, k)

    m0 = 0.0d0
    m1 = 0.0d0
    m2 = 0.0d0
    m3 = 0.0d0
    m4 = 0.0d0
    m5 = 0.0d0

    cidp0 = cid0(i, j, k)
    cidp1 = cid1(i, j, k)
    cidp2 = cid2(i, j, k)
    cidp3 = cid3(i, j, k)
    cidp4 = cid4(i, j, k)
    cidp5 = cid5(i, j, k)

    cidw0 = cid0(i-1, j, k)
    cide1 = cid1(i+1, j, k)
    cids2 = cid2(i, j-1, k)
    cidn3 = cid3(i, j+1, k)
    cidb4 = cid4(i, j, k-1)
    cidt5 = cid5(i, j, k+1)

    pidp = pid(i, j, k)

    if( cidp0 /= 0 ) then
			fw  = 0.0
			fww = 0.0
    endif
    if( cidp1 /= 0 ) then
			fe  = 0.0
			fee = 0.0
    endif
    if( cidp2 /= 0 ) then
			fs  = 0.0
			fss = 0.0
    endif
    if( cidp3 /= 0 ) then
			fn  = 0.0
			fnn = 0.0
    endif
    if( cidp4 /= 0 ) then
			fb  = 0.0
			fbb = 0.0
    endif
    if( cidp5 /= 0 ) then
			ft  = 0.0
			ftt = 0.0
    endif

    if( cidw0 /= 0 ) then
			fww = fw
    endif
    if( cide1 /= 0 ) then
			fee = fe
    endif
    if( cids2 /= 0 ) then
			fss = fs
    endif
    if( cidn3 /= 0 ) then
			fnn = fn
    endif
    if( cidb4 /= 0 ) then
			fbb = fb
    endif
    if( cidt5 /= 0 ) then
			ftt = ft
    endif


    if( cidp0 /= 0 ) then
      fw  = (1.0d0 - 1.0d0/d0)*fp + (1.0d0/d0)*fi
      fww = (1.0d0 - 2.0d0/d0)*fp + (2.0d0/d0)*fi
      m0 = 1.0d0
    endif
    if( cidp1 /= 0 ) then
      fe  = (1.0d0 - 1.0d0/d1)*fp + (1.0d0/d1)*fi
      fee = (1.0d0 - 2.0d0/d1)*fp + (2.0d0/d1)*fi
      m1 = 1.0d0
    endif
    if( cidp2 /= 0 ) then
      fs  = (1.0d0 - 1.0d0/d2)*fp + (1.0d0/d2)*fi
      fss = (1.0d0 - 2.0d0/d2)*fp + (2.0d0/d2)*fi
      m2 = 1.0d0
    endif
    if( cidp3 /= 0 ) then
      fn  = (1.0d0 - 1.0d0/d3)*fp + (1.0d0/d3)*fi
      fnn = (1.0d0 - 2.0d0/d3)*fp + (2.0d0/d3)*fi
      m3 = 1.0d0
    endif
    if( cidp4 /= 0 ) then
      fb  = (1.0d0 - 1.0d0/d4)*fp + (1.0d0/d4)*fi
      fbb = (1.0d0 - 2.0d0/d4)*fp + (2.0d0/d4)*fi
      m4 = 1.0d0
    endif
    if( cidp5 /= 0 ) then
      ft  = (1.0d0 - 1.0d0/d5)*fp + (1.0d0/d5)*fi
      ftt = (1.0d0 - 2.0d0/d5)*fp + (2.0d0/d5)*fi
      m5 = 1.0d0
    endif

    if( cidw0 /= 0 ) then
      fww = fi
    endif
    if( cide1 /= 0 ) then
      fee = fi
    endif
    if( cids2 /= 0 ) then
      fss = fi
    endif
    if( cidn3 /= 0 ) then
      fnn = fi
    endif
    if( cidb4 /= 0 ) then
      fbb = fi
    endif
    if( cidt5 /= 0 ) then
      ftt = fi
    endif

    f0 = bcut_getupwind(vx0, 0.125*(3.0*fp + 6.0*fw - fww), 0.125*(6.0*fp + 3.0*fw - fe) )
    f1 = bcut_getupwind(vx1, 0.125*(6.0*fp + 3.0*fe - fw) , 0.125*(3.0*fp + 6.0*fe - fee))
    f2 = bcut_getupwind(vy2, 0.125*(3.0*fp + 6.0*fs - fss), 0.125*(6.0*fp + 3.0*fs - fn) )
    f3 = bcut_getupwind(vy3, 0.125*(6.0*fp + 3.0*fn - fs) , 0.125*(3.0*fp + 6.0*fn - fnn))
    f4 = bcut_getupwind(vz4, 0.125*(3.0*fp + 6.0*fb - fbb), 0.125*(6.0*fp + 3.0*fb - ft) )
    f5 = bcut_getupwind(vz5, 0.125*(6.0*fp + 3.0*ft - fb) , 0.125*(3.0*fp + 6.0*ft - ftt))

    q0 = vx0*f0
    q1 = vx1*f1
    q2 = vy2*f2
    q3 = vy3*f3
    q4 = vz4*f4
    q5 = vz5*f5

    if( cidp0 /= 0 .and. cidp1 /= 0 ) then
      q0 = 0.0d0
      q1 = 0.0d0
    endif
    if( cidp2 /= 0 .and. cidp3 /= 0 ) then
      q2 = 0.0d0
      q3 = 0.0d0
    endif
    if( cidp4 /= 0 .and. cidp5 /= 0 ) then
      q4 = 0.0d0
      q5 = 0.0d0
    endif

    fc(i, j, k) = (q1 - q0)/dx &
                + (q3 - q2)/dx &
                + (q5 - q4)/dx

    if( pidp /= 1 ) then
      fc(i, j, k) = 0.0d0
    endif

  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_calc_c_u_quick

subroutine bcut_calc_d_u( &
                ud0_, &
                u0_, &
                c0, c1, c2, c3, c4, c5, &
                cid0, cid1, cid2, cid3, cid4, cid5, &
                pid, &
                rho, &
                mu, &
                dx, dt, &
                Uc, &
								eps, &
                sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: ud0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: u0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
  integer                  :: cidp
  integer                  :: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
  integer                  :: pidp, pidw, pide, pids, pidn, pidb, pidt
  real                    :: rho
  real                    :: mu
  real                    :: dx, dt
  real                    :: Uc
	real										:: eps
  real                    :: d0, d1, d2, d3, d4, d5
  real                    :: mu0, mu1, mu2, mu3, mu4, mu5
  real                    :: m0, m1, m2, m3, m4, m5
  real                    :: l0, l1, l2, l3, l4, l5
  real                    :: up, uw, ue, us, un, ub, ut
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp           private(cidp) &
!$omp           private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp           private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp           private(d0, d1, d2, d3, d4, d5) &
!$omp           private(mu0, mu1, mu2, mu3, mu4, mu5) &
!$omp           private(m0, m1, m2, m3, m4, m5) &
!$omp           private(l0, l1, l2, l3, l4, l5) &
!$omp           private(up, uw, ue, us, un, ub, ut) 
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
!dir$ simd
  do i=1, ix
    mu0 = mu
    mu1 = mu
    mu2 = mu
    mu3 = mu
    mu4 = mu
    mu5 = mu

    m0 = 0.0d0
    m1 = 0.0d0
    m2 = 0.0d0
    m3 = 0.0d0
    m4 = 0.0d0
    m5 = 0.0d0

    d0 = c0(i, j, k)
    d1 = c1(i, j, k)
    d2 = c2(i, j, k)
    d3 = c3(i, j, k)
    d4 = c4(i, j, k)
    d5 = c5(i, j, k)

    cidp0 = cid0(i, j, k)
    cidp1 = cid1(i, j, k)
    cidp2 = cid2(i, j, k)
    cidp3 = cid3(i, j, k)
    cidp4 = cid4(i, j, k)
    cidp5 = cid5(i, j, k)

    pidp = pid(i, j, k)

    if( cidp0 /= 0 ) then
      mu0 = mu/d0*2.0d0/(d0 + d1)
      mu1 = mu/d1*2.0d0/(d0 + d1)
			if( d0 < eps ) then
				mu0 = 0.0
				mu1 = mu
			end if
      m0 = 1.0d0
    endif

    if( cidp1 /= 0 ) then
      mu0 = mu/d0*2.0d0/(d0 + d1)
      mu1 = mu/d1*2.0d0/(d0 + d1)
			if( d1 < eps ) then
				mu0 = mu
				mu1 = 0.0
			end if
      m1 = 1.0d0
    endif

    if( cidp2 /= 0 ) then
      mu2 = mu/d2*2.0d0/(d2 + d3)
      mu3 = mu/d3*2.0d0/(d2 + d3)
			if( d2 < eps ) then
				mu2 = 0.0
				mu3 = mu
			end if
      m2 = 1.0d0
    endif

    if( cidp3 /= 0 ) then
      mu2 = mu/d2*2.0d0/(d2 + d3)
      mu3 = mu/d3*2.0d0/(d2 + d3)
			if( d3 < eps ) then
				mu2 = mu
				mu3 = 0.0
			end if
      m3 = 1.0d0
    endif

    if( cidp4 /= 0 ) then
      mu4 = mu/d4*2.0d0/(d4 + d5)
      mu5 = mu/d5*2.0d0/(d4 + d5)
			if( d4 < eps ) then
				mu4 = 0.0
				mu5 = mu
			end if
      m4 = 1.0d0
    endif

    if( cidp5 /= 0 ) then
      mu4 = mu/d4*2.0d0/(d4 + d5)
      mu5 = mu/d5*2.0d0/(d4 + d5)
			if( d5 < eps ) then
				mu4 = mu
				mu5 = 0.0
			end if
      m5 = 1.0d0
    endif

    l0 = mu0/(rho)/(dx*dx)
    l1 = mu1/(rho)/(dx*dx)
    l2 = mu2/(rho)/(dx*dx)
    l3 = mu3/(rho)/(dx*dx)
    l4 = mu4/(rho)/(dx*dx)
    l5 = mu5/(rho)/(dx*dx)

    up = u0_(i, j, k)
    uw = u0_(i-1, j, k)*(1.0 - m0) + Uc*m0
    ue = u0_(i+1, j, k)*(1.0 - m1) + Uc*m1
    us = u0_(i, j-1, k)*(1.0 - m2) + Uc*m2
    un = u0_(i, j+1, k)*(1.0 - m3) + Uc*m3
    ub = u0_(i, j, k-1)*(1.0 - m4) + Uc*m4
    ut = u0_(i, j, k+1)*(1.0 - m5) + Uc*m5

    ud0_(i, j, k) = ( &
                        l1*(ue - up) &
                      - l0*(up - uw) &
                      + l3*(un - up) &
                      - l2*(up - us) &
                      + l5*(ut - up) &
                      - l4*(up - ub) &
                    )

    if( pidp /= 1 ) then
      ud0_(i, j, k) = 0.0d0
    endif
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_calc_d_u

subroutine bcut_calc_d_u_e( &
                ud0_, &
                u0_, &
								nue_, &
                c0, c1, c2, c3, c4, c5, &
                cid0, cid1, cid2, cid3, cid4, cid5, &
                pid, &
                rho, &
                mu, &
                dx, dt, &
                Uc, &
                sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: ud0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: u0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: nue_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
  integer                  :: cidp
  integer                  :: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
  integer                  :: pidp, pidw, pide, pids, pidn, pidb, pidt
  real                    :: rho
  real                    :: mu
  real                    :: dx, dt
  real                    :: Uc
  real                    :: d0, d1, d2, d3, d4, d5
  real                    :: nup, nu0, nu1, nu2, nu3, nu4, nu5
  real                    :: ap, a0, a1, a2, a3, a4, a5
  real                    :: l0, l1, l2, l3, l4, l5
  real                    :: m0, m1, m2, m3, m4, m5
  real                    :: up, uw, ue, us, un, ub, ut
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp           private(cidp) &
!$omp           private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp           private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp           private(d0, d1, d2, d3, d4, d5) &
!$omp           private(nup, nu0, nu1, nu2, nu3, nu4, nu5) &
!$omp           private(m0, m1, m2, m3, m4, m5) &
!$omp           private(l0, l1, l2, l3, l4, l5) &
!$omp           private(up, uw, ue, us, un, ub, ut) 
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
!dir$ simd
  do i=1, ix
    d0 = c0(i, j, k)
    d1 = c1(i, j, k)
    d2 = c2(i, j, k)
    d3 = c3(i, j, k)
    d4 = c4(i, j, k)
    d5 = c5(i, j, k)

    cidp0 = cid0(i, j, k)
    cidp1 = cid1(i, j, k)
    cidp2 = cid2(i, j, k)
    cidp3 = cid3(i, j, k)
    cidp4 = cid4(i, j, k)
    cidp5 = cid5(i, j, k)

    pidp = pid(i, j, k)

    if( pidp /= 1 ) then
      ud0_(i, j, k) = 0.0d0
			cycle
    endif

		nup = nue_(i, j, k)
    nu0 = nue_(i-1, j, k)
    nu1 = nue_(i+1, j, k)
    nu2 = nue_(i, j-1, k)
    nu3 = nue_(i, j+1, k)
    nu4 = nue_(i, j, k-1)
    nu5 = nue_(i, j, k+1)

		a0 = nup*nu0/(nup + nu0)*2.0
		a1 = nup*nu1/(nup + nu1)*2.0
		a2 = nup*nu2/(nup + nu2)*2.0
		a3 = nup*nu3/(nup + nu3)*2.0
		a4 = nup*nu4/(nup + nu4)*2.0
		a5 = nup*nu5/(nup + nu5)*2.0

    m0 = 0.0d0
    m1 = 0.0d0
    m2 = 0.0d0
    m3 = 0.0d0
    m4 = 0.0d0
    m5 = 0.0d0
    if( cidp0 /= 0 ) then
      m0 = 1.0d0
			a0 = nup*(d0 - 0.5)/d0
    endif
    if( cidp1 /= 0 ) then
      m1 = 1.0d0
			a1 = nup*(d1 - 0.5)/d1
    endif
    if( cidp2 /= 0 ) then
      m2 = 1.0d0
			a2 = nup*(d2 - 0.5)/d2
    endif
    if( cidp3 /= 0 ) then
      m3 = 1.0d0
			a3 = nup*(d3 - 0.5)/d3
    endif
    if( cidp4 /= 0 ) then
      m4 = 1.0d0
			a4 = nup*(d4 - 0.5)/d4
    endif
    if( cidp5 /= 0 ) then
      m5 = 1.0d0
			a5 = nup*(d5 - 0.5)/d5
    endif

    l0 = a0/(dx*dx)
    l1 = a1/(dx*dx)
    l2 = a2/(dx*dx)
    l3 = a3/(dx*dx)
    l4 = a4/(dx*dx)
    l5 = a5/(dx*dx)

    up = u0_(i, j, k)
    uw = u0_(i-1, j, k)*(1.0 - m0) + Uc*m0
    ue = u0_(i+1, j, k)*(1.0 - m1) + Uc*m1
    us = u0_(i, j-1, k)*(1.0 - m2) + Uc*m2
    un = u0_(i, j+1, k)*(1.0 - m3) + Uc*m3
    ub = u0_(i, j, k-1)*(1.0 - m4) + Uc*m4
    ut = u0_(i, j, k+1)*(1.0 - m5) + Uc*m5

    ud0_(i, j, k) = ( &
                        l1*(ue - up) &
                      - l0*(up - uw) &
                      + l3*(un - up) &
                      - l2*(up - us) &
                      + l5*(ut - up) &
                      - l4*(up - ub) &
                    )

  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_calc_d_u_e

subroutine bcut_calc_nue( &
								nue_, &
								csm_c2, &
                ux_, &
								uy_, &
								uz_, &
                c0, c1, c2, c3, c4, c5, &
                cid0, cid1, cid2, cid3, cid4, cid5, &
                pid, &
                rho, &
                mu, &
                dx, dt, &
                Uc, Vc, Wc, &
                sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: nue_
	real										:: csm_c2
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: ux_, uy_, uz_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
  integer                  :: cidp
  integer                  :: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
  integer                  :: pidp, pidw, pide, pids, pidn, pidb, pidt
  real                    :: rho
  real                    :: mu
  real                    :: dx, dt
  real                    :: Uc, Vc, Wc
  real                    :: d0, d1, d2, d3, d4, d5
  real                    :: mu0, mu1, mu2, mu3, mu4, mu5
  real                    :: m0, m1, m2, m3, m4, m5
  real                    :: l0, l1, l2, l3, l4, l5
  real                    :: up, uw, ue, us, un, ub, ut
  real                    :: vp, vw, ve, vs, vn, vb, vt
  real                    :: wp, ww, we, ws, wn, wb, wt
	real										:: dudx, dudy, dudz
	real										:: dvdx, dvdy, dvdz
	real										:: dwdx, dwdy, dwdz
	real										:: sxx, sxy, sxz, syx, syy, syz, szx, szy, szz
	real										:: s2, sl
	real										:: wxx, wxy, wxz, wyx, wyy, wyz, wzx, wzy, wzz
	real										:: csm_e, csm_q, csm_fcs, csm_fom, csm_c
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp           private(cidp) &
!$omp           private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp           private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp           private(d0, d1, d2, d3, d4, d5) &
!$omp           private(mu0, mu1, mu2, mu3, mu4, mu5) &
!$omp           private(m0, m1, m2, m3, m4, m5) &
!$omp           private(l0, l1, l2, l3, l4, l5) &
!$omp           private(up, uw, ue, us, un, ub, ut) &
!$omp           private(vp, vw, ve, vs, vn, vb, vt) &
!$omp           private(wp, ww, we, ws, wn, wb, wt) &
!$omp						private(dudx, dudy, dudz) &
!$omp						private(dvdx, dvdy, dvdz) &
!$omp						private(dwdx, dwdy, dwdz) &
!$omp						private(sxx, sxy, sxz, syx, syy, syz, szx, szy, szz) &
!$omp						private(s2, sl) &
!$omp						private(wxx, wxy, wxz, wyx, wyy, wyz, wzx, wzy, wzz) &
!$omp						private(csm_e, csm_q, csm_fcs, csm_fom, csm_c) 
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
!dir$ simd
  do i=1, ix
    d0 = c0(i, j, k)
    d1 = c1(i, j, k)
    d2 = c2(i, j, k)
    d3 = c3(i, j, k)
    d4 = c4(i, j, k)
    d5 = c5(i, j, k)

    cidp0 = cid0(i, j, k)
    cidp1 = cid1(i, j, k)
    cidp2 = cid2(i, j, k)
    cidp3 = cid3(i, j, k)
    cidp4 = cid4(i, j, k)
    cidp5 = cid5(i, j, k)

    pidp = pid(i, j, k)

    m0 = 0.0d0
    m1 = 0.0d0
    m2 = 0.0d0
    m3 = 0.0d0
    m4 = 0.0d0
    m5 = 0.0d0
    if( cidp0 /= 0 ) then
      m0 = 1.0d0
    endif
    if( cidp1 /= 0 ) then
      m1 = 1.0d0
    endif
    if( cidp2 /= 0 ) then
      m2 = 1.0d0
    endif
    if( cidp3 /= 0 ) then
      m3 = 1.0d0
    endif
    if( cidp4 /= 0 ) then
      m4 = 1.0d0
    endif
    if( cidp5 /= 0 ) then
      m5 = 1.0d0
    endif

    up = ux_(i, j, k)
    uw = ux_(i-1, j, k)*(1.0 - m0) + Uc*m0
    ue = ux_(i+1, j, k)*(1.0 - m1) + Uc*m1
    us = ux_(i, j-1, k)*(1.0 - m2) + Uc*m2
    un = ux_(i, j+1, k)*(1.0 - m3) + Uc*m3
    ub = ux_(i, j, k-1)*(1.0 - m4) + Uc*m4
    ut = ux_(i, j, k+1)*(1.0 - m5) + Uc*m5

    vp = uy_(i, j, k)
    vw = uy_(i-1, j, k)*(1.0 - m0) + Vc*m0
    ve = uy_(i+1, j, k)*(1.0 - m1) + Vc*m1
    vs = uy_(i, j-1, k)*(1.0 - m2) + Vc*m2
    vn = uy_(i, j+1, k)*(1.0 - m3) + Vc*m3
    vb = uy_(i, j, k-1)*(1.0 - m4) + Vc*m4
    vt = uy_(i, j, k+1)*(1.0 - m5) + Vc*m5

    wp = uz_(i, j, k)
    ww = uz_(i-1, j, k)*(1.0 - m0) + Wc*m0
    we = uz_(i+1, j, k)*(1.0 - m1) + Wc*m1
    ws = uz_(i, j-1, k)*(1.0 - m2) + Wc*m2
    wn = uz_(i, j+1, k)*(1.0 - m3) + Wc*m3
    wb = uz_(i, j, k-1)*(1.0 - m4) + Wc*m4
    wt = uz_(i, j, k+1)*(1.0 - m5) + Wc*m5

		dudx = 2.0/(d1 + d0)*((ue - up)/(d1*dx) + (up - uw)/(d0*dx))
		dudy = 2.0/(d3 + d2)*((un - up)/(d3*dx) + (up - us)/(d2*dx))
		dudz = 2.0/(d5 + d4)*((ut - up)/(d5*dx) + (up - ub)/(d4*dx))

		dvdx = 2.0/(d1 + d0)*((ve - vp)/(d1*dx) + (vp - vw)/(d0*dx))
		dvdy = 2.0/(d3 + d2)*((vn - vp)/(d3*dx) + (vp - vs)/(d2*dx))
		dvdz = 2.0/(d5 + d4)*((vt - vp)/(d5*dx) + (vp - vb)/(d4*dx))

		dwdx = 2.0/(d1 + d0)*((we - wp)/(d1*dx) + (wp - ww)/(d0*dx))
		dwdy = 2.0/(d3 + d2)*((wn - wp)/(d3*dx) + (wp - ws)/(d2*dx))
		dwdz = 2.0/(d5 + d4)*((wt - wp)/(d5*dx) + (wp - wb)/(d4*dx))

		sxx = 0.5*(dudx + dudx)
		sxy = 0.5*(dvdx + dudy)
		sxz = 0.5*(dwdx + dudz)
		syx = 0.5*(dudy + dvdx)
		syy = 0.5*(dvdy + dvdy)
		syz = 0.5*(dwdy + dvdz)
		szx = 0.5*(dudz + dwdx)
		szy = 0.5*(dvdz + dwdy)
		szz = 0.5*(dwdz + dwdz)
		s2 = 2.0*(sxx*sxx + syy*syy + szz*szz) + 4.0*(sxy*sxy + syz*syz + szx*szx)
		sl = sqrt(s2)

!		wxx = 0.5*(dudx - dudx)
!		wxy = 0.5*(dvdx - dudy)
!		wxz = 0.5*(dwdx - dudz)
!		wyx = 0.5*(dudy - dvdx)
!		wyy = 0.5*(dvdy - dvdy)
!		wyz = 0.5*(dwdy - dvdz)
!		wzx = 0.5*(dudz - dwdx)
!		wzy = 0.5*(dvdz - dwdy)
!		wzz = 0.5*(dwdz - dwdz)

		csm_e =  0.5*(dudx*dudx + dudy*dudy + dudz*dudz &
								+ dvdx*dvdx + dvdy*dvdy + dvdz*dvdz &
								+ dwdx*dwdx + dwdy*dwdy + dwdz*dwdz)
		csm_q = -0.5*(dudx*dudx + dudy*dvdx + dudz*dwdx &
								+ dvdx*dudy + dvdy*dvdy + dvdz*dwdy &
								+ dwdx*dudz + dwdy*dvdz + dwdz*dwdz)
		csm_fcs = csm_q/csm_e
		csm_fom = 1.0 - csm_fcs
		csm_c = csm_c2*csm_fom*(abs(csm_fcs))**1.5

		nue_(i, j, k) = csm_c*dx*dx*sl
    if( pidp /= 1 ) then
			nue_(i, j, k) = 0.0
    endif
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_calc_nue

subroutine bcut_calc_d_t( &
                td0_, &
                t0_, &
                c0, c1, c2, c3, c4, c5, &
                cid0, cid1, cid2, cid3, cid4, cid5, &
                pid, &
								nnum, &
								nx_, ny_, nz_, &
								nidx0, nidx1, nidx2, nidx3, nidx4, nidx5, &
                rhof, rhos, &
                cpf, cps, &
                kf, ks, &
                bc_n, &
                bc_type, &
                bc_value, &
                org, &
                dx, dt, &
                sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: td0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: t0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
	integer									 :: nnum
	real, dimension(0:nnum-1) :: nx_, ny_, nz_
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: nidx0, nidx1, nidx2, nidx3, nidx4, nidx5
  integer                  :: cidp
  integer                  :: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
  integer                  :: pidp, pidw, pide, pids, pidn, pidb, pidt
  real                    :: rhof, rhos
  real                    :: cpf, cps
  real                    :: kf, ks
  integer                      :: bc_n
  integer, dimension(0:bc_n-1):: bc_type
  real, dimension(0:bc_n-1)   :: bc_value
  real                    :: dx, dt
  real, dimension(3)      :: org
  real                    :: d0, d1, d2, d3, d4, d5
  real                    :: k0, k1, k2, k3, k4, k5
	real										:: k_p, k_w, k_e, k_s, k_n, k_b, k_t
	real										:: rho_p
	real										:: cp_p
  real                    :: m0, m1, m2, m3, m4, m5
  real                    :: l0, l1, l2, l3, l4, l5
  real                    :: bp, b0, b1, b2, b3, b4, b5
  real                    :: tp, t0, t1, t2, t3, t4, t5
  real                    :: nx, ny, nz
  real                    :: x, y, z, r2, rl, xi, yi, zi, theta, phi
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp           private(cidp) &
!$omp           private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp           private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp           private(d0, d1, d2, d3, d4, d5) &
!$omp           private(k0, k1, k2, k3, k4, k5) &
!$omp						private(k_p, k_w, k_e, k_s, k_n, k_b, k_t) &
!$omp						private(rho_p, cp_p) &
!$omp           private(m0, m1, m2, m3, m4, m5) &
!$omp           private(l0, l1, l2, l3, l4, l5) &
!$omp           private(bp, b0, b1, b2, b3, b4, b5) &
!$omp           private(tp, t0, t1, t2, t3, t4, t5) &
!$omp           private(nx, ny, nz) &
!$omp           private(x, y, z, r2, rl, xi, yi, zi, theta, phi)
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
    x = org(1) + (real(i) - 0.5)*dx
    y = org(2) + (real(j) - 0.5)*dx
    z = org(3) + (real(k) - 0.5)*dx

    d0 = c0(i, j, k)
    d1 = c1(i, j, k)
    d2 = c2(i, j, k)
    d3 = c3(i, j, k)
    d4 = c4(i, j, k)
    d5 = c5(i, j, k)

    cidp0 = cid0(i, j, k)
    cidp1 = cid1(i, j, k)
    cidp2 = cid2(i, j, k)
    cidp3 = cid3(i, j, k)
    cidp4 = cid4(i, j, k)
    cidp5 = cid5(i, j, k)

    pidp = pid(i, j, k)

    m0 = 0.0
    m1 = 0.0
    m2 = 0.0
    m3 = 0.0
    m4 = 0.0
    m5 = 0.0
    if( cidp0 /= 0 ) then
      m0 = 1.0
    endif
    if( cidp1 /= 0 ) then
      m1 = 1.0
    endif
    if( cidp2 /= 0 ) then
      m2 = 1.0
    endif
    if( cidp3 /= 0 ) then
      m3 = 1.0
    endif
    if( cidp4 /= 0 ) then
      m4 = 1.0
    endif
    if( cidp5 /= 0 ) then
      m5 = 1.0
    endif

		rho_p = rhof
		cp_p  = cpf
		k_p   = kf
		k_w   = kf
		k_e   = kf
		k_s   = kf
		k_n   = kf
		k_b   = kf
		k_t   = kf
		if( pid(i, j, k) /= 0 ) then
			rho_p = rhos
			cp_p  = cps
			k_p   = ks
		end if
		if( pid(i-1, j, k) /= 0 ) then
			k_w   = ks
		end if
		if( pid(i+1, j, k) /= 0 ) then
			k_e   = ks
		end if
		if( pid(i, j-1, k) /= 0 ) then
			k_s   = ks
		end if
		if( pid(i, j+1, k) /= 0 ) then
			k_n   = ks
		end if
		if( pid(i, j, k-1) /= 0 ) then
			k_b   = ks
		end if
		if( pid(i, j, k+1) /= 0 ) then
			k_t   = ks
		end if

    k0 = kf
    k1 = kf
    k2 = kf
    k3 = kf
    k4 = kf
    k5 = kf
    tp = t0_(i, j, k)
    t0 = t0_(i-1, j, k) 
    t1 = t0_(i+1, j, k) 
    t2 = t0_(i, j-1, k) 
    t3 = t0_(i, j+1, k) 
    t4 = t0_(i, j, k-1) 
    t5 = t0_(i, j, k+1) 
    b0 = 0.0
    b1 = 0.0
    b2 = 0.0
    b3 = 0.0
    b4 = 0.0
    b5 = 0.0
    if( bc_type(cidp0) == 0 ) then
      k0 = kf/d0*2.0/(d0 + d1)
      k1 = kf/d1*2.0/(d0 + d1)
      t0 = bc_value(cidp0)
			b0 = 0.0
    else if( bc_type(cidp0) == 1 ) then
			nx = nx_( nidx0(i, j, k) )

      k0 = 0.0
      k1 = kf/(d0 + 0.5)
      t0 = 0.0
      b0 = - nx*bc_value(cidp0)*dx/(d0 + 0.5)
    else if( bc_type(cidp0) == 2 ) then
      k0 = k_p*k_w/(k_p*(1.0 - d0) + k_w*d0)
      k1 = k_p*k_e/(k_p*(1.0 - d1) + k_e*d1)
      t0 = 0.0
      b0 = 0.0
    endif
    if( bc_type(cidp1) == 0 ) then
      k0 = kf/d0*2.0/(d0 + d1)
      k1 = kf/d1*2.0/(d0 + d1)
      t1 = bc_value(cidp1)
    else if( bc_type(cidp1) == 1 ) then
			nx = nx_( nidx1(i, j, k) )

      k0 = kf/(d1 + 0.5)
      k1 = 0.0
      t1 = 0.0
      b1 = + nx*bc_value(cidp1)*dx/(d1 + 0.5)
    else if( bc_type(cidp1) == 2 ) then
      k0 = k_p*k_w/(k_p*(1.0 - d0) + k_w*d0)
      k1 = k_p*k_e/(k_p*(1.0 - d1) + k_e*d1)
      t1 = 0.0
      b1 = 0.0
    endif
    if( bc_type(cidp2) == 0 ) then
      k2 = kf/d2*2.0/(d2 + d3)
      k3 = kf/d3*2.0/(d2 + d3)
      t2 = bc_value(cidp2)
    else if( bc_type(cidp2) == 1 ) then
			ny = ny_( nidx2(i, j, k) )

      k2 = 0.0
      k3 = kf/(d2 + 0.5)
      t2 = 0.0
      b2 = - ny*bc_value(cidp2)*dx/(d2 + 0.5)
    else if( bc_type(cidp2) == 2 ) then
      k2 = k_p*k_s/(k_p*(1.0 - d2) + k_s*d2)
      k3 = k_p*k_n/(k_p*(1.0 - d3) + k_n*d3)
      t2 = 0.0
      b2 = 0.0
    endif
    if( bc_type(cidp3) == 0 ) then
      k2 = kf/d2*2.0/(d2 + d3)
      k3 = kf/d3*2.0/(d2 + d3)
      t3 = bc_value(cidp3)
    else if( bc_type(cidp3) == 1 ) then
			ny = ny_( nidx3(i, j, k) )

      k2 = kf/(d3 + 0.5)
      k3 = 0.0
      t3 = 0.0
      b3 = + ny*bc_value(cidp3)*dx/(d3 + 0.5)
    else if( bc_type(cidp3) == 2 ) then
      k2 = k_p*k_s/(k_p*(1.0 - d2) + k_s*d2)
      k3 = k_p*k_n/(k_p*(1.0 - d3) + k_n*d3)
      t3 = 0.0
      b3 = 0.0
    endif
    if( bc_type(cidp4) == 0 ) then
      k4 = kf/d4*2.0/(d4 + d5)
      k5 = kf/d5*2.0/(d4 + d5)
      t4 = bc_value(cidp4)
    else if( bc_type(cidp4) == 1 ) then
			nz = nz_( nidx4(i, j, k) )

      k4 = 0.0
      k5 = kf/(d4 + 0.5)
      t4 = 0.0
      b4 = - nz*bc_value(cidp4)*dx/(d4 + 0.5)
    else if( bc_type(cidp4) == 2 ) then
      k4 = k_p*k_b/(k_p*(1.0 - d4) + k_b*d4)
      k5 = k_p*k_t/(k_p*(1.0 - d5) + k_t*d5)
      t4 = 0.0
      b4 = 0.0
    endif
    if( bc_type(cidp5) == 0 ) then
      k4 = kf/d4*2.0/(d4 + d5)
      k5 = kf/d5*2.0/(d4 + d5)
      t5 = bc_value(cidp5)
    else if( bc_type(cidp5) == 1 ) then
			nz = nz_( nidx5(i, j, k) )

      k4 = kf/(d5 + 0.5)
      k5 = 0.0
      t5 = 0.0
      b5 = + nz*bc_value(cidp5)*dx/(d5 + 0.5)
    else if( bc_type(cidp5) == 2 ) then
      k4 = k_p*k_b/(k_p*(1.0 - d4) + k_b*d4)
      k5 = k_p*k_t/(k_p*(1.0 - d5) + k_t*d5)
      t5 = 0.0
      b5 = 0.0
    endif

    l0 = k0/(rho_p*cp_p)/(dx*dx)*dt
    l1 = k1/(rho_p*cp_p)/(dx*dx)*dt
    l2 = k2/(rho_p*cp_p)/(dx*dx)*dt
    l3 = k3/(rho_p*cp_p)/(dx*dx)*dt
    l4 = k4/(rho_p*cp_p)/(dx*dx)*dt
    l5 = k5/(rho_p*cp_p)/(dx*dx)*dt

    td0_(i, j, k) = ( k1*(t1 - tp) &
                    - k0*(tp - t0) &
                    + k3*(t3 - tp) &
                    - k2*(tp - t2) &
                    + k5*(t5 - tp) &
                    - k4*(tp - t4) &
                    )/(rhof*cpf)/(dx*dx)

    if( pidp /= 1 ) then
      td0_(i, j, k) = 0.0
    endif

  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_calc_d_t

subroutine bcut_calc_ab_u( &
                Ap, Aw, Ae, As, An, Ab, At, b, &
                u0_, &
                uc0_, ucp_, &
                ud0_, &
                p0_, &
                c0, c1, c2, c3, c4, c5, &
                cid0, cid1, cid2, cid3, cid4, cid5, &
                pid, &
                dir, &
                rhof, &
                mu, &
                dx, dt, &
                Uc, &
								eps, &
                gx, gy, gz, &
                sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: Ap, Aw, Ae, As, An, Ab, At, b
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: u0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: uc0_, ucp_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: ud0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: p0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
  integer                  :: cidp
  integer                  :: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
  integer                  :: pidp, pidw, pide, pids, pidn, pidb, pidt
  integer                  :: dir
  real                    :: rhof
  real                    :: mu
  real                    :: dx, dt
  real                    :: Uc
	real										:: eps
  real                    :: gx, gy, gz
  real                    :: d0, d1, d2, d3, d4, d5
  real                    :: mu0, mu1, mu2, mu3, mu4, mu5
  real                    :: m0, m1, m2, m3, m4, m5
  real                    :: l0, l1, l2, l3, l4, l5
  real                    :: r0, r1, r2, r3, r4, r5
  real                    :: dpx0, dpx1, dpy2, dpy3, dpz4, dpz5
  real                    :: rdpx0, rdpx1, rdpy2, rdpy3, rdpz4, rdpz5
  real                    :: rdpx, rdpy, rdpz, rdp
  real                    :: u0
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp           private(cidp) &
!$omp           private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp           private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp           private(d0, d1, d2, d3, d4, d5) &
!$omp           private(mu0, mu1, mu2, mu3, mu4, mu5) &
!$omp           private(m0, m1, m2, m3, m4, m5) &
!$omp           private(l0, l1, l2, l3, l4, l5) &
!$omp           private(r0, r1, r2, r3, r4, r5) &
!$omp           private(dpx0, dpx1, dpy2, dpy3, dpz4, dpz5) &
!$omp           private(rdpx0, rdpx1, rdpy2, rdpy3, rdpz4, rdpz5) &
!$omp           private(rdpx, rdpy, rdpz, rdp) &
!$omp           private(u0) 
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
    d0 = c0(i, j, k)
    d1 = c1(i, j, k)
    d2 = c2(i, j, k)
    d3 = c3(i, j, k)
    d4 = c4(i, j, k)
    d5 = c5(i, j, k)

    cidp0 = cid0(i, j, k)
    cidp1 = cid1(i, j, k)
    cidp2 = cid2(i, j, k)
    cidp3 = cid3(i, j, k)
    cidp4 = cid4(i, j, k)
    cidp5 = cid5(i, j, k)

    pidp = pid(i, j, k)

    mu0 = mu
    mu1 = mu
    mu2 = mu
    mu3 = mu
    mu4 = mu
    mu5 = mu

    r0 = 1.0d0/rhof
    r1 = 1.0d0/rhof
    r2 = 1.0d0/rhof
    r3 = 1.0d0/rhof
    r4 = 1.0d0/rhof
    r5 = 1.0d0/rhof

    dpx0 = (p0_(i, j, k) - p0_(i-1, j, k))/dx
    dpx1 = (p0_(i+1, j, k) - p0_(i, j, k))/dx
    dpy2 = (p0_(i, j, k) - p0_(i, j-1, k))/dx
    dpy3 = (p0_(i, j+1, k) - p0_(i, j, k))/dx
    dpz4 = (p0_(i, j, k) - p0_(i, j, k-1))/dx
    dpz5 = (p0_(i, j, k+1) - p0_(i, j, k))/dx

    rdpx0 = r0*dpx0 - gx*0.0
    rdpx1 = r1*dpx1 - gx*0.0
    rdpy2 = r2*dpy2 - gy*0.0
    rdpy3 = r3*dpy3 - gy*0.0
    rdpz4 = r4*dpz4 - gz*0.0
    rdpz5 = r5*dpz5 - gz*0.0

    m0 = 0.0d0
    m1 = 0.0d0
    m2 = 0.0d0
    m3 = 0.0d0
    m4 = 0.0d0
    m5 = 0.0d0

    if( cidp0 /= 0 ) then
      mu0 = mu/d0*2.0d0/(d0 + d1)
      mu1 = mu/d1*2.0d0/(d0 + d1)
			if( d0 < eps ) then
				mu0 = 0.0
				mu1 = mu
			end if
      rdpx0 = rdpx1*(d0 - 0.5d0)/(d0 + 0.5d0)
      dpx0 = dpx1*(d0 - 0.5d0)/(d0 + 0.5d0)
      m0 = 1.0d0
    endif

    if( cidp1 /= 0 ) then
      mu0 = mu/d0*2.0d0/(d0 + d1)
      mu1 = mu/d1*2.0d0/(d0 + d1)
			if( d1 < eps ) then
				mu0 = mu
				mu1 = 0.0
			end if
      rdpx1 = rdpx0*(d1 - 0.5d0)/(d1 + 0.5d0)
      dpx1 = dpx0*(d1 - 0.5d0)/(d1 + 0.5d0)
      m1 = 1.0d0
    endif

    if( cidp2 /= 0 ) then
      mu2 = mu/d2*2.0d0/(d2 + d3)
      mu3 = mu/d3*2.0d0/(d2 + d3)
			if( d2 < eps ) then
				mu2 = 0.0
				mu3 = mu
			end if
      rdpy2 = rdpy3*(d2 - 0.5d0)/(d2 + 0.5d0)
      dpy2 = dpy3*(d2 - 0.5d0)/(d2 + 0.5d0)
      m2 = 1.0d0
    endif

    if( cidp3 /= 0 ) then
      mu2 = mu/d2*2.0d0/(d2 + d3)
      mu3 = mu/d3*2.0d0/(d2 + d3)
			if( d3 < eps ) then
				mu2 = mu
				mu3 = 0.0
			end if
      rdpy3 = rdpy2*(d3 - 0.5d0)/(d3 + 0.5d0)
      dpy3 = dpy2*(d3 - 0.5d0)/(d3 + 0.5d0)
      m3 = 1.0d0
    endif

    if( cidp4 /= 0 ) then
      mu4 = mu/d4*2.0d0/(d4 + d5)
      mu5 = mu/d5*2.0d0/(d4 + d5)
			if( d4 < eps ) then
				mu4 = 0.0
				mu5 = mu
			end if
      rdpz4 = rdpz5*(d4 - 0.5d0)/(d4 + 0.5d0)
      dpz4 = dpz5*(d4 - 0.5d0)/(d4 + 0.5d0)
      m4 = 1.0d0
    endif

    if( cidp5 /= 0 ) then
      mu4 = mu/d4*2.0d0/(d4 + d5)
      mu5 = mu/d5*2.0d0/(d4 + d5)
			if( d5 < eps ) then
				mu4 = mu
				mu5 = 0.0
			end if
      rdpz5 = rdpz4*(d5 - 0.5d0)/(d5 + 0.5d0)
      dpz5 = dpz4*(d5 - 0.5d0)/(d5 + 0.5d0)
      m5 = 1.0d0
    endif

    if( cidp0 /= 0 .and. cidp1 /= 0 ) then
      rdpx0 = 0.0
      rdpx1 = 0.0
      dpx0 = 0.0
      dpx1 = 0.0
      m0 = 1.0d0
      m1 = 1.0d0
    endif

    if( cidp2 /= 0 .and. cidp3 /= 0 ) then
      rdpy2 = 0.0
      rdpy3 = 0.0
      dpy2 = 0.0
      dpy3 = 0.0
      m2 = 1.0d0
      m3 = 1.0d0
    endif

    if( cidp4 /= 0 .and. cidp5 /= 0 ) then
      rdpz4 = 0.0
      rdpz5 = 0.0
      dpz4 = 0.0
      dpz5 = 0.0
      m4 = 1.0d0
      m5 = 1.0d0
    endif

    rdpx = 0.5d0*(rdpx0 + rdpx1)
    rdpy = 0.5d0*(rdpy2 + rdpy3)
    rdpz = 0.5d0*(rdpz4 + rdpz5)

    rdp = 0.0
    if( dir==0 ) then
      rdp = rdpx 
    else if( dir==1 ) then
      rdp = rdpy
    else if( dir==2 ) then
      rdp = rdpz
    end if

    l0 = mu0/(rhof)*dt/(dx*dx)
    l1 = mu1/(rhof)*dt/(dx*dx)
    l2 = mu2/(rhof)*dt/(dx*dx)
    l3 = mu3/(rhof)*dt/(dx*dx)
    l4 = mu4/(rhof)*dt/(dx*dx)
    l5 = mu5/(rhof)*dt/(dx*dx)

    u0 = u0_(i, j, k)

    Ap(i, j, k) = 1.0d0 + 0.5d0*(l0 + l1) &
                        + 0.5d0*(l2 + l3) &
                        + 0.5d0*(l4 + l5)
    Aw(i, j, k) = -0.5d0*l0*(1.0d0 - m0)
    Ae(i, j, k) = -0.5d0*l1*(1.0d0 - m1)
    As(i, j, k) = -0.5d0*l2*(1.0d0 - m2)
    An(i, j, k) = -0.5d0*l3*(1.0d0 - m3)
    Ab(i, j, k) = -0.5d0*l4*(1.0d0 - m4)
    At(i, j, k) = -0.5d0*l5*(1.0d0 - m5)
     b(i, j, k) = u0 &
                  - (1.5d0*uc0_(i, j, k) - 0.5d0*ucp_(i, j, k) )*dt &
                  + 0.5d0*ud0_(i, j, k)*dt &
                  - rdp*dt &
                  + 0.5d0*l0*Uc*m0 &
                  + 0.5d0*l1*Uc*m1 &
                  + 0.5d0*l2*Uc*m2 &
                  + 0.5d0*l3*Uc*m3 &
                  + 0.5d0*l4*Uc*m4 &
                  + 0.5d0*l5*Uc*m5 

    if( pidp /= 1 ) then
      Ap(i, j, k) = 1.0d0
      Aw(i, j, k) = 0.0d0
      Ae(i, j, k) = 0.0d0
      As(i, j, k) = 0.0d0
      An(i, j, k) = 0.0d0
      Ab(i, j, k) = 0.0d0
      At(i, j, k) = 0.0d0
      b (i, j, k) = Uc
      u0_(i, j, k) = 0.0d0
    endif

  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_calc_ab_u

subroutine bcut_update_u( &
                u0_, &
                uc0_, ucp_, &
                ud0_, udp_, &
                p0_, &
                c0, c1, c2, c3, c4, c5, &
                cid0, cid1, cid2, cid3, cid4, cid5, &
                pid, &
                dir, &
                rhof, &
                mu, &
                dx, dt, &
                Uc, &
                gx, gy, gz, &
                sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: u0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: uc0_, ucp_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: ud0_, udp_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: p0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
  integer                  :: cidp
  integer                  :: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
  integer                  :: pidp, pidw, pide, pids, pidn, pidb, pidt
  integer                  :: dir
  real                    :: rhof
  real                    :: mu
  real                    :: dx, dt
  real                    :: Uc
  real                    :: gx, gy, gz
  real                    :: d0, d1, d2, d3, d4, d5
  real                    :: mu0, mu1, mu2, mu3, mu4, mu5
  real                    :: m0, m1, m2, m3, m4, m5
  real                    :: l0, l1, l2, l3, l4, l5
  real                    :: r0, r1, r2, r3, r4, r5
  real                    :: dpx0, dpx1, dpy2, dpy3, dpz4, dpz5
  real                    :: rdpx0, rdpx1, rdpy2, rdpy3, rdpz4, rdpz5
  real                    :: rdpx, rdpy, rdpz, rdp
  real                    :: u0
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp           private(cidp) &
!$omp           private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp           private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp           private(d0, d1, d2, d3, d4, d5) &
!$omp           private(mu0, mu1, mu2, mu3, mu4, mu5) &
!$omp           private(m0, m1, m2, m3, m4, m5) &
!$omp           private(l0, l1, l2, l3, l4, l5) &
!$omp           private(r0, r1, r2, r3, r4, r5) &
!$omp           private(dpx0, dpx1, dpy2, dpy3, dpz4, dpz5) &
!$omp           private(rdpx0, rdpx1, rdpy2, rdpy3, rdpz4, rdpz5) &
!$omp           private(rdpx, rdpy, rdpz, rdp) &
!$omp           private(u0) 
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
    d0 = c0(i, j, k)
    d1 = c1(i, j, k)
    d2 = c2(i, j, k)
    d3 = c3(i, j, k)
    d4 = c4(i, j, k)
    d5 = c5(i, j, k)

    cidp0 = cid0(i, j, k)
    cidp1 = cid1(i, j, k)
    cidp2 = cid2(i, j, k)
    cidp3 = cid3(i, j, k)
    cidp4 = cid4(i, j, k)
    cidp5 = cid5(i, j, k)

    pidp = pid(i, j, k)

    mu0 = mu
    mu1 = mu
    mu2 = mu
    mu3 = mu
    mu4 = mu
    mu5 = mu

    r0 = 1.0d0/rhof
    r1 = 1.0d0/rhof
    r2 = 1.0d0/rhof
    r3 = 1.0d0/rhof
    r4 = 1.0d0/rhof
    r5 = 1.0d0/rhof

    dpx0 = (p0_(i, j, k) - p0_(i-1, j, k))/dx
    dpx1 = (p0_(i+1, j, k) - p0_(i, j, k))/dx
    dpy2 = (p0_(i, j, k) - p0_(i, j-1, k))/dx
    dpy3 = (p0_(i, j+1, k) - p0_(i, j, k))/dx
    dpz4 = (p0_(i, j, k) - p0_(i, j, k-1))/dx
    dpz5 = (p0_(i, j, k+1) - p0_(i, j, k))/dx

    rdpx0 = r0*dpx0 - gx*0.0
    rdpx1 = r1*dpx1 - gx*0.0
    rdpy2 = r2*dpy2 - gy*0.0
    rdpy3 = r3*dpy3 - gy*0.0
    rdpz4 = r4*dpz4 - gz*0.0
    rdpz5 = r5*dpz5 - gz*0.0

    m0 = 0.0d0
    m1 = 0.0d0
    m2 = 0.0d0
    m3 = 0.0d0
    m4 = 0.0d0
    m5 = 0.0d0

    if( cidp0 /= 0 ) then
      mu0 = mu/d0*2.0d0/(d0 + d1)
      mu1 = mu/d1*2.0d0/(d0 + d1)
      rdpx0 = rdpx1*(d0 - 0.5d0)/(d0 + 0.5d0)
      dpx0 = dpx1*(d0 - 0.5d0)/(d0 + 0.5d0)
      m0 = 1.0d0
    endif

    if( cidp1 /= 0 ) then
      mu0 = mu/d0*2.0d0/(d0 + d1)
      mu1 = mu/d1*2.0d0/(d0 + d1)
      rdpx1 = rdpx0*(d1 - 0.5d0)/(d1 + 0.5d0)
      dpx1 = dpx0*(d1 - 0.5d0)/(d1 + 0.5d0)
      m1 = 1.0d0
    endif

    if( cidp2 /= 0 ) then
      mu2 = mu/d2*2.0d0/(d2 + d3)
      mu3 = mu/d3*2.0d0/(d2 + d3)
      rdpy2 = rdpy3*(d2 - 0.5d0)/(d2 + 0.5d0)
      dpy2 = dpy3*(d2 - 0.5d0)/(d2 + 0.5d0)
      m2 = 1.0d0
    endif

    if( cidp3 /= 0 ) then
      mu2 = mu/d2*2.0d0/(d2 + d3)
      mu3 = mu/d3*2.0d0/(d2 + d3)
      rdpy3 = rdpy2*(d3 - 0.5d0)/(d3 + 0.5d0)
      dpy3 = dpy2*(d3 - 0.5d0)/(d3 + 0.5d0)
      m3 = 1.0d0
    endif

    if( cidp4 /= 0 ) then
      mu4 = mu/d4*2.0d0/(d4 + d5)
      mu5 = mu/d5*2.0d0/(d4 + d5)
      rdpz4 = rdpz5*(d4 - 0.5d0)/(d4 + 0.5d0)
      dpz4 = dpz5*(d4 - 0.5d0)/(d4 + 0.5d0)
      m4 = 1.0d0
    endif

    if( cidp5 /= 0 ) then
      mu4 = mu/d4*2.0d0/(d4 + d5)
      mu5 = mu/d5*2.0d0/(d4 + d5)
      rdpz5 = rdpz4*(d5 - 0.5d0)/(d5 + 0.5d0)
      dpz5 = dpz4*(d5 - 0.5d0)/(d5 + 0.5d0)
      m5 = 1.0d0
    endif

    if( cidp0 /= 0 .and. cidp1 /= 0 ) then
      rdpx0 = 0.0
      rdpx1 = 0.0
      dpx0 = 0.0
      dpx1 = 0.0
      m0 = 1.0d0
      m1 = 1.0d0
    endif

    if( cidp2 /= 0 .and. cidp3 /= 0 ) then
      rdpy2 = 0.0
      rdpy3 = 0.0
      dpy2 = 0.0
      dpy3 = 0.0
      m2 = 1.0d0
      m3 = 1.0d0
    endif

    if( cidp4 /= 0 .and. cidp5 /= 0 ) then
      rdpz4 = 0.0
      rdpz5 = 0.0
      dpz4 = 0.0
      dpz5 = 0.0
      m4 = 1.0d0
      m5 = 1.0d0
    endif

    rdpx = 0.5d0*(rdpx0 + rdpx1)
    rdpy = 0.5d0*(rdpy2 + rdpy3)
    rdpz = 0.5d0*(rdpz4 + rdpz5)

    rdp = 0.0
    if( dir==0 ) then
      rdp = rdpx 
    else if( dir==1 ) then
      rdp = rdpy
    else if( dir==2 ) then
      rdp = rdpz
    end if

    u0_(i, j, k) = u0_(i, j, k) &
                  - (1.5*uc0_(i, j, k) - 0.5*ucp_(i, j, k))*dt &
                  + (1.5*ud0_(i, j, k) - 0.5*udp_(i, j, k))*dt &
                  - rdp*dt

    if( pidp /= 1 ) then
      u0_(i, j, k) = Uc
    endif
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_update_u

subroutine bcut_add_g( &
                ux_, uy_, uz_, &
								t_, &
                gx, gy, gz, &
								betag, tr, &
								dx, dt, &
                sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: ux_, uy_, uz_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: t_
	real										:: gx, gy, gz
	real										:: betag, tr
	real										:: dx, dt
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) 
!$omp do
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
		ux_(i, j, k) = ux_(i, j, k) + gx*dt
		uy_(i, j, k) = uy_(i, j, k) + gy*dt
		uz_(i, j, k) = uz_(i, j, k) + gz*dt + betag*(t_(i, j, k) - tr)*dt
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_add_g

subroutine bcut_remove_p( &
                ux, uy, uz, &
                p0_, &
                c0, c1, c2, c3, c4, c5, &
                cid0, cid1, cid2, cid3, cid4, cid5, &
                pid, &
                rhof, &
                dx, dt, &
                sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: ux, uy, uz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: p0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
  integer                  :: cidp
  integer                  :: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
  integer                  :: pidp, pidw, pide, pids, pidn, pidb, pidt
  real                    :: rhof
  real                    :: dx, dt
  real                    :: d0, d1, d2, d3, d4, d5
  real                    :: r0, r1, r2, r3, r4, r5
  real                    :: m0, m1, m2, m3, m4, m5
  real                    :: l0, l1, l2, l3, l4, l5
  real                    :: dpx0, dpx1, dpy2, dpy3, dpz4, dpz5
  real                    :: rdpx0, rdpx1, rdpy2, rdpy3, rdpz4, rdpz5
  real                    :: rdpx, rdpy, rdpz, rdp
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp           private(cidp) &
!$omp           private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp           private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp           private(d0, d1, d2, d3, d4, d5) &
!$omp           private(r0, r1, r2, r3, r4, r5) &
!$omp           private(m0, m1, m2, m3, m4, m5) &
!$omp           private(l0, l1, l2, l3, l4, l5) &
!$omp           private(dpx0, dpx1, dpy2, dpy3, dpz4, dpz5) &
!$omp           private(rdpx0, rdpx1, rdpy2, rdpy3, rdpz4, rdpz5) &
!$omp           private(rdpx, rdpy, rdpz, rdp)
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
    r0 = 1.0d0/rhof
    r1 = 1.0d0/rhof
    r2 = 1.0d0/rhof
    r3 = 1.0d0/rhof
    r4 = 1.0d0/rhof
    r5 = 1.0d0/rhof

    dpx0 = (p0_(i, j, k) - p0_(i-1, j, k))/dx
    dpx1 = (p0_(i+1, j, k) - p0_(i, j, k))/dx
    dpy2 = (p0_(i, j, k) - p0_(i, j-1, k))/dx
    dpy3 = (p0_(i, j+1, k) - p0_(i, j, k))/dx
    dpz4 = (p0_(i, j, k) - p0_(i, j, k-1))/dx
    dpz5 = (p0_(i, j, k+1) - p0_(i, j, k))/dx

    rdpx0 = r0*dpx0
    rdpx1 = r1*dpx1
    rdpy2 = r2*dpy2
    rdpy3 = r3*dpy3
    rdpz4 = r4*dpz4
    rdpz5 = r5*dpz5

    m0 = 0.0d0
    m1 = 0.0d0
    m2 = 0.0d0
    m3 = 0.0d0
    m4 = 0.0d0
    m5 = 0.0d0

    d0 = c0(i, j, k)
    d1 = c1(i, j, k)
    d2 = c2(i, j, k)
    d3 = c3(i, j, k)
    d4 = c4(i, j, k)
    d5 = c5(i, j, k)

    cidp0 = cid0(i, j, k)
    cidp1 = cid1(i, j, k)
    cidp2 = cid2(i, j, k)
    cidp3 = cid3(i, j, k)
    cidp4 = cid4(i, j, k)
    cidp5 = cid5(i, j, k)

    pidp = pid(i, j, k)

    if( cidp0 /= 0 ) then
      rdpx0 = rdpx1*(d0 - 0.5d0)/(d0 + 0.5d0)
      dpx0 = dpx1*(d0 - 0.5d0)/(d0 + 0.5d0)
      m0 = 1.0d0
    endif

    if( cidp1 /= 0 ) then
      rdpx1 = rdpx0*(d1 - 0.5d0)/(d1 + 0.5d0)
      dpx1 = dpx0*(d1 - 0.5d0)/(d1 + 0.5d0)
      m1 = 1.0d0
    endif

    if( cidp2 /= 0 ) then
      rdpy2 = rdpy3*(d2 - 0.5d0)/(d2 + 0.5d0)
      dpy2 = dpy3*(d2 - 0.5d0)/(d2 + 0.5d0)
      m2 = 1.0d0
    endif

    if( cidp3 /= 0 ) then
      rdpy3 = rdpy2*(d3 - 0.5d0)/(d3 + 0.5d0)
      dpy3 = dpy2*(d3 - 0.5d0)/(d3 + 0.5d0)
      m3 = 1.0d0
    endif

    if( cidp4 /= 0 ) then
      rdpz4 = rdpz5*(d4 - 0.5d0)/(d4 + 0.5d0)
      dpz4 = dpz5*(d4 - 0.5d0)/(d4 + 0.5d0)
      m4 = 1.0d0
    endif

    if( cidp5 /= 0 ) then
      rdpz5 = rdpz4*(d5 - 0.5d0)/(d5 + 0.5d0)
      dpz5 = dpz4*(d5 - 0.5d0)/(d5 + 0.5d0)
      m5 = 1.0d0
    endif

    if( cidp0 /= 0 .and. cidp1 /= 0 ) then
      rdpx0 = 0.0
      rdpx1 = 0.0
      dpx0 = 0.0
      dpx1 = 0.0
      m0 = 1.0d0
      m1 = 1.0d0
    endif

    if( cidp2 /= 0 .and. cidp3 /= 0 ) then
      rdpy2 = 0.0
      rdpy3 = 0.0
      dpy2 = 0.0
      dpy3 = 0.0
      m2 = 1.0d0
      m3 = 1.0d0
    endif

    if( cidp4 /= 0 .and. cidp5 /= 0 ) then
      rdpz4 = 0.0
      rdpz5 = 0.0
      dpz4 = 0.0
      dpz5 = 0.0
      m4 = 1.0d0
      m5 = 1.0d0
    endif

    rdpx = 0.5d0*(rdpx0 + rdpx1)
    rdpy = 0.5d0*(rdpy2 + rdpy3)
    rdpz = 0.5d0*(rdpz4 + rdpz5)

    ux(i, j, k) = ux(i, j, k) + rdpx*dt
    uy(i, j, k) = uy(i, j, k) + rdpy*dt
    uz(i, j, k) = uz(i, j, k) + rdpz*dt

    if( pidp /= 1 ) then
      ux(i, j, k) = 0.0
      uy(i, j, k) = 0.0
      uz(i, j, k) = 0.0
    endif

  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_remove_p

subroutine bcut_calc_ab_p( &
                Ap, Aw, Ae, As, An, Ab, At, b, &
                vw, ve, vs, vn, vb, vt, &
                p0_, &
                ux, uy, uz, &
                c0, c1, c2, c3, c4, c5, &
                cid0, cid1, cid2, cid3, cid4, cid5, &
                pid, &
                rhof, &
                dx, dt, &
                sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: Ap, Aw, Ae, As, An, Ab, At, b
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: vw, ve, vs, vn, vb, vt
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: p0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: ux, uy, uz
  real                                                    :: vx0, vx1, vy2, vy3, vz4, vz5
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
  integer                  :: cidp
  integer                  :: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
  integer                  :: pidp, pidw, pide, pids, pidn, pidb, pidt
  real                    :: rhof
  real                    :: dx, dt
  real                    :: d0, d1, d2, d3, d4, d5
  real                    :: r0, r1, r2, r3, r4, r5
  real                    :: m0, m1, m2, m3, m4, m5
  real                    :: l0, l1, l2, l3, l4, l5
  real                    :: divv
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp           private(vx0, vx1, vy2, vy3, vz4, vz5) &
!$omp           private(cidp) &
!$omp           private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp           private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp           private(d0, d1, d2, d3, d4, d5) &
!$omp           private(r0, r1, r2, r3, r4, r5) &
!$omp           private(m0, m1, m2, m3, m4, m5) &
!$omp           private(l0, l1, l2, l3, l4, l5) &
!$omp           private(divv)
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
    r0 = 1.0d0/rhof
    r1 = 1.0d0/rhof
    r2 = 1.0d0/rhof
    r3 = 1.0d0/rhof
    r4 = 1.0d0/rhof
    r5 = 1.0d0/rhof

    vx0 = 0.5d0*(ux(i, j, k) + ux(i-1, j, k))
    vx1 = 0.5d0*(ux(i, j, k) + ux(i+1, j, k))
    vy2 = 0.5d0*(uy(i, j, k) + uy(i, j-1, k))
    vy3 = 0.5d0*(uy(i, j, k) + uy(i, j+1, k))
    vz4 = 0.5d0*(uz(i, j, k) + uz(i, j, k-1))
    vz5 = 0.5d0*(uz(i, j, k) + uz(i, j, k+1))

    m0 = 0.0d0
    m1 = 0.0d0
    m2 = 0.0d0
    m3 = 0.0d0
    m4 = 0.0d0
    m5 = 0.0d0

    d0 = c0(i, j, k)
    d1 = c1(i, j, k)
    d2 = c2(i, j, k)
    d3 = c3(i, j, k)
    d4 = c4(i, j, k)
    d5 = c5(i, j, k)

    cidp0 = cid0(i, j, k)
    cidp1 = cid1(i, j, k)
    cidp2 = cid2(i, j, k)
    cidp3 = cid3(i, j, k)
    cidp4 = cid4(i, j, k)
    cidp5 = cid5(i, j, k)

    pidp = pid(i, j, k)

    if( cidp0 /= 0 ) then
      r0 = 0.0d0
      vx0 = 0.0d0
      vx0 = vx1*(d0 - 0.5d0)/(d0 + 0.5d0)
      r1 = (1.0d0/rhof)/(d0 + 0.5d0)
      m0 = 1.0d0
    endif

    if( cidp1 /= 0 ) then
      r1 = 0.0d0
      vx1 = 0.0d0
      vx1 = vx0*(d1 - 0.5d0)/(d1 + 0.5d0)
      r0 = (1.0d0/rhof)/(d1 + 0.5d0)
      m1 = 1.0d0
    endif

    if( cidp2 /= 0 ) then
      r2 = 0.0d0
      vy2 = 0.0d0
      vy2 = vy3*(d2 - 0.5d0)/(d2 + 0.5d0)
      r3 = (1.0d0/rhof)/(d2 + 0.5d0)
      m2 = 1.0d0
    endif

    if( cidp3 /= 0 ) then
      r3 = 0.0d0
      vy3 = 0.0d0
      vy3 = vy2*(d3 - 0.5d0)/(d3 + 0.5d0)
      r2 = (1.0d0/rhof)/(d3 + 0.5d0)
      m3 = 1.0d0
    endif

    if( cidp4 /= 0 ) then
      r4 = 0.0d0
      vz4 = 0.0d0
      vz4 = vz5*(d4 - 0.5d0)/(d4 + 0.5d0)
      r5 = (1.0d0/rhof)/(d4 + 0.5d0)
      m4 = 1.0d0
    endif

    if( cidp5 /= 0 ) then
      r5 = 0.0d0
      vz5 = 0.0d0
      vz5 = vz4*(d5 - 0.5d0)/(d5 + 0.5d0)
      r4 = (1.0d0/rhof)/(d5 + 0.5d0)
      m5 = 1.0d0
    endif

    if( cidp0 /= 0 .and. cidp1 /= 0 ) then
      r0 = 0.0d0
      r1 = 0.0d0
      vx0 = 0.5d0*(2.0d0 - 1.0d0/d0)*ux(i, j, k)
      vx1 = 0.5d0*(2.0d0 - 1.0d0/d1)*ux(i, j, k)
      vx0 = 0.0d0
      vx1 = 0.0d0
      m0 = 1.0d0
      m1 = 1.0d0
    endif

    if( cidp2 /= 0 .and. cidp3 /= 0 ) then
      r2 = 0.0d0
      r3 = 0.0d0
      vy2 = 0.5d0*(2.0d0 - 1.0d0/d2)*uy(i, j, k)
      vy3 = 0.5d0*(2.0d0 - 1.0d0/d3)*uy(i, j, k)
      vy2 = 0.0d0
      vy3 = 0.0d0
      m2 = 1.0d0
      m3 = 1.0d0
    endif

    if( cidp4 /= 0 .and. cidp5 /= 0 ) then
      r4 = 0.0d0
      r5 = 0.0d0
      vz4 = 0.5d0*(2.0d0 - 1.0d0/d4)*uz(i, j, k)
      vz5 = 0.5d0*(2.0d0 - 1.0d0/d5)*uz(i, j, k)
      vz4 = 0.0d0
      vz5 = 0.0d0
      m4 = 1.0d0
      m5 = 1.0d0
    endif

    l0 = r0/(dx*dx)
    l1 = r1/(dx*dx)
    l2 = r2/(dx*dx)
    l3 = r3/(dx*dx)
    l4 = r4/(dx*dx)
    l5 = r5/(dx*dx)

    divv = (vx1 - vx0 + vy3 - vy2 + vz5 - vz4)/dx

    Ap(i, j, k) = - (l0 + l1) - (l2 + l3) - (l4 + l5)
    Aw(i, j, k) = l0
    Ae(i, j, k) = l1
    As(i, j, k) = l2
    An(i, j, k) = l3
    Ab(i, j, k) = l4
    At(i, j, k) = l5

     b(i, j, k) = divv/dt

    vw(i, j, k) = vx0
    ve(i, j, k) = vx1
    vs(i, j, k) = vy2
    vn(i, j, k) = vy3
    vb(i, j, k) = vz4
    vt(i, j, k) = vz5

    if( pidp /= 1 ) then
      Ap(i, j, k) =-6.0d0
      Aw(i, j, k) = 1.0d0
      Ae(i, j, k) = 1.0d0
      As(i, j, k) = 1.0d0
      An(i, j, k) = 1.0d0
      Ab(i, j, k) = 1.0d0
      At(i, j, k) = 1.0d0
      b (i, j, k) = 0.0d0
      vw(i, j, k) = 0.0d0
      ve(i, j, k) = 0.0d0
      vs(i, j, k) = 0.0d0
      vn(i, j, k) = 0.0d0
      vb(i, j, k) = 0.0d0
      vt(i, j, k) = 0.0d0
    endif

  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_calc_ab_p

subroutine bcut_corr_u( &
                ux0_, uy0_, uz0_, &
                vw, ve, vs, vn, vb, vt, &
                lapp, &
                p0_, &
                c0, c1, c2, c3, c4, c5, &
                cid0, cid1, cid2, cid3, cid4, cid5, &
                pid, &
                rhof, &
                dx, dt, &
                sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: ux0_, uy0_, uz0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: vw, ve, vs, vn, vb, vt
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: lapp
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: p0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
  integer                  :: cidp
  integer                  :: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
  integer                  :: pidp, pidw, pide, pids, pidn, pidb, pidt
  real                    :: rhof
  real                    :: dx, dt
  real                    :: d0, d1, d2, d3, d4, d5
  real                    :: r0, r1, r2, r3, r4, r5
  real                    :: m0, m1, m2, m3, m4, m5
  real                    :: l0, l1, l2, l3, l4, l5
  real                    :: dpx0, dpx1, dpy2, dpy3, dpz4, dpz5
  real                    :: rdpx0, rdpx1, rdpy2, rdpy3, rdpz4, rdpz5
  real                    :: rdpx, rdpy, rdpz
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp           private(cidp) &
!$omp           private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp           private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp           private(d0, d1, d2, d3, d4, d5) &
!$omp           private(r0, r1, r2, r3, r4, r5) &
!$omp           private(m0, m1, m2, m3, m4, m5) &
!$omp           private(l0, l1, l2, l3, l4, l5) &
!$omp           private(dpx0, dpx1, dpy2, dpy3, dpz4, dpz5) &
!$omp           private(rdpx0, rdpx1, rdpy2, rdpy3, rdpz4, rdpz5) &
!$omp           private(rdpx, rdpy, rdpz)
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
    r0 = 1.0d0/rhof
    r1 = 1.0d0/rhof
    r2 = 1.0d0/rhof
    r3 = 1.0d0/rhof
    r4 = 1.0d0/rhof
    r5 = 1.0d0/rhof

    dpx0 = (p0_(i, j, k) - p0_(i-1, j, k))/dx
    dpx1 = (p0_(i+1, j, k) - p0_(i, j, k))/dx
    dpy2 = (p0_(i, j, k) - p0_(i, j-1, k))/dx
    dpy3 = (p0_(i, j+1, k) - p0_(i, j, k))/dx
    dpz4 = (p0_(i, j, k) - p0_(i, j, k-1))/dx
    dpz5 = (p0_(i, j, k+1) - p0_(i, j, k))/dx

    rdpx0 = r0*dpx0
    rdpx1 = r1*dpx1
    rdpy2 = r2*dpy2
    rdpy3 = r3*dpy3
    rdpz4 = r4*dpz4
    rdpz5 = r5*dpz5

    d0 = c0(i, j, k)
    d1 = c1(i, j, k)
    d2 = c2(i, j, k)
    d3 = c3(i, j, k)
    d4 = c4(i, j, k)
    d5 = c5(i, j, k)

    cidp0 = cid0(i, j, k)
    cidp1 = cid1(i, j, k)
    cidp2 = cid2(i, j, k)
    cidp3 = cid3(i, j, k)
    cidp4 = cid4(i, j, k)
    cidp5 = cid5(i, j, k)

    pidp = pid(i, j, k)

    m0 = 0.0d0
    m1 = 0.0d0
    m2 = 0.0d0
    m3 = 0.0d0
    m4 = 0.0d0
    m5 = 0.0d0
    if( cidp0 /= 0 ) then
      rdpx0 = rdpx1
      rdpx0 = 0.0d0
      rdpx0 = rdpx1*(d0 - 0.5d0)/(d0 + 0.5d0)
      dpx0 = 0.0d0
      dpx0 = dpx1*(d0 - 0.5d0)/(d0 + 0.5d0)
      m0 = 1.0d0
    endif

    if( cidp1 /= 0 ) then
      rdpx1 = rdpx0
      rdpx1 = 0.0d0
      rdpx1 = rdpx0*(d1 - 0.5d0)/(d1 + 0.5d0)
      dpx1 = 0.0d0
      dpx1 = dpx0*(d1 - 0.5d0)/(d1 + 0.5d0)
      m1 = 1.0d0
    endif

    if( cidp2 /= 0 ) then
      rdpy2 = rdpy3
      rdpy2 = 0.0d0
      rdpy2 = rdpy3*(d2 - 0.5d0)/(d2 + 0.5d0)
      dpy2 = 0.0d0
      dpy2 = dpy3*(d2 - 0.5d0)/(d2 + 0.5d0)
      m2 = 1.0d0
    endif

    if( cidp3 /= 0 ) then
      rdpy3 = rdpy2
      rdpy3 = 0.0d0
      rdpy3 = rdpy2*(d3 - 0.5d0)/(d3 + 0.5d0)
      dpy3 = 0.0d0
      dpy3 = dpy2*(d3 - 0.5d0)/(d3 + 0.5d0)
      m3 = 1.0d0
    endif

    if( cidp4 /= 0 ) then
      rdpz4 = rdpz5
      rdpz4 = 0.0d0
      rdpz4 = rdpz5*(d4 - 0.5d0)/(d4 + 0.5d0)
      dpz4 = 0.0d0
      dpz4 = dpz5*(d4 - 0.5d0)/(d4 + 0.5d0)
      m4 = 1.0d0
    endif

    if( cidp5 /= 0 ) then
      rdpz5 = rdpz4
      rdpz5 = 0.0d0
      rdpz5 = rdpz4*(d5 - 0.5d0)/(d5 + 0.5d0)
      dpz5 = 0.0d0
      dpz5 = dpz4*(d5 - 0.5d0)/(d5 + 0.5d0)
      m5 = 1.0d0
    endif

    if( cidp0 /= 0 .and. cidp1 /= 0 ) then
      rdpx0 = 0.0
      rdpx1 = 0.0
      dpx0 = 0.0
      dpx1 = 0.0
      m0 = 1.0d0
      m1 = 1.0d0
    endif

    if( cidp2 /= 0 .and. cidp3 /= 0 ) then
      rdpy2 = 0.0
      rdpy3 = 0.0
      dpy2 = 0.0
      dpy3 = 0.0
      m2 = 1.0d0
      m3 = 1.0d0
    endif

    if( cidp4 /= 0 .and. cidp5 /= 0 ) then
      rdpz4 = 0.0
      rdpz5 = 0.0
      dpz4 = 0.0
      dpz5 = 0.0
      m4 = 1.0d0
      m5 = 1.0d0
    endif

    rdpx = 0.5d0*(rdpx0 + rdpx1)
    rdpy = 0.5d0*(rdpy2 + rdpy3)
    rdpz = 0.5d0*(rdpz4 + rdpz5)

    ux0_(i, j, k) = ux0_(i, j, k) - rdpx*dt
    uy0_(i, j, k) = uy0_(i, j, k) - rdpy*dt
    uz0_(i, j, k) = uz0_(i, j, k) - rdpz*dt

    vw(i, j, k) = vw(i, j, k) - rdpx0*dt
    ve(i, j, k) = ve(i, j, k) - rdpx1*dt
    vs(i, j, k) = vs(i, j, k) - rdpy2*dt
    vn(i, j, k) = vn(i, j, k) - rdpy3*dt
    vb(i, j, k) = vb(i, j, k) - rdpz4*dt
    vt(i, j, k) = vt(i, j, k) - rdpz5*dt

    lapp(i, j, k) = (rdpx1 - rdpx0 + rdpy3 - rdpy2 + rdpz5 - rdpz4)/dx*rhof
    lapp(i, j, k) = (dpx1 - dpx0 + dpy3 - dpy2 + dpz5 - dpz4)/dx
    lapp(i, j, k) = ( p0_(i+1, j, k) + p0_(i-1, j, k) &
                    + p0_(i, j+1, k) + p0_(i, j-1, k) &
                    + p0_(i, j, k+1) + p0_(i, j, k-1) &
                    - 6.0*p0_(i, j, k) )/(dx*dx)

    if( pidp /= 1 ) then
      ux0_(i, j, k) = 0.0d0
      uy0_(i, j, k) = 0.0d0
      uz0_(i, j, k) = 0.0d0

      vw(i, j, k) = 0.0d0
      ve(i, j, k) = 0.0d0
      vs(i, j, k) = 0.0d0
      vn(i, j, k) = 0.0d0
      vb(i, j, k) = 0.0d0
      vt(i, j, k) = 0.0d0
      lapp(i, j, k) = 0.0d0
    endif

    if( cidp0 /= 0 .or. &
        cidp1 /= 0 .or. &
        cidp2 /= 0 .or. &
        cidp3 /= 0 .or. &
        cidp4 /= 0 .or. &
        cidp5 /= 0 ) then
      lapp(i, j, k) = 0.0d0
    endif

  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_corr_u

subroutine bcut_calc_ab_t( &
                Ap, Aw, Ae, As, An, Ab, At, b, &
                t0_, &
                tc0_, tcp_, &
                td0_, &
                c0, c1, c2, c3, c4, c5, &
                cid0, cid1, cid2, cid3, cid4, cid5, &
                pid, &
                rhof, rhos, &
                cpf, cps, &
                kf, ks, &
                dx, dt, &
                Tc, &
                sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: Ap, Aw, Ae, As, An, Ab, At, b
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: t0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: tc0_, tcp_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: td0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
  integer                  :: cidp
  integer                  :: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
  integer                  :: pidp, pidw, pide, pids, pidn, pidb, pidt
  real                    :: rhof, rhos
  real                    :: cpf, cps
  real                    :: kf, ks
  real                    :: dx, dt
  real                    :: Tc
  real                    :: d0, d1, d2, d3, d4, d5
  real                    :: k0, k1, k2, k3, k4, k5
  real                    :: m0, m1, m2, m3, m4, m5
  real                    :: l0, l1, l2, l3, l4, l5
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp           private(cidp) &
!$omp           private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp           private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp           private(d0, d1, d2, d3, d4, d5) &
!$omp           private(k0, k1, k2, k3, k4, k5) &
!$omp           private(m0, m1, m2, m3, m4, m5) &
!$omp           private(l0, l1, l2, l3, l4, l5) 
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
    k0 = kf
    k1 = kf
    k2 = kf
    k3 = kf
    k4 = kf
    k5 = kf

    m0 = 0.0
    m1 = 0.0
    m2 = 0.0
    m3 = 0.0
    m4 = 0.0
    m5 = 0.0

    d0 = c0(i, j, k)
    d1 = c1(i, j, k)
    d2 = c2(i, j, k)
    d3 = c3(i, j, k)
    d4 = c4(i, j, k)
    d5 = c5(i, j, k)

    cidp0 = cid0(i, j, k)
    cidp1 = cid1(i, j, k)
    cidp2 = cid2(i, j, k)
    cidp3 = cid3(i, j, k)
    cidp4 = cid4(i, j, k)
    cidp5 = cid5(i, j, k)

    pidp = pid(i, j, k)

    if( cidp0 /= 0 ) then
      k0 = kf/d0*2.0/(d0 + d1)
      k1 = kf/d1*2.0/(d0 + d1)
      m0 = 1.0
    endif

    if( cidp1 /= 0 ) then
      k0 = kf/d0*2.0/(d0 + d1)
      k1 = kf/d1*2.0/(d0 + d1)
      m1 = 1.0
    endif

    if( cidp2 /= 0 ) then
      k2 = kf/d2*2.0/(d2 + d3)
      k3 = kf/d3*2.0/(d2 + d3)
      m2 = 1.0
    endif

    if( cidp3 /= 0 ) then
      k2 = kf/d2*2.0/(d2 + d3)
      k3 = kf/d3*2.0/(d2 + d3)
      m3 = 1.0
    endif

    if( cidp4 /= 0 ) then
      k4 = kf/d4*2.0/(d4 + d5)
      k5 = kf/d5*2.0/(d4 + d5)
      m4 = 1.0
    endif

    if( cidp5 /= 0 ) then
      k4 = kf/d4*2.0/(d4 + d5)
      k5 = kf/d5*2.0/(d4 + d5)
      m5 = 1.0
    endif

    l0 = k0/(rhof*cpf)*dt/(dx*dx)
    l1 = k1/(rhof*cpf)*dt/(dx*dx)
    l2 = k2/(rhof*cpf)*dt/(dx*dx)
    l3 = k3/(rhof*cpf)*dt/(dx*dx)
    l4 = k4/(rhof*cpf)*dt/(dx*dx)
    l5 = k5/(rhof*cpf)*dt/(dx*dx)

    Ap(i, j, k) = 1.0d0 + 0.5d0*(l0 + l1) &
                        + 0.5d0*(l2 + l3) &
                        + 0.5d0*(l4 + l5)
    Aw(i, j, k) = -0.5d0*l0*(1.0d0 - m0)
    Ae(i, j, k) = -0.5d0*l1*(1.0d0 - m1)
    As(i, j, k) = -0.5d0*l2*(1.0d0 - m2)
    An(i, j, k) = -0.5d0*l3*(1.0d0 - m3)
    Ab(i, j, k) = -0.5d0*l4*(1.0d0 - m4)
    At(i, j, k) = -0.5d0*l5*(1.0d0 - m5)
     b(i, j, k) = t0_(i, j, k) &
                  + 0.5d0*l0*Tc*m0 &
                  + 0.5d0*l1*Tc*m1 &
                  + 0.5d0*l2*Tc*m2 &
                  + 0.5d0*l3*Tc*m3 &
                  + 0.5d0*l4*Tc*m4 &
                  + 0.5d0*l5*Tc*m5 &
                  - (1.5d0*tc0_(i, j, k) - 0.5d0*tcp_(i, j, k))*dt &
                  + 0.5d0*td0_(i, j, k)*dt

    if( pidp /= 1 ) then
      Ap(i, j, k) = 1.0
      Aw(i, j, k) = 0.0
      Ae(i, j, k) = 0.0
      As(i, j, k) = 0.0
      An(i, j, k) = 0.0
      Ab(i, j, k) = 0.0
      At(i, j, k) = 0.0
      b (i, j, k) = Tc
    endif

  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_calc_ab_t

subroutine bcut_calc_abd_t( &
                Ap, Aw, Ae, As, An, Ab, At, b, &
                t0_, &
                tc0_, tcp_, &
                td0_, &
                c0, c1, c2, c3, c4, c5, &
                cid0, cid1, cid2, cid3, cid4, cid5, &
                pid, &
								nnum, &
								nx_, ny_, nz_, &
								nidx0, nidx1, nidx2, nidx3, nidx4, nidx5, &
                rhof, rhos, &
                cpf, cps, &
                kf, ks, &
                bc_n, &
                bc_type, &
                bc_value, &
                org, &
                dx, dt, &
                sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: Ap, Aw, Ae, As, An, Ab, At, b
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: t0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: tc0_, tcp_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: td0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
	integer									 :: nnum
	real, dimension(0:nnum-1) :: nx_, ny_, nz_
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: nidx0, nidx1, nidx2, nidx3, nidx4, nidx5
  integer                  :: cidp
  integer                  :: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
  integer                  :: pidp, pidw, pide, pids, pidn, pidb, pidt
  real                    :: rhof, rhos
  real                    :: cpf, cps
  real                    :: kf, ks
  integer                      :: bc_n
  integer, dimension(0:bc_n-1):: bc_type
  real, dimension(0:bc_n-1)   :: bc_value
	real										:: eps
  real                    :: dx, dt
  real, dimension(3)      :: org
  real                    :: d0, d1, d2, d3, d4, d5
  real                    :: k0, k1, k2, k3, k4, k5
	real										:: k_p, k_w, k_e, k_s, k_n, k_b, k_t
	real										:: rho_p
	real										:: cp_p
  real                    :: m0, m1, m2, m3, m4, m5
  real                    :: l0, l1, l2, l3, l4, l5
  real                    :: bp, b0, b1, b2, b3, b4, b5
  real                    :: tp, t0, t1, t2, t3, t4, t5
  real                    :: nx, ny, nz
  real                    :: x, y, z, r2, rl, xi, yi, zi, theta, phi
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
	eps = 0.01
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp           private(cidp) &
!$omp           private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp           private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp           private(d0, d1, d2, d3, d4, d5) &
!$omp           private(k0, k1, k2, k3, k4, k5) &
!$omp						private(k_p, k_w, k_e, k_s, k_n, k_b, k_t) &
!$omp						private(rho_p, cp_p) &
!$omp           private(m0, m1, m2, m3, m4, m5) &
!$omp           private(l0, l1, l2, l3, l4, l5) &
!$omp           private(bp, b0, b1, b2, b3, b4, b5) &
!$omp           private(tp, t0, t1, t2, t3, t4, t5) &
!$omp           private(nx, ny, nz) &
!$omp           private(x, y, z, r2, rl, xi, yi, zi, theta, phi)
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
    d0 = c0(i, j, k)
    d1 = c1(i, j, k)
    d2 = c2(i, j, k)
    d3 = c3(i, j, k)
    d4 = c4(i, j, k)
    d5 = c5(i, j, k)

    cidp0 = cid0(i, j, k)
    cidp1 = cid1(i, j, k)
    cidp2 = cid2(i, j, k)
    cidp3 = cid3(i, j, k)
    cidp4 = cid4(i, j, k)
    cidp5 = cid5(i, j, k)

    pidp = pid(i, j, k)

    m0 = 0.0
    m1 = 0.0
    m2 = 0.0
    m3 = 0.0
    m4 = 0.0
    m5 = 0.0
    if( cidp0 /= 0 ) then
      m0 = 1.0
    endif
    if( cidp1 /= 0 ) then
      m1 = 1.0
    endif
    if( cidp2 /= 0 ) then
      m2 = 1.0
    endif
    if( cidp3 /= 0 ) then
      m3 = 1.0
    endif
    if( cidp4 /= 0 ) then
      m4 = 1.0
    endif
    if( cidp5 /= 0 ) then
      m5 = 1.0
    endif

		rho_p = rhof
		cp_p  = cpf
		k_p   = kf
		k_w   = kf
		k_e   = kf
		k_s   = kf
		k_n   = kf
		k_b   = kf
		k_t   = kf
		if( pid(i, j, k) /= 1 ) then
			rho_p = rhos
			cp_p  = cps
			k_p   = ks
		end if
		if( pid(i-1, j, k) /= 1 ) then
			k_w   = ks
		end if
		if( pid(i+1, j, k) /= 1 ) then
			k_e   = ks
		end if
		if( pid(i, j-1, k) /= 1 ) then
			k_s   = ks
		end if
		if( pid(i, j+1, k) /= 1 ) then
			k_n   = ks
		end if
		if( pid(i, j, k-1) /= 1 ) then
			k_b   = ks
		end if
		if( pid(i, j, k+1) /= 1 ) then
			k_t   = ks
		end if

    k0 = k_p
    k1 = k_p
    k2 = k_p
    k3 = k_p
    k4 = k_p
    k5 = k_p
    tp = t0_(i, j, k)
    t0 = t0_(i-1, j, k) 
    t1 = t0_(i+1, j, k) 
    t2 = t0_(i, j-1, k) 
    t3 = t0_(i, j+1, k) 
    t4 = t0_(i, j, k-1) 
    t5 = t0_(i, j, k+1) 
    b0 = 0.0
    b1 = 0.0
    b2 = 0.0
    b3 = 0.0
    b4 = 0.0
    b5 = 0.0
    if( bc_type(cidp0) == 0 ) then
      k0 = k_p/d0*2.0/(d0 + d1)
      k1 = k_p/d1*2.0/(d0 + d1)
			if( d0 < eps ) then
				k0 = 0.0
				k1 = k_p
			end if
      t0 = bc_value(cidp0)
			b0 = 0.0
    else if( bc_type(cidp0) == 1 ) then
			nx = nx_( nidx0(i, j, k) )

      k0 = 0.0
      k1 = k_p/(d0 + 0.5)
      t0 = 0.0
      b0 = - nx*bc_value(cidp0)*dx/(d0 + 0.5)
    else if( bc_type(cidp0) == 2 ) then
      k0 = k_p*k_w/(k_p*(1.0 - d0) + k_w*d0)
    end if

    if( bc_type(cidp1) == 0 ) then
      k0 = k_p/d0*2.0/(d0 + d1)
      k1 = k_p/d1*2.0/(d0 + d1)
			if( d1 < eps ) then
				k0 = k_p
				k1 = 0.0
			end if
      t1 = bc_value(cidp1)
			b1 = 0.0
    else if( bc_type(cidp1) == 1 ) then
			nx = nx_( nidx1(i, j, k) )

      k0 = k_p/(d1 + 0.5)
      k1 = 0.0
      t1 = 0.0
      b1 = + nx*bc_value(cidp1)*dx/(d1 + 0.5)
    else if( bc_type(cidp1) == 2 ) then
      k1 = k_p*k_e/(k_p*(1.0 - d1) + k_e*d1)
    end if

    if( bc_type(cidp2) == 0 ) then
      k2 = k_p/d2*2.0/(d2 + d3)
      k3 = k_p/d3*2.0/(d2 + d3)
			if( d2 < eps ) then
				k2 = 0.0
				k3 = k_p
			end if
      t2 = bc_value(cidp2)
			b2 = 0.0
    else if( bc_type(cidp2) == 1 ) then
			ny = ny_( nidx2(i, j, k) )

      k2 = 0.0
      k3 = k_p/(d2 + 0.5)
      t2 = 0.0
      b2 = - ny*bc_value(cidp2)*dx/(d2 + 0.5)
    else if( bc_type(cidp2) == 2 ) then
      k2 = k_p*k_s/(k_p*(1.0 - d2) + k_s*d2)
    end if

    if( bc_type(cidp3) == 0 ) then
      k2 = k_p/d2*2.0/(d2 + d3)
      k3 = k_p/d3*2.0/(d2 + d3)
			if( d3 < eps ) then
				k2 = k_p
				k3 = 0.0
			end if
      t3 = bc_value(cidp3)
			b3 = 0.0
    else if( bc_type(cidp3) == 1 ) then
			ny = ny_( nidx3(i, j, k) )

      k2 = k_p/(d3 + 0.5)
      k3 = 0.0
      t3 = 0.0
      b3 = + ny*bc_value(cidp3)*dx/(d3 + 0.5)
    else if( bc_type(cidp3) == 2 ) then
      k3 = k_p*k_n/(k_p*(1.0 - d3) + k_n*d3)
    end if

    if( bc_type(cidp4) == 0 ) then
      k4 = k_p/d4*2.0/(d4 + d5)
      k5 = k_p/d5*2.0/(d4 + d5)
			if( d4 < eps ) then
				k4 = 0.0
				k5 = k_p
			end if
      t4 = bc_value(cidp4)
			b4 = 0.0
    else if( bc_type(cidp4) == 1 ) then
			nz = nz_( nidx4(i, j, k) )

      k4 = 0.0
      k5 = k_p/(d4 + 0.5)
      t4 = 0.0
      b4 = - nz*bc_value(cidp4)*dx/(d4 + 0.5)
    else if( bc_type(cidp4) == 2 ) then
      k4 = k_p*k_b/(k_p*(1.0 - d4) + k_b*d4)
    end if

    if( bc_type(cidp5) == 0 ) then
      k4 = k_p/d4*2.0/(d4 + d5)
      k5 = k_p/d5*2.0/(d4 + d5)
			if( d5 < eps ) then
				k4 = k_p
				k5 = 0.0
			end if
      t5 = bc_value(cidp5)
			b5 = 0.0
    else if( bc_type(cidp5) == 1 ) then
			nz = nz_( nidx5(i, j, k) )

      k4 = k_p/(d5 + 0.5)
      k5 = 0.0
      t5 = 0.0
      b5 = + nz*bc_value(cidp5)*dx/(d5 + 0.5)
    else if( bc_type(cidp5) == 2 ) then
      k5 = k_p*k_t/(k_p*(1.0 - d5) + k_t*d5)
    end if

    l0 = k0/(rho_p*cp_p)/(dx*dx)*dt
    l1 = k1/(rho_p*cp_p)/(dx*dx)*dt
    l2 = k2/(rho_p*cp_p)/(dx*dx)*dt
    l3 = k3/(rho_p*cp_p)/(dx*dx)*dt
    l4 = k4/(rho_p*cp_p)/(dx*dx)*dt
    l5 = k5/(rho_p*cp_p)/(dx*dx)*dt

    td0_(i, j, k) = ( k1*(t1 - tp) &
                    - k0*(tp - t0) &
                    + k3*(t3 - tp) &
                    - k2*(tp - t2) &
                    + k5*(t5 - tp) &
                    - k4*(tp - t4) &
                    )/(rho_p*cp_p)/(dx*dx)

    bp = (b0 + b1 + b2 + b3 + b4 + b5)*k_p/(rho_p*cp_p)/(dx*dx)

    if( pidp /= 1 ) then
			tc0_(i, j, k) = 0.0
			tcp_(i, j, k) = 0.0
		end if

    Ap(i, j, k) = 1.0d0 + 0.5d0*(l0 + l1) &
                        + 0.5d0*(l2 + l3) &
                        + 0.5d0*(l4 + l5)
    Aw(i, j, k) = -0.5d0*l0*(1.0d0 - m0)
    Ae(i, j, k) = -0.5d0*l1*(1.0d0 - m1)
    As(i, j, k) = -0.5d0*l2*(1.0d0 - m2)
    An(i, j, k) = -0.5d0*l3*(1.0d0 - m3)
    Ab(i, j, k) = -0.5d0*l4*(1.0d0 - m4)
    At(i, j, k) = -0.5d0*l5*(1.0d0 - m5)
     b(i, j, k) = t0_(i, j, k) &
                  + 0.5d0*l0*t0*m0 &
                  + 0.5d0*l1*t1*m1 &
                  + 0.5d0*l2*t2*m2 &
                  + 0.5d0*l3*t3*m3 &
                  + 0.5d0*l4*t4*m4 &
                  + 0.5d0*l5*t5*m5 &
                  - (1.5d0*tc0_(i, j, k) - 0.5d0*tcp_(i, j, k))*dt &
                  + 0.5d0*td0_(i, j, k)*dt &
                  + bp*dt

    if( pidp /= 1 ) then
!      Ap(i, j, k) = 1.0
!      Aw(i, j, k) = 0.0
!      Ae(i, j, k) = 0.0
!      As(i, j, k) = 0.0
!      An(i, j, k) = 0.0
!      Ab(i, j, k) = 0.0
!      At(i, j, k) = 0.0
!      b (i, j, k) = 0.0
    endif

  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_calc_abd_t

subroutine bcut_update_t( &
                t0_, &
                tc0_, tcp_, &
                td0_, tdp_, &
                c0, c1, c2, c3, c4, c5, &
                cid0, cid1, cid2, cid3, cid4, cid5, &
                pid, &
                rhof, rhos, &
                cpf, cps, &
                kf, ks, &
                dx, dt, &
                sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: t0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: tc0_, tcp_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: td0_, tdp_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
  integer                  :: cidp
  integer                  :: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
  integer                  :: pidp, pidw, pide, pids, pidn, pidb, pidt
  real                    :: rhof, rhos
  real                    :: cpf, cps
  real                    :: kf, ks
  real                    :: dx, dt
  real                    :: Tc
  real                    :: d0, d1, d2, d3, d4, d5
  real                    :: k0, k1, k2, k3, k4, k5
  real                    :: m0, m1, m2, m3, m4, m5
  real                    :: l0, l1, l2, l3, l4, l5
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp           private(cidp) &
!$omp           private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp           private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp           private(d0, d1, d2, d3, d4, d5) &
!$omp           private(k0, k1, k2, k3, k4, k5) &
!$omp           private(m0, m1, m2, m3, m4, m5) &
!$omp           private(l0, l1, l2, l3, l4, l5) 
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
    pidp = pid(i, j, k)

    t0_(i, j, k) = t0_(i, j, k) &
                  - (1.5*tc0_(i, j, k) - 0.5*tcp_(i, j, k))*dt &
                  + (1.5*td0_(i, j, k) - 0.5*tdp_(i, j, k))*dt 

    if( pidp /= 1 ) then
      t0_(i, j, k) = 0.0
    endif

  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_update_t

subroutine bcut_calc_f_p_0( &
                fspx, &
                fspy, &
                fspz, &
                fsp, &
                p, &
                c0, c1, c2, c3, c4, c5, &
                cid0, cid1, cid2, cid3, cid4, cid5, &
                pid, &
                dx, dt, &
                sz, g)
  implicit none
  integer                  :: i, j, k
  integer                  :: ix, jx, kx
  integer                  :: g
  integer, dimension(3)    :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: fspx, fspy, fspz
  real, dimension(1:3)    :: fsp
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: p
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
  integer                  :: cidp
  integer                  :: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
  integer                  :: pidp, pidw, pide, pids, pidn, pidb, pidt
  real                    :: dx, dt
  real                    :: fsx, fsy, fsz
  real                    :: d0, d1, d2, d3, d4, d5
  real                    :: m0, m1, m2, m3, m4, m5
  real                    :: p0, p1, p2, p3, p4, p5
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  fsx = 0.0
  fsy = 0.0
  fsz = 0.0
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp           private(cidp) &
!$omp           private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp           private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp           private(d0, d1, d2, d3, d4, d5) &
!$omp           private(m0, m1, m2, m3, m4, m5) &
!$omp           private(p0, p1, p2, p3, p4, p5) 
!$omp do schedule(static, 1), &
!$omp     reduction(+:fsx, fsy, fsz)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
    p0 = 0.5*(p(i, j, k) + p(i-1, j, k))
    p1 = 0.5*(p(i, j, k) + p(i+1, j, k))
    p2 = 0.5*(p(i, j, k) + p(i, j-1, k))
    p3 = 0.5*(p(i, j, k) + p(i, j+1, k))
    p4 = 0.5*(p(i, j, k) + p(i, j, k-1))
    p5 = 0.5*(p(i, j, k) + p(i, j, k+1))

    m0 = 0.0d0
    m1 = 0.0d0
    m2 = 0.0d0
    m3 = 0.0d0
    m4 = 0.0d0
    m5 = 0.0d0

    d0 = c0(i, j, k)
    d1 = c1(i, j, k)
    d2 = c2(i, j, k)
    d3 = c3(i, j, k)
    d4 = c4(i, j, k)
    d5 = c5(i, j, k)

    cidp0 = cid0(i, j, k)
    cidp1 = cid1(i, j, k)
    cidp2 = cid2(i, j, k)
    cidp3 = cid3(i, j, k)
    cidp4 = cid4(i, j, k)
    cidp5 = cid5(i, j, k)

    pidp = pid(i, j, k)

    if( cidp0 /= 0 ) then
      p0 = p(i, j, k)
      p0 = p(i, j, k) + (p(i, j, k) - p(i+1, j, k))*d0*d0/(d0 + 0.5d0)
      p0 = p(i, j, k) + (p(i, j, k) - p(i+1, j, k))*d0*(d0 - 0.5d0)/(d0 + 0.5d0)
      p0 = p(i, j, k) + (p(i, j, k) - p(i+1, j, k))*0.5d0*(d0 - 0.5d0)/(d0 + 0.5d0)
      p0 = p(i, j, k) + (p(i, j, k) - p(i+1, j, k))*d0*d0/(d0 + 0.5d0)
      p0 = p(i, j, k) + (p(i, j, k) - p(i+1, j, k))*d0
      p0 = p(i, j, k) + (p(i, j, k) - p(i+1, j, k))*0.5d0
      m0 = 1.0d0
    endif

    if( cidp1 /= 0 ) then
      p1 = p(i, j, k)
      p1 = p(i, j, k) + (p(i, j, k) - p(i-1, j, k))*d1*(d1 - 0.5d0)/(d1 + 0.5d0) 
      p1 = p(i, j, k) + (p(i, j, k) - p(i-1, j, k))*0.5d0*(d1 - 0.5d0)/(d1 + 0.5d0) 
      p1 = p(i, j, k) + (p(i, j, k) - p(i-1, j, k))*d1*d1/(d1 + 0.5d0)
      p1 = p(i, j, k) + (p(i, j, k) - p(i-1, j, k))*d1 
      p1 = p(i, j, k) + (p(i, j, k) - p(i-1, j, k))*0.5d0
      m1 = 1.0d0
    endif

    if( cidp2 /= 0 ) then
      p2 = p(i, j, k)
      p2 = p(i, j, k) + (p(i, j, k) - p(i, j+1, k))*d2*(d2 - 0.5d0)/(d2 + 0.5d0)
      p2 = p(i, j, k) + (p(i, j, k) - p(i, j+1, k))*0.5d0*(d2 - 0.5d0)/(d2 + 0.5d0)
      p2 = p(i, j, k) + (p(i, j, k) - p(i, j+1, k))*d2*d2/(d2 + 0.5d0)
      p2 = p(i, j, k) + (p(i, j, k) - p(i, j+1, k))*d2
      p2 = p(i, j, k) + (p(i, j, k) - p(i, j+1, k))*0.5d0
      m2 = 1.0d0
    endif

    if( cidp3 /= 0 ) then
      p3 = p(i, j, k)
      p3 = p(i, j, k) + (p(i, j, k) - p(i, j-1, k))*d3*(d3 - 0.5d0)/(d3 + 0.5d0)
      p3 = p(i, j, k) + (p(i, j, k) - p(i, j-1, k))*0.5d0*(d3 - 0.5d0)/(d3 + 0.5d0)
      p3 = p(i, j, k) + (p(i, j, k) - p(i, j-1, k))*d3*d3/(d3 + 0.5d0)
      p3 = p(i, j, k) + (p(i, j, k) - p(i, j-1, k))*d3
      p3 = p(i, j, k) + (p(i, j, k) - p(i, j-1, k))*0.5d0
      m3 = 1.0d0
    endif

    if( cidp4 /= 0 ) then
      p4 = p(i, j, k)
      p4 = p(i, j, k) + (p(i, j, k) - p(i, j, k+1))*d4*(d4 - 0.5d0)/(d4 + 0.5d0)
      p4 = p(i, j, k) + (p(i, j, k) - p(i, j, k+1))*0.5d0*(d4 - 0.5d0)/(d4 + 0.5d0)
      p4 = p(i, j, k) + (p(i, j, k) - p(i, j, k+1))*d4*d4/(d4 + 0.5d0)
      p4 = p(i, j, k) + (p(i, j, k) - p(i, j, k+1))*d4
      p4 = p(i, j, k) + (p(i, j, k) - p(i, j, k+1))*0.5d0
      m4 = 1.0d0
    endif

    if( cidp5 /= 0 ) then
      p5 = p(i, j, k)
      p5 = p(i, j, k) + (p(i, j, k) - p(i, j, k-1))*d5*(d5 - 0.5d0)/(d5 + 0.5d0)
      p5 = p(i, j, k) + (p(i, j, k) - p(i, j, k-1))*0.5d0*(d5 - 0.5d0)/(d5 + 0.5d0)
      p5 = p(i, j, k) + (p(i, j, k) - p(i, j, k-1))*d5*d5/(d5 + 0.5d0)
      p5 = p(i, j, k) + (p(i, j, k) - p(i, j, k-1))*d5
      p5 = p(i, j, k) + (p(i, j, k) - p(i, j, k-1))*0.5d0
      m5 = 1.0d0
    endif

    if( cidp0 /= 0 .and. cidp1 /= 0 ) then
      p0 = p(i, j, k)
      p1 = p(i, j, k)
      m0 = 1.0d0
      m1 = 1.0d0
    endif

    if( cidp2 /= 0 .and. cidp3 /= 0 ) then
      p2 = p(i, j, k)
      p3 = p(i, j, k)
      m2 = 1.0d0
      m3 = 1.0d0
    endif

    if( cidp4 /= 0 .and. cidp5 /= 0 ) then
      p4 = p(i, j, k)
      p5 = p(i, j, k)
      m4 = 1.0d0
      m5 = 1.0d0
    endif


    if( pidp /= 1 ) then
      m0 = 0.0d0
      m1 = 0.0d0
      m2 = 0.0d0
      m3 = 0.0d0
      m4 = 0.0d0
      m5 = 0.0d0
    endif

    fspx(i, j, k) = (m1*p1 - m0*p0)
    fspy(i, j, k) = (m3*p3 - m2*p2)
    fspz(i, j, k) = (m5*p5 - m4*p4)

    fsx = fsx + fspx(i, j, k)*dx*dx
    fsy = fsy + fspy(i, j, k)*dx*dx
    fsz = fsz + fspz(i, j, k)*dx*dx
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
  fsp(1) = fsx
  fsp(2) = fsy
  fsp(3) = fsz
end subroutine bcut_calc_f_p_0

subroutine bcut_calc_f_v_0( &
                fsvx, &
                fsvy, &
                fsvz, &
                fsv, &
                ux, &
                uy, &
                uz, &
                c0, c1, c2, c3, c4, c5, &
                cid0, cid1, cid2, cid3, cid4, cid5, &
                pid, &
                rho, &
                mu, &
                dx, dt, &
                Uc, &
                sz, g)
  implicit none
  integer                  :: i, j, k
  integer                  :: ix, jx, kx
  integer                  :: g
  integer, dimension(3)    :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: fsvx, fsvy, fsvz
  real, dimension(1:3)    :: fsv
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: ux, uy, uz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
  real                    :: rho
  real                    :: mu
  real                    :: dx, dt
  real                    :: Uc
  integer                  :: cidp
  integer                  :: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
  integer                  :: pidp, pidw, pide, pids, pidn, pidb, pidt
  real                    :: d0, d1, d2, d3, d4, d5
  real                    :: m0, m1, m2, m3, m4, m5
  real                    :: uxp, uxw, uxe, uxs, uxn, uxb, uxt
  real                    :: uyp, uyw, uye, uys, uyn, uyb, uyt
  real                    :: uzp, uzw, uze, uzs, uzn, uzb, uzt
  real                    :: fsx, fsy, fsz
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  fsx = 0.0
  fsy = 0.0
  fsz = 0.0
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp           private(cidp) &
!$omp           private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp           private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp           private(d0, d1, d2, d3, d4, d5) &
!$omp           private(m0, m1, m2, m3, m4, m5) &
!$omp           private(uxp, uxw, uxe, uxs, uxn, uxb, uxt) &
!$omp           private(uyp, uyw, uye, uys, uyn, uyb, uyt) &
!$omp           private(uzp, uzw, uze, uzs, uzn, uzb, uzt) 
!$omp do schedule(dynamic, 1), &
!$omp     reduction(+:fsx, fsy, fsz)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
    uxp = ux(i, j, k)
    uxw = ux(i-1, j, k)
    uxe = ux(i+1, j, k)
    uxs = ux(i, j-1, k)
    uxn = ux(i, j+1, k)
    uxb = ux(i, j, k-1)
    uxt = ux(i, j, k+1)

    uyp = uy(i, j, k)
    uyw = uy(i-1, j, k)
    uye = uy(i+1, j, k)
    uys = uy(i, j-1, k)
    uyn = uy(i, j+1, k)
    uyb = uy(i, j, k-1)
    uyt = uy(i, j, k+1)

    uzp = uz(i, j, k)
    uzw = uz(i-1, j, k)
    uze = uz(i+1, j, k)
    uzs = uz(i, j-1, k)
    uzn = uz(i, j+1, k)
    uzb = uz(i, j, k-1)
    uzt = uz(i, j, k+1)

    m0 = 0.0d0
    m1 = 0.0d0
    m2 = 0.0d0
    m3 = 0.0d0
    m4 = 0.0d0
    m5 = 0.0d0

    d0 = c0(i, j, k)
    d1 = c1(i, j, k)
    d2 = c2(i, j, k)
    d3 = c3(i, j, k)
    d4 = c4(i, j, k)
    d5 = c5(i, j, k)

    cidp0 = cid0(i, j, k)
    cidp1 = cid1(i, j, k)
    cidp2 = cid2(i, j, k)
    cidp3 = cid3(i, j, k)
    cidp4 = cid4(i, j, k)
    cidp5 = cid5(i, j, k)

    pidp = pid(i, j, k)

    if( cidp0 /= 0 ) then
      uxw  = (1.0d0 - 1.0d0/d0)*uxp 
      uyw  = (1.0d0 - 1.0d0/d0)*uyp 
      uzw  = (1.0d0 - 1.0d0/d0)*uzp 
      m0 = 1.0d0
    endif
    if( cidp1 /= 0 ) then
      uxe  = (1.0d0 - 1.0d0/d1)*uxp
      uye  = (1.0d0 - 1.0d0/d1)*uyp
      uze  = (1.0d0 - 1.0d0/d1)*uzp
      m1 = 1.0d0
    endif
    if( cidp2 /= 0 ) then
      uxs  = (1.0d0 - 1.0d0/d2)*uxp 
      uys  = (1.0d0 - 1.0d0/d2)*uyp 
      uzs  = (1.0d0 - 1.0d0/d2)*uzp 
      m2 = 1.0d0
    endif
    if( cidp3 /= 0 ) then
      uxn  = (1.0d0 - 1.0d0/d3)*uxp
      uyn  = (1.0d0 - 1.0d0/d3)*uyp
      uzn  = (1.0d0 - 1.0d0/d3)*uzp
      m3 = 1.0d0
    endif
    if( cidp4 /= 0 ) then
      uxb  = (1.0d0 - 1.0d0/d4)*uxp
      uyb  = (1.0d0 - 1.0d0/d4)*uyp
      uzb  = (1.0d0 - 1.0d0/d4)*uzp
      m4 = 1.0d0
    endif
    if( cidp5 /= 0 ) then
      uxt  = (1.0d0 - 1.0d0/d5)*uxp
      uyt  = (1.0d0 - 1.0d0/d5)*uyp
      uzt  = (1.0d0 - 1.0d0/d5)*uzp
      m5 = 1.0d0
    endif

    if( pidp /= 1 ) then
      m0 = 0.0d0
      m1 = 0.0d0
      m2 = 0.0d0
      m3 = 0.0d0
      m4 = 0.0d0
      m5 = 0.0d0
    endif

    fsvx(i, j, k) = + mu*((uxp - uxw))/dx*m0*2.0 &
                    - mu*((uxe - uxp))/dx*m1*2.0 &
                    + mu*((uxp - uxs) + (uye - uyw)*0.5)/dx*m2 &
                    - mu*((uxn - uxp) + (uye - uyw)*0.5)/dx*m3 &
                    + mu*((uxp - uxb) + (uze - uzw)*0.5)/dx*m4 &
                    - mu*((uxt - uxp) + (uze - uzw)*0.5)/dx*m5
    fsvy(i, j, k) = + mu*((uyp - uyw) + (uxn - uxs)*0.5)/dx*m0 &
                    - mu*((uye - uyp) + (uxn - uxs)*0.5)/dx*m1 &
                    + mu*((uyp - uys))/dx*m2*2.0 &
                    - mu*((uyn - uyp))/dx*m3*2.0 &
                    + mu*((uyp - uyb) + (uzn - uzs)*0.5)/dx*m4 &
                    - mu*((uyt - uyp) + (uzn - uzs)*0.5)/dx*m5
    fsvz(i, j, k) = + mu*((uzp - uzw) + (uxt - uxb)*0.5)/dx*m0 &
                    - mu*((uze - uzp) + (uxt - uxb)*0.5)/dx*m1 &
                    + mu*((uzp - uzs) + (uyt - uyb)*0.5)/dx*m2 &
                    - mu*((uzn - uzp) + (uyt - uyb)*0.5)/dx*m3 &
                    + mu*((uzp - uzb))/dx*m4*2.0 &
                    - mu*((uzt - uzp))/dx*m5*2.0

    fsx = fsx + fsvx(i, j, k)*dx*dx
    fsy = fsy + fsvy(i, j, k)*dx*dx
    fsz = fsz + fsvz(i, j, k)*dx*dx
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
  fsv(1) = fsx
  fsv(2) = fsy
  fsv(3) = fsz
end subroutine bcut_calc_f_v_0

subroutine bcut_calc_f_p( &
                fspx, &
                fspy, &
                fspz, &
                fsp, &
								cid_target, &
                p, &
                c0, c1, c2, c3, c4, c5, &
                cid0, cid1, cid2, cid3, cid4, cid5, &
                pid, &
                dx, dt, &
                sz, g)
  implicit none
  integer                  :: i, j, k
  integer                  :: ix, jx, kx
  integer                  :: g
  integer, dimension(3)    :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: fspx, fspy, fspz
  real, dimension(1:3)    :: fsp
	integer									:: cid_target
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: p
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
  integer                  :: cidp
  integer                  :: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
  integer                  :: pidp, pidw, pide, pids, pidn, pidb, pidt
  real                    :: dx, dt
  real                    :: fsx, fsy, fsz
  real                    :: d0, d1, d2, d3, d4, d5
  real                    :: m0, m1, m2, m3, m4, m5
  real                    :: p0, p1, p2, p3, p4, p5
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  fsx = 0.0
  fsy = 0.0
  fsz = 0.0
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp           private(cidp) &
!$omp           private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp           private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp           private(d0, d1, d2, d3, d4, d5) &
!$omp           private(m0, m1, m2, m3, m4, m5) &
!$omp           private(p0, p1, p2, p3, p4, p5) 
!$omp do schedule(static, 1), &
!$omp     reduction(+:fsx, fsy, fsz)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
    p0 = 0.5*(p(i, j, k) + p(i-1, j, k))
    p1 = 0.5*(p(i, j, k) + p(i+1, j, k))
    p2 = 0.5*(p(i, j, k) + p(i, j-1, k))
    p3 = 0.5*(p(i, j, k) + p(i, j+1, k))
    p4 = 0.5*(p(i, j, k) + p(i, j, k-1))
    p5 = 0.5*(p(i, j, k) + p(i, j, k+1))

    m0 = 0.0d0
    m1 = 0.0d0
    m2 = 0.0d0
    m3 = 0.0d0
    m4 = 0.0d0
    m5 = 0.0d0

    d0 = c0(i, j, k)
    d1 = c1(i, j, k)
    d2 = c2(i, j, k)
    d3 = c3(i, j, k)
    d4 = c4(i, j, k)
    d5 = c5(i, j, k)

    cidp0 = cid0(i, j, k)
    cidp1 = cid1(i, j, k)
    cidp2 = cid2(i, j, k)
    cidp3 = cid3(i, j, k)
    cidp4 = cid4(i, j, k)
    cidp5 = cid5(i, j, k)

    pidp = pid(i, j, k)

    if( cidp0 == cid_target ) then
      p0 = p(i, j, k)
      p0 = p(i, j, k) + (p(i, j, k) - p(i+1, j, k))*d0*d0/(d0 + 0.5d0)
      p0 = p(i, j, k) + (p(i, j, k) - p(i+1, j, k))*d0*(d0 - 0.5d0)/(d0 + 0.5d0)
      p0 = p(i, j, k) + (p(i, j, k) - p(i+1, j, k))*0.5d0*(d0 - 0.5d0)/(d0 + 0.5d0)
      p0 = p(i, j, k) + (p(i, j, k) - p(i+1, j, k))*d0*d0/(d0 + 0.5d0)
      p0 = p(i, j, k) + (p(i, j, k) - p(i+1, j, k))*d0
      p0 = p(i, j, k) + (p(i, j, k) - p(i+1, j, k))*0.5d0
      m0 = 1.0d0
    endif

    if( cidp1 == cid_target ) then
      p1 = p(i, j, k)
      p1 = p(i, j, k) + (p(i, j, k) - p(i-1, j, k))*d1*(d1 - 0.5d0)/(d1 + 0.5d0) 
      p1 = p(i, j, k) + (p(i, j, k) - p(i-1, j, k))*0.5d0*(d1 - 0.5d0)/(d1 + 0.5d0) 
      p1 = p(i, j, k) + (p(i, j, k) - p(i-1, j, k))*d1*d1/(d1 + 0.5d0)
      p1 = p(i, j, k) + (p(i, j, k) - p(i-1, j, k))*d1 
      p1 = p(i, j, k) + (p(i, j, k) - p(i-1, j, k))*0.5d0
      m1 = 1.0d0
    endif

    if( cidp2 == cid_target ) then
      p2 = p(i, j, k)
      p2 = p(i, j, k) + (p(i, j, k) - p(i, j+1, k))*d2*(d2 - 0.5d0)/(d2 + 0.5d0)
      p2 = p(i, j, k) + (p(i, j, k) - p(i, j+1, k))*0.5d0*(d2 - 0.5d0)/(d2 + 0.5d0)
      p2 = p(i, j, k) + (p(i, j, k) - p(i, j+1, k))*d2*d2/(d2 + 0.5d0)
      p2 = p(i, j, k) + (p(i, j, k) - p(i, j+1, k))*d2
      p2 = p(i, j, k) + (p(i, j, k) - p(i, j+1, k))*0.5d0
      m2 = 1.0d0
    endif

    if( cidp3 == cid_target ) then
      p3 = p(i, j, k)
      p3 = p(i, j, k) + (p(i, j, k) - p(i, j-1, k))*d3*(d3 - 0.5d0)/(d3 + 0.5d0)
      p3 = p(i, j, k) + (p(i, j, k) - p(i, j-1, k))*0.5d0*(d3 - 0.5d0)/(d3 + 0.5d0)
      p3 = p(i, j, k) + (p(i, j, k) - p(i, j-1, k))*d3*d3/(d3 + 0.5d0)
      p3 = p(i, j, k) + (p(i, j, k) - p(i, j-1, k))*d3
      p3 = p(i, j, k) + (p(i, j, k) - p(i, j-1, k))*0.5d0
      m3 = 1.0d0
    endif

    if( cidp4 == cid_target ) then
      p4 = p(i, j, k)
      p4 = p(i, j, k) + (p(i, j, k) - p(i, j, k+1))*d4*(d4 - 0.5d0)/(d4 + 0.5d0)
      p4 = p(i, j, k) + (p(i, j, k) - p(i, j, k+1))*0.5d0*(d4 - 0.5d0)/(d4 + 0.5d0)
      p4 = p(i, j, k) + (p(i, j, k) - p(i, j, k+1))*d4*d4/(d4 + 0.5d0)
      p4 = p(i, j, k) + (p(i, j, k) - p(i, j, k+1))*d4
      p4 = p(i, j, k) + (p(i, j, k) - p(i, j, k+1))*0.5d0
      m4 = 1.0d0
    endif

    if( cidp5 == cid_target ) then
      p5 = p(i, j, k)
      p5 = p(i, j, k) + (p(i, j, k) - p(i, j, k-1))*d5*(d5 - 0.5d0)/(d5 + 0.5d0)
      p5 = p(i, j, k) + (p(i, j, k) - p(i, j, k-1))*0.5d0*(d5 - 0.5d0)/(d5 + 0.5d0)
      p5 = p(i, j, k) + (p(i, j, k) - p(i, j, k-1))*d5*d5/(d5 + 0.5d0)
      p5 = p(i, j, k) + (p(i, j, k) - p(i, j, k-1))*d5
      p5 = p(i, j, k) + (p(i, j, k) - p(i, j, k-1))*0.5d0
      m5 = 1.0d0
    endif

    if( pidp /= 1 ) then
      m0 = 0.0d0
      m1 = 0.0d0
      m2 = 0.0d0
      m3 = 0.0d0
      m4 = 0.0d0
      m5 = 0.0d0
    endif

    fspx(i, j, k) = (m1*p1 - m0*p0)
    fspy(i, j, k) = (m3*p3 - m2*p2)
    fspz(i, j, k) = (m5*p5 - m4*p4)

    fsx = fsx + fspx(i, j, k)*dx*dx
    fsy = fsy + fspy(i, j, k)*dx*dx
    fsz = fsz + fspz(i, j, k)*dx*dx
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
  fsp(1) = fsx
  fsp(2) = fsy
  fsp(3) = fsz
end subroutine bcut_calc_f_p

subroutine bcut_calc_f_v( &
                fsvx, &
                fsvy, &
                fsvz, &
                fsv, &
								cid_target, &
                ux, &
                uy, &
                uz, &
                c0, c1, c2, c3, c4, c5, &
                cid0, cid1, cid2, cid3, cid4, cid5, &
                pid, &
                rho, &
                mu, &
                dx, dt, &
                Uc, &
                sz, g)
  implicit none
  integer                  :: i, j, k
  integer                  :: ix, jx, kx
  integer                  :: g
  integer, dimension(3)    :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: fsvx, fsvy, fsvz
  real, dimension(1:3)    :: fsv
	integer									:: cid_target
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: ux, uy, uz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
  real                    :: rho
  real                    :: mu
  real                    :: dx, dt
  real                    :: Uc
  integer                  :: cidp
  integer                  :: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
  integer                  :: pidp, pidw, pide, pids, pidn, pidb, pidt
  real                    :: d0, d1, d2, d3, d4, d5
  real                    :: m0, m1, m2, m3, m4, m5
  real                    :: uxp, uxw, uxe, uxs, uxn, uxb, uxt
  real                    :: uyp, uyw, uye, uys, uyn, uyb, uyt
  real                    :: uzp, uzw, uze, uzs, uzn, uzb, uzt
  real                    :: fsx, fsy, fsz
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  fsx = 0.0
  fsy = 0.0
  fsz = 0.0
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp           private(cidp) &
!$omp           private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp           private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp           private(d0, d1, d2, d3, d4, d5) &
!$omp           private(m0, m1, m2, m3, m4, m5) &
!$omp           private(uxp, uxw, uxe, uxs, uxn, uxb, uxt) &
!$omp           private(uyp, uyw, uye, uys, uyn, uyb, uyt) &
!$omp           private(uzp, uzw, uze, uzs, uzn, uzb, uzt) 
!$omp do schedule(dynamic, 1), &
!$omp     reduction(+:fsx, fsy, fsz)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
    uxp = ux(i, j, k)
    uxw = ux(i-1, j, k)
    uxe = ux(i+1, j, k)
    uxs = ux(i, j-1, k)
    uxn = ux(i, j+1, k)
    uxb = ux(i, j, k-1)
    uxt = ux(i, j, k+1)

    uyp = uy(i, j, k)
    uyw = uy(i-1, j, k)
    uye = uy(i+1, j, k)
    uys = uy(i, j-1, k)
    uyn = uy(i, j+1, k)
    uyb = uy(i, j, k-1)
    uyt = uy(i, j, k+1)

    uzp = uz(i, j, k)
    uzw = uz(i-1, j, k)
    uze = uz(i+1, j, k)
    uzs = uz(i, j-1, k)
    uzn = uz(i, j+1, k)
    uzb = uz(i, j, k-1)
    uzt = uz(i, j, k+1)

    m0 = 0.0d0
    m1 = 0.0d0
    m2 = 0.0d0
    m3 = 0.0d0
    m4 = 0.0d0
    m5 = 0.0d0

    d0 = c0(i, j, k)
    d1 = c1(i, j, k)
    d2 = c2(i, j, k)
    d3 = c3(i, j, k)
    d4 = c4(i, j, k)
    d5 = c5(i, j, k)

    cidp0 = cid0(i, j, k)
    cidp1 = cid1(i, j, k)
    cidp2 = cid2(i, j, k)
    cidp3 = cid3(i, j, k)
    cidp4 = cid4(i, j, k)
    cidp5 = cid5(i, j, k)

    pidp = pid(i, j, k)

    if( cidp0 == cid_target ) then
      uxw  = (1.0d0 - 1.0d0/d0)*uxp 
      uyw  = (1.0d0 - 1.0d0/d0)*uyp 
      uzw  = (1.0d0 - 1.0d0/d0)*uzp 
      m0 = 1.0d0
    endif
    if( cidp1 == cid_target ) then
      uxe  = (1.0d0 - 1.0d0/d1)*uxp
      uye  = (1.0d0 - 1.0d0/d1)*uyp
      uze  = (1.0d0 - 1.0d0/d1)*uzp
      m1 = 1.0d0
    endif
    if( cidp2 == cid_target ) then
      uxs  = (1.0d0 - 1.0d0/d2)*uxp 
      uys  = (1.0d0 - 1.0d0/d2)*uyp 
      uzs  = (1.0d0 - 1.0d0/d2)*uzp 
      m2 = 1.0d0
    endif
    if( cidp3 == cid_target ) then
      uxn  = (1.0d0 - 1.0d0/d3)*uxp
      uyn  = (1.0d0 - 1.0d0/d3)*uyp
      uzn  = (1.0d0 - 1.0d0/d3)*uzp
      m3 = 1.0d0
    endif
    if( cidp4 == cid_target ) then
      uxb  = (1.0d0 - 1.0d0/d4)*uxp
      uyb  = (1.0d0 - 1.0d0/d4)*uyp
      uzb  = (1.0d0 - 1.0d0/d4)*uzp
      m4 = 1.0d0
    endif
    if( cidp5 == cid_target ) then
      uxt  = (1.0d0 - 1.0d0/d5)*uxp
      uyt  = (1.0d0 - 1.0d0/d5)*uyp
      uzt  = (1.0d0 - 1.0d0/d5)*uzp
      m5 = 1.0d0
    endif

    if( pidp /= 1 ) then
      m0 = 0.0d0
      m1 = 0.0d0
      m2 = 0.0d0
      m3 = 0.0d0
      m4 = 0.0d0
      m5 = 0.0d0
    endif

    fsvx(i, j, k) = + mu*((uxp - uxw))/dx*m0*2.0 &
                    - mu*((uxe - uxp))/dx*m1*2.0 &
                    + mu*((uxp - uxs) + (uye - uyw)*0.5)/dx*m2 &
                    - mu*((uxn - uxp) + (uye - uyw)*0.5)/dx*m3 &
                    + mu*((uxp - uxb) + (uze - uzw)*0.5)/dx*m4 &
                    - mu*((uxt - uxp) + (uze - uzw)*0.5)/dx*m5
    fsvy(i, j, k) = + mu*((uyp - uyw) + (uxn - uxs)*0.5)/dx*m0 &
                    - mu*((uye - uyp) + (uxn - uxs)*0.5)/dx*m1 &
                    + mu*((uyp - uys))/dx*m2*2.0 &
                    - mu*((uyn - uyp))/dx*m3*2.0 &
                    + mu*((uyp - uyb) + (uzn - uzs)*0.5)/dx*m4 &
                    - mu*((uyt - uyp) + (uzn - uzs)*0.5)/dx*m5
    fsvz(i, j, k) = + mu*((uzp - uzw) + (uxt - uxb)*0.5)/dx*m0 &
                    - mu*((uze - uzp) + (uxt - uxb)*0.5)/dx*m1 &
                    + mu*((uzp - uzs) + (uyt - uyb)*0.5)/dx*m2 &
                    - mu*((uzn - uzp) + (uyt - uyb)*0.5)/dx*m3 &
                    + mu*((uzp - uzb))/dx*m4*2.0 &
                    - mu*((uzt - uzp))/dx*m5*2.0

    fsx = fsx + fsvx(i, j, k)*dx*dx
    fsy = fsy + fsvy(i, j, k)*dx*dx
    fsz = fsz + fsvz(i, j, k)*dx*dx
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
  fsv(1) = fsx
  fsv(2) = fsy
  fsv(3) = fsz
end subroutine bcut_calc_f_v

subroutine bcut_calc_q( &
                qx, &
                qy, &
                qz, &
                q, &
                sa, &
                cid_target, &
                t0_, &
                c0, c1, c2, c3, c4, c5, &
                cid0, cid1, cid2, cid3, cid4, cid5, &
                pid, &
								nnum, &
								nx_, ny_, nz_, &
								nidx0, nidx1, nidx2, nidx3, nidx4, nidx5, &
                rhof, rhos, &
                cpf, cps, &
                kf, ks, &
                bc_n, &
                bc_type, &
                bc_value, &
                org, &
                dx, dt, &
                sz, g)
  implicit none
  integer                  :: i, j, k
  integer                  :: ix, jx, kx
  integer                  :: g
  integer, dimension(3)    :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: qx, qy, qz
  real, dimension(1:3)    :: q
  real                    :: sa
  integer                  :: cid_target
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: t0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
	integer									 :: nnum
	real, dimension(0:nnum-1) :: nx_, ny_, nz_
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: nidx0, nidx1, nidx2, nidx3, nidx4, nidx5
  integer                  :: cidp
  integer                  :: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
  integer                  :: pidp, pidw, pide, pids, pidn, pidb, pidt
  real                    :: rhof, rhos
  real                    :: cpf, cps
  real                    :: kf, ks
  integer                      :: bc_n
  integer, dimension(0:bc_n-1):: bc_type
  real, dimension(0:bc_n-1)   :: bc_value
  real                    :: dx, dt
  real, dimension(3)      :: org
  real                    :: d0, d1, d2, d3, d4, d5
  real                    :: k0, k1, k2, k3, k4, k5
  real                    :: m0, m1, m2, m3, m4, m5
  real                    :: l0, l1, l2, l3, l4, l5
  real                    :: bp, b0, b1, b2, b3, b4, b5
  real                    :: tp, t0, t1, t2, t3, t4, t5
  real                    :: nx, ny, nz
  real                    :: x, y, z, r2, rl, xi, yi, zi, theta, phi
  real                    :: qx0, qy0, qz0
  real                    :: qxt, qyt, qzt
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  qxt = 0.0
  qyt = 0.0
  qzt = 0.0
  sa = 0.0
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp           private(cidp) &
!$omp           private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp           private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp           private(d0, d1, d2, d3, d4, d5) &
!$omp           private(k0, k1, k2, k3, k4, k5) &
!$omp           private(m0, m1, m2, m3, m4, m5) &
!$omp           private(l0, l1, l2, l3, l4, l5) &
!$omp           private(bp, b0, b1, b2, b3, b4, b5) &
!$omp           private(tp, t0, t1, t2, t3, t4, t5) &
!$omp           private(nx, ny, nz) &
!$omp           private(qx0, qy0, qz0) &
!$omp           private(x, y, z, r2, rl, xi, yi, zi, theta, phi)
!$omp do schedule(static, 1), &
!$omp     reduction(+:qxt, qyt, qzt, sa)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
    d0 = c0(i, j, k)
    d1 = c1(i, j, k)
    d2 = c2(i, j, k)
    d3 = c3(i, j, k)
    d4 = c4(i, j, k)
    d5 = c5(i, j, k)

    cidp0 = cid0(i, j, k)
    cidp1 = cid1(i, j, k)
    cidp2 = cid2(i, j, k)
    cidp3 = cid3(i, j, k)
    cidp4 = cid4(i, j, k)
    cidp5 = cid5(i, j, k)

    pidp = pid(i, j, k)

    if( pidp /= 1 ) then
      cycle
    endif

    m0 = 0.0d0
    m1 = 0.0d0
    m2 = 0.0d0
    m3 = 0.0d0
    m4 = 0.0d0
    m5 = 0.0d0
    if( cidp0 == cid_target ) then
      m0 = 1.0d0
    endif
    if( cidp1 == cid_target ) then
      m1 = 1.0d0
    endif
    if( cidp2 == cid_target ) then
      m2 = 1.0d0
    endif
    if( cidp3 == cid_target ) then
      m3 = 1.0d0
    endif
    if( cidp4 == cid_target ) then
      m4 = 1.0d0
    endif
    if( cidp5 == cid_target ) then
      m5 = 1.0d0
    endif

    k0 = kf
    k1 = kf
    k2 = kf
    k3 = kf
    k4 = kf
    k5 = kf
    tp = t0_(i, j, k)
    t0 = t0_(i-1, j, k) 
    t1 = t0_(i+1, j, k) 
    t2 = t0_(i, j-1, k) 
    t3 = t0_(i, j+1, k) 
    t4 = t0_(i, j, k-1) 
    t5 = t0_(i, j, k+1) 
    b0 = 0.0
    b1 = 0.0
    b2 = 0.0
    b3 = 0.0
    b4 = 0.0
    b5 = 0.0

    if( cidp0 == cid_target ) then
      if( bc_type(cidp0) == 0 ) then
				nx = nx_( nidx0(i, j, k) )

        k0 = kf/d0*2.0/(d0 + d1)
        k1 = kf/d1*2.0/(d0 + d1)
        t0 = bc_value(cidp0)
        qx0 = -(tp - t0)/(d0*dx)

        qx(i, j, k) = qx0
        qxt = qxt + qx0*dx*dx
        sa = sa + abs(nx)*dx*dx
      else if( bc_type(cidp0) == 1 ) then
 				nx = nx_( nidx0(i, j, k) )

        k0 = 0.0
        k1 = kf/(d0 + 0.5)
        t0 = 0.0
        b0 = - nx*bc_value(cidp0)*dx/(d0 + 0.5)
        t0 = tp - nx*bc_value(cidp0)*d0*dx
        t0 = tp + abs(nx)*bc_value(cidp0)*d0*dx
        qx0 = 1.0/t0
        qx0 = t0

        qx(i, j, k) = qx0
        qxt = qxt + qx0*abs(nx)*dx*dx
        sa = sa + abs(nx)*dx*dx
      endif
    endif

    if( cidp1 == cid_target ) then
      if( bc_type(cidp1) == 0 ) then
				nx = nx_( nidx1(i, j, k) )

        k0 = kf/d0*2.0/(d0 + d1)
        k1 = kf/d1*2.0/(d0 + d1)
        t1 = bc_value(cidp1)
        qx0 = -(t1 - tp)/(d1*dx)

        qx(i, j, k) = qx0
        qxt = qxt - qx0*dx*dx
        sa = sa + abs(nx)*dx*dx
      else if( bc_type(cidp1) == 1 ) then
				nx = nx_( nidx1(i, j, k) )

        k0 = kf/(d1 + 0.5)
        k1 = 0.0
        t1 = 0.0
        b1 = + nx*bc_value(cidp1)*dx/(d1 + 0.5)
        t1 = tp + nx*bc_value(cidp0)*d1*dx
        t1 = tp + abs(nx)*bc_value(cidp0)*d1*dx
        qx0 = 1.0/t1
        qx0 = t1

        qx(i, j, k) = qx0
        qxt = qxt + qx0*abs(nx)*dx*dx
        sa = sa + abs(nx)*dx*dx
      endif
    endif

    if( cidp2 == cid_target ) then
      if( bc_type(cidp2) == 0 ) then
				ny = ny_( nidx2(i, j, k) )

        k2 = kf/d2*2.0/(d2 + d3)
        k3 = kf/d3*2.0/(d2 + d3)
        t2 = bc_value(cidp2)
        qy0 = -(tp - t2)/(d2*dx)

        qy(i, j, k) = qy0
        qyt = qyt + qy0*dx*dx
        sa = sa + abs(ny)*dx*dx
      else if( bc_type(cidp2) == 1 ) then
				ny = ny_( nidx2(i, j, k) )

        k2 = 0.0
        k3 = kf/(d2 + 0.5)
        t2 = 0.0
        b2 = - ny*bc_value(cidp2)*dx/(d2 + 0.5)
        t2 = tp - ny*bc_value(cidp2)*d2*dx
        t2 = tp + abs(ny)*bc_value(cidp2)*d2*dx
        qy0 = 1.0/t2
        qy0 = t2

        qy(i, j, k) = qy0
        qyt = qyt + qy0*abs(ny)*dx*dx
        sa = sa + abs(ny)*dx*dx
      endif
    endif

    if( cidp3 == cid_target ) then
      if( bc_type(cidp3) == 0 ) then
				ny = ny_( nidx3(i, j, k) )

        k2 = kf/d2*2.0/(d2 + d3)
        k3 = kf/d3*2.0/(d2 + d3)
        t3 = bc_value(cidp3)
        qy0 = -(t3 - tp)/(d3*dx)

        qy(i, j, k) = qy0
        qyt = qyt - qy0*dx*dx
        sa = sa + abs(ny)*dx*dx
      else if( bc_type(cidp3) == 1 ) then
				ny = ny_( nidx3(i, j, k) )

        k2 = kf/(d3 + 0.5)
        k3 = 0.0
        t3 = 0.0
        b3 = + ny*bc_value(cidp3)*dx/(d3 + 0.5)
        t3 = tp + ny*bc_value(cidp3)*d3*dx
        t3 = tp + abs(ny)*bc_value(cidp3)*d3*dx
				qy0 = 1.0/t3
				qy0 = t3

        qy(i, j, k) = qy0
        qyt = qyt + qy0*abs(ny)*dx*dx
        sa = sa + abs(ny)*dx*dx
      endif
    endif

    if( cidp4 == cid_target ) then
      if( bc_type(cidp4) == 0 ) then
				nz = nz_( nidx4(i, j, k) )

        k4 = kf/d4*2.0/(d4 + d5)
        k5 = kf/d5*2.0/(d4 + d5)
        t4 = bc_value(cidp4)
        qz0 = -(tp - t4)/(d4*dx)

        qz(i, j, k) = qz0
        qzt = qzt + qz0*dx*dx
        sa = sa + abs(nz)*dx*dx
      else if( bc_type(cidp4) == 1 ) then
				nz = nz_( nidx4(i, j, k) )

        k4 = 0.0
        k5 = kf/(d4 + 0.5)
        t4 = 0.0
        b4 = - nz*bc_value(cidp4)*dx/(d4 + 0.5)
        t4 = tp - nz*bc_value(cidp4)*d4*dx
        t4 = tp + abs(nz)*bc_value(cidp4)*d4*dx
        qz0 = 1.0/t4
        qz0 = t4

        qz(i, j, k) = qz0
        qzt = qzt + qz0*abs(nz)*dx*dx
        sa = sa + abs(nz)*dx*dx
      endif
    endif

    if( cidp5 == cid_target ) then
      if( bc_type(cidp5) == 0 ) then
				nz = nz_( nidx5(i, j, k) )

        k4 = kf/d4*2.0/(d4 + d5)
        k5 = kf/d5*2.0/(d4 + d5)
        t5 = bc_value(cidp5)
        qz0 = -(t5 - tp)/(d5*dx)

        qz(i, j, k) = qz0
        qzt = qzt - qz0*dx*dx
        sa = sa + abs(nz)*dx*dx
      else if( bc_type(cidp5) == 1 ) then
				nz = nz_( nidx5(i, j, k) )

        k4 = kf/(d5 + 0.5)
        k5 = 0.0
        t5 = 0.0
        b5 = + nz*bc_value(cidp5)*dx/(d5 + 0.5)
        t5 = tp + nz*bc_value(cidp5)*d5*dx
        t5 = tp + abs(nz)*bc_value(cidp5)*d5*dx
        qz0 = 1.0/t5
        qz0 = t5

        qz(i, j, k) = qz0
        qzt = qzt + qz0*abs(nz)*dx*dx
        sa = sa + abs(nz)*dx*dx
      endif
    endif
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
  q(1) = qxt
  q(2) = qyt
  q(3) = qzt
end subroutine bcut_calc_q

subroutine bcut_set_fluidseed( &
                pid, &
                xs, ys, zs, &
                dx, &
                org, &
                sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
  real                    :: xs, ys, zs
  real                    :: dx
  real, dimension(3)      :: org
  real                    :: x0, y0, z0
  real                    :: x1, y1, z1
  integer                  :: pidp
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp           private(x0, y0, z0) &
!$omp           private(x1, y1, z1) &
!$omp           private(pidp)
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
    x0 = org(1) + (real(i-1))*dx
    y0 = org(2) + (real(j-1))*dx
    z0 = org(3) + (real(k-1))*dx
    x1 = org(1) + (real(i))*dx
    y1 = org(2) + (real(j))*dx
    z1 = org(3) + (real(k))*dx

    if( (x0 <= xs .and. xs <= x1) .and. &
        (y0 <= ys .and. ys <= y1) .and. &
        (z0 <= zs .and. zs <= z1) ) then
      pid(i, j, k) = 1
!      write(*, *) i, j, k
!      write(*, *) x0, xs, x1
!      write(*, *) y0, ys, y1
!      write(*, *) z0, zs, z1
    endif
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_set_fluidseed

subroutine bcut_set_reference_value( &
                Ap, Aw, Ae, As, An, Ab, At, b, &
                xr, yr, zr, &
                pr, &
                dx, &
                org, &
                sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: Ap, Aw, Ae, As, An, Ab, At, b
  real                    :: xr, yr, zr
  real                    :: pr
  real                    :: dx
  real, dimension(3)      :: org
  real                    :: x0, y0, z0
  real                    :: x1, y1, z1
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp           private(x0, y0, z0) &
!$omp           private(x1, y1, z1) 
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
    x0 = org(1) + (real(i-1))*dx
    y0 = org(2) + (real(j-1))*dx
    z0 = org(3) + (real(k-1))*dx
    x1 = org(1) + (real(i))*dx
    y1 = org(2) + (real(j))*dx
    z1 = org(3) + (real(k))*dx

    if( (x0 <= xr .and. xr <= x1) .and. &
        (y0 <= yr .and. yr <= y1) .and. &
        (z0 <= zr .and. zr <= z1) ) then
      Ap(i, j, k) = 1.0d0
      Aw(i, j, k) = 0.0d0
      Ae(i, j, k) = 0.0d0
      As(i, j, k) = 0.0d0
      An(i, j, k) = 0.0d0
      Ab(i, j, k) = 0.0d0
      At(i, j, k) = 0.0d0
      b (i, j, k) = pr

      Ap(i+1, j, k) = Ap(i+1, j, k) - Aw(i+1, j, k)
			b (i+1, j, k) = b (i+1, j, k) - Aw(i+1, j, k)*pr*2.0
      Aw(i+1, j, k) = 0.0d0

      Ap(i-1, j, k) = Ap(i-1, j, k) - Ae(i+1, j, k)
			b (i-1, j, k) = b (i-1, j, k) - Ae(i-1, j, k)*pr*2.0
      Ae(i-1, j, k) = 0.0d0

      Ap(i, j+1, k) = Ap(i, j+1, k) - As(i, j+1, k)
			b (i, j+1, k) = b (i, j+1, k) - As(i, j+1, k)*pr*2.0
      As(i, j+1, k) = 0.0d0

      Ap(i, j-1, k) = Ap(i, j-1, k) - An(i, j-1, k)
			b (i, j-1, k) = b (i, j-1, k) - An(i, j-1, k)*pr*2.0
      An(i, j-1, k) = 0.0d0

      Ap(i, j, k+1) = Ap(i, j, k+1) - Ab(i, j, k+1)
			b (i, j, k+1) = b (i, j, k+1) - Ab(i, j, k+1)*pr*2.0
      Ab(i, j, k+1) = 0.0d0

      Ap(i, j, k-1) = Ap(i, j, k-1) - At(i, j, k-1)
			b (i, j, k-1) = b (i, j, k-1) - At(i, j, k-1)*pr*2.0
      At(i, j, k-1) = 0.0d0
    endif
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_set_reference_value

subroutine bcut_set_a( &
                ux, uy, uz, &
                p, &
                pid, &
                rho, &
                mu, &
                a, &
                U, &
                dx, dt, &
                org, &
                sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: ux, uy, uz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: p
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
  real                    :: rho
  real                    :: mu
  real                    :: a
  real                    :: U
  real                    :: dx, dt
  real, dimension(3)      :: org
  real                    :: x, y, z
  real                    :: r
  integer                  :: pidp
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp           private(x, y, z) &
!$omp           private(r) &
!$omp           private(pidp)
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
    x = org(1) + (real(i) - 0.5)*dx
    y = org(2) + (real(j) - 0.5)*dx
    z = org(3) + (real(k) - 0.5)*dx
    r = sqrt(x*x + y*y + z*z)
    pidp = pid(i, j, k)
!    if( pidp == 1 ) then
    if( r > 0.3 ) then
      p(i, j, k) = -1.5*mu*a*U*x/(r*r*r)
      ux(i, j, k) = U - 0.25*a*U/r*(3.0 + a*a/(r*r)) - 0.75*a*U*x*x/(r*r*r)*(1.0 - a*a/(r*r))
      uy(i, j, k) =   - 0.75*a*U*x*y/(r*r*r)*(1.0 - a*a/(r*r))
      uz(i, j, k) =   - 0.75*a*U*x*z/(r*r*r)*(1.0 - a*a/(r*r))
    endif
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_set_a

