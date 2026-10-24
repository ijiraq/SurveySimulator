program test_cosdec_footprint
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! Load data/cosdec_survey/pointings.list and check that the footprint
! area and the extent tested by point_in_polygon do not depend on
! declination, i.e. that the cos(dec) term is applied when the
! footprints are built. Also checks RA wrap, a far-away point, and the
! survey.conf edges setting: long edges follow the RA/Dec line by default
! (data/cosdec_survey) and are single great circles with
! edges = great_circle (data/cosdec_gc).
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

  use parameters
  use datadec
  use poly_lib
  use getsur

  implicit none

  integer, parameter :: n_expected = 10, i_ears = 7, i_strip = 9, &
       i_slant = 10
  type(t_pointing), allocatable :: points(:)
  real (kind=8), allocatable :: sur_mm(:)
  real (kind=8) :: ra0(n_expected), dec0(n_expected)
  real (kind=8) :: area, expected, w, h, e_w, e_h, p(2), f, dra_t
  real (kind=8) :: sx(4), sy(4), vx(4), vy(4), off
  integer :: n_sur, ierr, nfail, i, k
  integer :: n_epochs, epoch_n(max_epochs), epoch_off(max_epochs), det_req
  character(name_len) :: sname

  data ra0  / 30.d0, 60.d0, 90.d0, 120.d0, 150.d0, 0.1d0, 200.d0, 250.d0, &
       300.d0, 213.87d0 /
  data dec0 /  0.d0, 30.d0, 60.d0, -60.d0,  75.d0, 30.d0,  60.d0,  60.d0, &
       60.d0, -12.54d0 /
  data sx / -3.5d0, -3.5d0, 3.5d0, 3.5d0 /
  data sy / -0.434889d0, 2.565111d0, 0.473417d0, -2.526583d0 /

  nfail = 0
  allocate (points(n_sur_max), sur_mm(n_sur_max))
  call GetSurvey('data/cosdec_survey', 21, n_sur, points, sur_mm, ierr)
  if ((ierr .ne. 0) .or. (n_sur .ne. n_expected)) then
     write (6, *) 'FAIL: GetSurvey ierr=', ierr, ' n_sur=', n_sur
     stop 1
  end if
  write (6, *) 'PASS: loaded n_sur=', n_sur

! Vertex counts: single fields and the ears keep their vertices, edges
! longer than ~1 deg are densified
  do i = 1, n_expected
     if (i .eq. i_ears) then
        call check_count('ears vertices', i, points(i)%poly%n, 12, 12)
     else if ((i .eq. i_strip) .or. (i .eq. i_slant)) then
        call check_count('densified vertices', i, points(i)%poly%n, 5, &
             n_e_max - 1)
     else
        call check_count('vertices', i, points(i)%poly%n, 4, 4)
     end if
  end do

! Area: 1 deg x 1 deg at every declination
  do i = 1, n_expected
     if ((i .eq. i_ears) .or. (i .eq. i_strip) .or. (i .eq. i_slant)) cycle
     area = polygon_area(points(i)%poly)
     call check_close('area', i, area, 1.d0, 1.d-3)
  end do
  area = polygon_area(points(i_strip)%poly)
  call check_close('strip area', i_strip, area, 8.d0, 1.d-3)

! Ears area from the create_ears constants
  w = 0.008545d0
  h = 0.008678d0
  e_w = 0.001913d0
  e_h = 0.004169d0
  expected = (4.d0*w*h + 4.d0*e_w*e_h)*(180.d0/Pi)**2
  area = polygon_area(points(i_ears)%poly)
  call check_close('ears area', i_ears, area, expected, 1.d-3)

! Membership from tangent-plane offsets [deg] about each centre
  do i = 1, n_expected
     if ((i .eq. i_ears) .or. (i .eq. i_strip) .or. (i .eq. i_slant)) cycle
     call check_in(i,  0.49d0,  0.d0,   .true.)
     call check_in(i, -0.49d0,  0.d0,   .true.)
     call check_in(i,  0.d0,    0.49d0, .true.)
     call check_in(i,  0.d0,   -0.49d0, .true.)
     call check_in(i,  0.49d0,  0.49d0, .true.)
     call check_in(i,  0.51d0,  0.d0,   .false.)
     call check_in(i, -0.51d0,  0.d0,   .false.)
     call check_in(i,  0.d0,    0.51d0, .false.)
     call check_in(i,  0.d0,   -0.51d0, .false.)
     call check_in(i,  0.51d0,  0.49d0, .false.)
  end do

! Ears: half-width 0.4896 deg, ears out to 0.5992 deg, ear half-height
! 0.2389 deg, half-height 0.4972 deg
  call check_in(i_ears,  0.59d0,  0.d0,   .true.)
  call check_in(i_ears, -0.59d0,  0.d0,   .true.)
  call check_in(i_ears,  0.61d0,  0.d0,   .false.)
  call check_in(i_ears, -0.61d0,  0.d0,   .false.)
  call check_in(i_ears,  0.48d0,  0.45d0, .true.)
  call check_in(i_ears,  0.50d0,  0.45d0, .false.)
  call check_in(i_ears,  0.d0,    0.49d0, .true.)
  call check_in(i_ears,  0.d0,   -0.51d0, .false.)

! Strip: top and bottom edges follow Dec = 60 +/- 0.5 along the whole
! strip, including midway between field boundaries. A single great
! circle between the strip corners would bulge ~0.25 deg poleward.
  dra_t = 4.d0/dcos((dec0(i_strip) + 0.5d0)*drad)
  do k = 0, 7
     f = (dble(k) + 0.5d0)/4.d0 - 1.d0
     p(1) = (ra0(i_strip) + f*dra_t)*drad
     p(2) = (dec0(i_strip) + 0.49d0)*drad
     call check_point('strip top in', i_strip, p, .true.)
     p(2) = (dec0(i_strip) + 0.51d0)*drad
     call check_point('strip top out', i_strip, p, .false.)
     p(2) = (dec0(i_strip) - 0.49d0)*drad
     call check_point('strip bottom in', i_strip, p, .true.)
     p(2) = (dec0(i_strip) - 0.51d0)*drad
     call check_point('strip bottom out', i_strip, p, .false.)
  end do

! Slanted block: the long top (2->3) and bottom (4->1) edges follow the
! straight RA/Dec line between the vertices, to 10 arcsec.
  do k = 1, 4
     vy(k) = dec0(i_slant) + sy(k)
     vx(k) = ra0(i_slant) + sx(k)/dcos(vy(k)*drad)
  end do
  off = 10.d0/3600.d0
  do k = 1, 3
     f = 0.25d0*dble(k)
     p(1) = (vx(2) + f*(vx(3) - vx(2)))*drad
     p(2) = (vy(2) + f*(vy(3) - vy(2)) - off)*drad
     call check_point('slant top in', i_slant, p, .true.)
     p(2) = (vy(2) + f*(vy(3) - vy(2)) + off)*drad
     call check_point('slant top out', i_slant, p, .false.)
     p(1) = (vx(4) + f*(vx(1) - vx(4)))*drad
     p(2) = (vy(4) + f*(vy(1) - vy(4)) + off)*drad
     call check_point('slant bottom in', i_slant, p, .true.)
     p(2) = (vy(4) + f*(vy(1) - vy(4)) - off)*drad
     call check_point('slant bottom out', i_slant, p, .false.)
  end do

! RA wrap: RA 359.8 deg is 0.3 deg of RA west of the RA 0.1 deg field
  p(1) = 359.8d0*drad
  p(2) = 30.d0*drad
  call check_point('RA wrap', 6, p, .true.)

! Far point: 120 deg from the centre of pointing 1
  p(1) = 150.d0*drad
  p(2) = 0.d0
  call check_point('far point', 1, p, .false.)

! survey.conf edges = great_circle: the same strip keeps 4 vertices, its
! top edge is one great circle bulging ~0.25 deg poleward, and the
! setting does not leak into later direct GetSurvey calls.
  call LoadSurveyRoot('data/cosdec_gc', 21, n_sur, points, sur_mm, &
       n_epochs, epoch_n, epoch_off, det_req, sname, ierr)
  if ((ierr .ne. 0) .or. (n_sur .ne. 1)) then
     write (6, *) 'FAIL: LoadSurveyRoot cosdec_gc ierr=', ierr, ' n_sur=', n_sur
     nfail = nfail + 1
  else
     call check_count('great_circle vertices', 1, points(1)%poly%n, 4, 4)
     area = polygon_area(points(1)%poly)
     if (area .lt. 8.03d0) then
        write (6, *) 'FAIL: great_circle strip area', area, ' expected ~8.05'
        nfail = nfail + 1
     else
        write (6, *) 'PASS: great_circle strip area =', area
     end if
     p(1) = ra0(i_strip)*drad
     p(2) = (dec0(i_strip) + 0.51d0)*drad
     call check_point('great_circle bulge', 1, p, .true.)
  end if
  if (.not. radec_edges) then
     write (6, *) 'FAIL: radec_edges not reset after LoadSurveyRoot'
     nfail = nfail + 1
  else
     write (6, *) 'PASS: radec_edges reset after LoadSurveyRoot'
  end if

  if (nfail .eq. 0) then
     write (6, *) 'ALL TESTS PASSED'
     stop 0
  else
     write (6, *) 'FAILURES:', nfail
     stop 1
  end if

contains

  subroutine check_count(label, i, n, nmin, nmax)
    character(*), intent(in) :: label
    integer, intent(in) :: i, n, nmin, nmax

    if ((n .lt. nmin) .or. (n .gt. nmax)) then
       write (6, *) 'FAIL: ', label, ' pointing', i, ' n=', n, &
            ' expected', nmin, ' to', nmax
       nfail = nfail + 1
    else
       write (6, *) 'PASS: ', label, ' pointing', i, ' n=', n
    end if
  end subroutine check_count

  subroutine check_close(label, i, val, ref, tol)
    character(*), intent(in) :: label
    integer, intent(in) :: i
    real (kind=8), intent(in) :: val, ref, tol

    if (dabs(val - ref) .gt. tol*dabs(ref)) then
       write (6, *) 'FAIL: ', label, ' pointing', i, ' got', val, &
            ' expected', ref
       nfail = nfail + 1
    else
       write (6, *) 'PASS: ', label, ' pointing', i, ' =', val
    end if
  end subroutine check_close

  subroutine check_in(i, xi_deg, eta_deg, inside)
! Deproject the tangent-plane offset (east, north) about the centre of
! pointing i and test it against that pointing's polygon.
    integer, intent(in) :: i
    real (kind=8), intent(in) :: xi_deg, eta_deg
    logical, intent(in) :: inside
    real (kind=8) :: a, d, xi, eta, u(3), q(2)
    character(len=64) :: label

    a = ra0(i)*drad
    d = dec0(i)*drad
    xi = xi_deg*drad
    eta = eta_deg*drad
    u(1) = dcos(d)*dcos(a) - xi*dsin(a) - eta*dsin(d)*dcos(a)
    u(2) = dcos(d)*dsin(a) + xi*dcos(a) - eta*dsin(d)*dsin(a)
    u(3) = dsin(d) + eta*dcos(d)
    q(1) = modulo(datan2(u(2), u(1)), TwoPi)
    q(2) = dasin(u(3)/dsqrt(u(1)**2 + u(2)**2 + u(3)**2))
    write (label, '(a,f6.2,a,f6.2,a)') 'offset (', xi_deg, ',', eta_deg, ')'
    call check_point(label(1:len_trim(label)), i, q, inside)
  end subroutine check_in

  subroutine check_point(label, i, q, inside)
    character(*), intent(in) :: label
    integer, intent(in) :: i
    real (kind=8), intent(in) :: q(2)
    logical, intent(in) :: inside
    integer :: res
    logical :: got

    res = point_in_polygon(q, points(i)%poly)
    got = res .gt. 1
    if ((got .neqv. inside) .or. (res .eq. 1)) then
       write (6, *) 'FAIL: ', label, ' pointing', i, ' dec', dec0(i), &
            ' inside expected', inside, ' point_in_polygon=', res
       nfail = nfail + 1
    else
       write (6, *) 'PASS: ', label, ' pointing', i, ' dec', dec0(i), &
            ' inside=', got
    end if
  end subroutine check_point

end program test_cosdec_footprint
