module poly_lib

  use poly_dec
  use parameters

contains
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
! File Polygon-lib.f
!
! J.-M. Petit  Observatoire de Besacon
! Version 1 :  February 2016
!
! The purpose of this library is to provide polygon-oriented routines.
! The first and most important one is:
!     point_in_polygon (p, poly, n)
! which tells if the point "p" is inside, outside or touching the
! polygon "poly". Points and vertices are RA/Dec [rad]; polygon edges are
! great circles, and RA wrap at 0/2pi is handled.
!
! Also provides:
!     polygon_area (poly)
! spherical area of a simple RA/Dec footprint, in square degrees, and
!     densify_polygon (poly)
! which makes long edges follow the straight RA/Dec line.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

  integer function point_in_polygon(p, poly)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This function checks if point "p" is inside the polygon "poly" using
! the quadrant method.
!
! The polygon vertices are gnomonically projected about the given point,
! which puts the point at the origin and maps the great-circle edges onto
! straight lines. Then for each edge of the polygon, calc_walk_summand is
! called. If the sum of all returned values from these calls is +4 or
! -4, the point lies indeed inside the polygon. Otherwise, if a
! PolygonsTouching exception has been caught, the point lies on one of
! he edges of the polygon.
!
! Returns the number of nodes of the polygon, if the point lies inside,
! otherwise 1 if the point lies on the polygon and if not, 0.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J.-M. Petit  Observatoire de Besancon
! Version 1 : February 2016
! Version 2 : October 2026
!             Great-circle edges via gnomonic projection about p
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     p     : Point (array (ra,dec)) [rad] (2*R8)
!     poly  : Polygon structure (polygon), RA/Dec [rad]
!
! OUTPUT
!     point_in_polygon: result of the call (I4)
!                Point inside polygon : n
!                Point on polygon     : 1
!                Point outside polygon: 0
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! Set of F2PY directive to create a Python module
!
!f2py intent(in) p
!f2py intent(in) poly
    implicit none

    type(t_polygon), intent(in) :: poly
    real (kind=8), intent(in) :: p(2)
    integer :: i, walk_sum, walk
    real (kind=8) :: moved(2,n_e_max)
    real (kind=8) :: sd0, cd0, dra, sy, cy, cosc

! Gnomonic projection about p. A vertex 90 deg or more away cannot belong
! to a footprint containing p (footprints are smaller than a hemisphere).
    sd0 = dsin(p(2))
    cd0 = dcos(p(2))
    do i = 1, poly%n+1
       dra = poly%x(i) - p(1)
       sy = dsin(poly%y(i))
       cy = dcos(poly%y(i))
       cosc = sd0*sy + cd0*cy*dcos(dra)
       if (cosc .le. 0.d0) then
          point_in_polygon = 0
          return
       end if
       moved(1,i) = cy*dsin(dra)/cosc
       moved(2,i) = (cd0*sy - sd0*cy*dcos(dra))/cosc
    end do

! Computing the running sum
    walk_sum = 0
    do i = 1, poly%n
       walk = calc_walk_summand(moved(1,i), moved(1,i+1))
       if (walk .eq. -100) then
! Point is touching the polygon
          point_in_polygon = 1
          return
       end if
! Point is not on polygon
       walk_sum = walk_sum + walk
    end do

! Final check
    if (abs(walk_sum) .eq. 4) then
       point_in_polygon = poly%n
    else
       point_in_polygon = 0
    end if
    return

  end function point_in_polygon

  integer function calc_walk_summand(p1, p2)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This function calculates the summand along one edge depending on axis
! crossings.
!
! Follows the edge between two points and checks if one or both axes
! are being crossed. If They are crossed in clockwise sense, it returns
! +1 otherwise -1. Going through the origin raises the PolygonsTouching
! exception (returns -100).
!
! Returns one of -2, -1, 0, +1, +2 or raises PolygonsTouching
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J.-M. Petit  Observatoire de Besancon
! Version 1 : February 2016
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     p1    : First point (array (x,y)) of edge (2*R8)
!     p2    : Second point of edge (2*R8)
!
! OUTPUT
!     calc_walk_summand: result of the call (I4)
!                Clockwise crossing        : +1
!                No crossing               : 0
!                Counter-clockwise crossing: -1
!                Diagonal crossing         : +/-2
!                Point on edge             : -100
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! Set of F2PY directive to create a Python module
!
!f2py intent(in) p1
!f2py intent(in) p2
    implicit none

    real (kind=8), intent(in) :: p1(2), p2(2)
! Indices for better readability
    integer, parameter :: x=1, y=2
    integer :: summand
    real (kind=8) :: tx, ty, x_y0, y_x0

    summand = 0
! Here, we assume the 2 points are different!
!
! Checking for vertical line
    if (p1(x) .ne. p2(x)) then
       ty = p1(x)/(p1(x) - p2(x))
    else
       ty = p1(y)/(p1(y) - p2(y))
    end if

! Checking for horizontal line
    if (p1(y) .ne. p2(y)) then
       tx = p1(y)/(p1(y) - p2(y))
    else
       tx = ty
    end if

! Compute position of axis intersection
    x_y0 = p1(x) + tx*(p2(x) - p1(x))
    y_x0 = p1(y) + ty*(p2(y) - p1(y))

! Check if crossing x axis
    if ((tx .ge. 0.d0) .and. (tx .lt. 1.d0)) then
! Check if origin on edge
       if ((x_y0 .eq. 0.d0) .and. (y_x0 .eq. 0.d0)) then
          calc_walk_summand = -100
          return
       end if
       x_y0 = x_y0*(p2(y) - p1(y))
       if (x_y0 .ne. 0.d0) summand = summand + ceiling(sign(1.d0, x_y0))
    end if

! Check if crossing y axis
    if ((ty .ge. 0.d0) .and. (ty .lt. 1.d0)) then
! Check if origin on edge
       if ((x_y0 .eq. 0.d0) .and. (y_x0 .eq. 0.d0)) then
          calc_walk_summand = -100
          return
       end if
       y_x0 = y_x0*(p1(x) - p2(x))
       if (y_x0 .ne. 0.d0) summand = summand + ceiling(sign(1.d0, y_x0))
    end if
    calc_walk_summand = summand
    return

  end function calc_walk_summand

  real (kind=8) function polygon_area(poly)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! Spherical area of a simple closed footprint on the sky.
!
! Vertices are RA (x) and Dec (y) in radians, as stored after read_sur.
! The polygon is fan-triangulated from the first vertex; each spherical
! triangle area uses the van Oosterom & Strackee formula on unit vectors.
! Result is absolute area in square degrees.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! Version 1 : September 2026
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     poly  : Polygon structure (polygon), RA/Dec [rad]
!
! OUTPUT
!     polygon_area : area in square degrees (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) poly
    implicit none

    type(t_polygon), intent(in) :: poly
    integer :: i, n
    real (kind=8) :: area_sr, tri
    real (kind=8) :: v1(3), v2(3), v3(3)

    polygon_area = 0.d0
    n = poly%n
    if (n .lt. 3) return

    call ra_dec_to_xyz(poly%x(1), poly%y(1), v1)
    area_sr = 0.d0
    do i = 2, n - 1
       call ra_dec_to_xyz(poly%x(i), poly%y(i), v2)
       call ra_dec_to_xyz(poly%x(i+1), poly%y(i+1), v3)
       tri = spherical_triangle_area(v1, v2, v3)
       area_sr = area_sr + tri
    end do
    polygon_area = dabs(area_sr)*(180.d0/Pi)**2
    return
  end function polygon_area

  subroutine ra_dec_to_xyz(ra, dec, v)
! Convert RA/Dec [rad] to a unit vector.
!f2py intent(in) ra
!f2py intent(in) dec
!f2py intent(out) v
    implicit none
    real (kind=8), intent(in) :: ra, dec
    real (kind=8), intent(out) :: v(3)
    real (kind=8) :: cd

    cd = dcos(dec)
    v(1) = cd*dcos(ra)
    v(2) = cd*dsin(ra)
    v(3) = dsin(dec)
    return
  end subroutine ra_dec_to_xyz

  real (kind=8) function spherical_triangle_area(a, b, c)
! Signed spherical triangle area [sr] for unit vectors a,b,c.
! van Oosterom & Strackee: 2*atan2(det, 1 + a·b + b·c + c·a)
!f2py intent(in) a
!f2py intent(in) b
!f2py intent(in) c
    implicit none
    real (kind=8), intent(in) :: a(3), b(3), c(3)
    real (kind=8) :: det, denom

    det = a(1)*(b(2)*c(3) - b(3)*c(2)) &
         + a(2)*(b(3)*c(1) - b(1)*c(3)) &
         + a(3)*(b(1)*c(2) - b(2)*c(1))
    denom = 1.d0 + dot3(a, b) + dot3(b, c) + dot3(c, a)
    if ((det .eq. 0.d0) .and. (denom .eq. 0.d0)) then
       spherical_triangle_area = 0.d0
    else
       spherical_triangle_area = 2.d0*datan2(det, denom)
    end if
    return
  end function spherical_triangle_area

  real (kind=8) function angle3(u, v)
! Angle [rad] between unit vectors u and v.
    implicit none
    real (kind=8), intent(in) :: u(3), v(3)
    real (kind=8) :: c(3)
    c(1) = u(2)*v(3) - u(3)*v(2)
    c(2) = u(3)*v(1) - u(1)*v(3)
    c(3) = u(1)*v(2) - u(2)*v(1)
    angle3 = datan2(dsqrt(dot3(c, c)), dot3(u, v))
    return
  end function angle3

  subroutine densify_polygon(poly)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! Split every edge longer than edge_long into short great-circle segments
! whose vertices lie on the straight RA/Dec line between the original
! vertices, so that point_in_polygon and polygon_area follow that line to
! within edge_tol. Shorter edges stay single great circles.
!
! A great circle's departure from the RA/Dec line grows as the square of
! its length, so an edge whose midpoint departs by dev needs
! ceiling(sqrt(dev/edge_tol)) segments.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT/OUTPUT
!     poly  : Closed polygon structure (polygon), RA/Dec [rad]
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in,out) poly
    implicit none

    type(t_polygon), intent(inout) :: poly
    type(t_polygon) :: dense
    integer :: nseg(n_e_max), i, j, k, n, total
    real (kind=8) :: dev(n_e_max), tol, u1(3), u2(3), g(3), m(3)
    logical, save :: warned = .false.

    n = poly%n
    if (n .lt. 2) return
    do i = 1, n
       call ra_dec_to_xyz(poly%x(i), poly%y(i), u1)
       call ra_dec_to_xyz(poly%x(i+1), poly%y(i+1), u2)
       dev(i) = -1.d0
       if (angle3(u1, u2) .gt. edge_long) then
          g = u1 + u2
          g = g/dsqrt(dot3(g, g))
          call ra_dec_to_xyz(0.5d0*(poly%x(i) + poly%x(i+1)), &
               0.5d0*(poly%y(i) + poly%y(i+1)), m)
          dev(i) = angle3(g, m)
       end if
    end do

    tol = edge_tol
100 continue
    total = 0
    do i = 1, n
       nseg(i) = 1
       if (dev(i) .gt. tol) nseg(i) = ceiling(dsqrt(dev(i)/tol))
       total = total + nseg(i)
    end do
    if (total .gt. n_e_max - 1) then
       if (.not. warned) then
          write (6, *) 'densify_polygon: too many vertices, loosening ', &
               'edge tolerance beyond', edge_tol*206264.8062d0, ' arcsec'
          warned = .true.
       end if
       tol = 2.d0*tol
       goto 100
    end if
    if (total .eq. n) return

    k = 0
    do i = 1, n
       do j = 0, nseg(i) - 1
          k = k + 1
          dense%x(k) = poly%x(i) + (poly%x(i+1) - poly%x(i))*dble(j)/dble(nseg(i))
          dense%y(k) = poly%y(i) + (poly%y(i+1) - poly%y(i))*dble(j)/dble(nseg(i))
       end do
    end do
    dense%n = k
    dense%x(k+1) = dense%x(1)
    dense%y(k+1) = dense%y(1)
    poly = dense
    return
  end subroutine densify_polygon

  real (kind=8) function dot3(u, v)
    implicit none
    real (kind=8), intent(in) :: u(3), v(3)
    dot3 = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)
    return
  end function dot3

  subroutine check_polygon(poly)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! NOTE: user must make sure the input array is at least (2,n+1) long.
! This subroutine will copy first point into index n+1, if last point
! not already same as first. It will also check that there are no two
! consecutive points the same.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J.-M. Petit  Observatoire de Besancon
! Version 1 : February 2016
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     poly  : Polygon structure (polygon)
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in,out) poly
    implicit none

    type(t_polygon), intent(inout) :: poly
    integer :: i, j
 
    if ((poly%x(poly%n) .eq. poly%x(1)) .and. (poly%y(poly%n) .eq. poly%y(1))) &
         then
       poly%n = poly%n - 1
    else
       poly%x(poly%n+1) = poly%x(1)
       poly%y(poly%n+1) = poly%y(1)
    end if
    j = 1
1000 continue
    if ((poly%x(j) .eq. poly%x(j+1)) .and. (poly%y(j) .eq. poly%y(j+1))) then
       do i = j+1, poly%n
          poly%x(i) = poly%x(i+1)
          poly%y(i) = poly%y(i+1)
       end do
       poly%n = poly%n - 1
       j = j - 1
    end if
    j = j + 1
    if (j .le. poly%n) goto 1000
    return

  end subroutine check_polygon
end module poly_lib
