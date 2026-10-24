module poly_dec
  ! define length of array parameters
  integer, parameter :: n_e_max=101

  ! Footprint edges longer than edge_long [rad] are split by densify_polygon
  ! so that the great-circle segments stay within edge_tol [rad] of the
  ! straight RA/Dec line between the original vertices.
  real (kind=8), parameter :: edge_long = 1.1d0*3.141592653589793238d0/180.d0
  real (kind=8), parameter :: edge_tol = &
       2.d0/3600.d0*3.141592653589793238d0/180.d0

  ! define derived types
  type t_polygon
     integer :: n
     real (kind=8), dimension(n_e_max) :: x, y
  end type t_polygon

end module poly_dec
