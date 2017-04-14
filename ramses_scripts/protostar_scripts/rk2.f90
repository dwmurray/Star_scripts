
subroutine derivs( x, n, y, dydx) 
  implicit none
  integer, parameter :: nvars = 2
  double precision, intent(in) :: x, n 
  double precision, intent(in), dimension(nvars) :: y
  double precision, intent(out), dimension(nvars) :: dydx
  double precision :: eps, dphideps, phi

  phi = y(1)
  dphideps = y(2)
  eps = x

  dydx(1) = dphideps
  dydx(2) = -phi**n - 2./eps*dphideps  

end subroutine derivs

subroutine rk2( func, x1, x2, n, y1, y2) 
  implicit none
  integer, parameter :: nvars = 2
  double precision, intent(in) :: x1, x2, n
  double precision, intent(in), dimension(nvars) :: y1
  double precision, intent(out), dimension(nvars) :: y2
  external :: func

  double precision, dimension(nvars) :: yhalf, dydx
  double precision :: dx

  dx = x2-x1

  ! compute the half step
  call func( x1, n, y1, dydx)

  yhalf = y1 + dx*0.5*dydx

  ! compute the full step
  call func( x1+dx*0.5, n, yhalf, dydx)

  y2 = y1 + dx*dydx
  
  return
end subroutine rk2


module rk2_polytrope
  implicit none
  PUBLIC :: polytrope
!
contains
!
subroutine polytrope( n, eps1, dphideps1, phi1)
  implicit none
  integer, parameter :: nvars = 2
  double precision, intent(in) :: n
  double precision :: phi, dphideps, eps, deps
  double precision, dimension(nvars) :: y1, y2, dydx
  external :: derivs, rk2
  double precision, intent(out) :: eps1, dphideps1, phi1
!  integer, intent(out) :: array_length

!  array_length = 0

  ! set initial conditions
!  eps = 0.2D0
  eps = 1.0D-5
  deps = 0.00001D0
  phi = 1. - eps*eps/6. + n/120.*eps**4.
  dphideps = -eps/3. + n/30.*eps**3.

  do while( phi > 1e-7) 
     y1(1) = phi
     y1(2) = dphideps

     call rk2( derivs, eps, eps+deps, n, y1, y2)
     !calculate the next step

     eps = eps+deps
     phi = y2(1)
     dphideps = y2(2)
     deps = max(1.e-6, min(0.01, abs(phi/dphideps)*0.1))
     !call derivs(eps, n, y2, dydx)
     !deps = max(0.0001,0.1*min(abs(y2(1)/dydx(1)), abs(y2(2)/dydx(2))))
     !Pass out phi, eps, dphideps
     eps1 = eps
     dphideps1 = dphideps
     phi1 = phi

!     array_length = array_length + 1
  enddo

!  write(*,*) array_length
!  write(*,*) 'eps1 is', eps1(array_length)
!  write(*,*) 'eps1(1) is', eps1(1), 'eps1(2) is ', eps1(2)
  return

end subroutine polytrope


!
subroutine phi_of_r( n, eps, deps, phi, dphideps, eps_r, dphideps_r, phi_r)
  implicit none
  integer, parameter :: nvars = 2
  double precision, intent(in) :: n, eps, deps, phi, dphideps
!  double precision :: phi, dphideps, eps, deps
  double precision, dimension(nvars) :: y1, y2, dydx
  external :: derivs, rk2
  double precision, intent(out) :: eps_r, dphideps_r, phi_r
!  integer, intent(in) :: maxiter
  integer :: i

  y1(1) = phi
  y1(2) = dphideps
!  write(*,*) n, eps, deps, phi, dphideps
  call rk2( derivs, eps, eps+deps, n, y1, y2)

  !Pass out phi, eps, dphideps
  eps_r = eps + deps
  phi_r =y2(1) !phi
  dphideps_r = y2(2)!dphideps
!  write(*,*) 'eps_r is', eps_r


!  do i=1, maxiter !while( phi > 1e-7) 
!     y1(1) = phi
!     y1(2) = dphideps
!
!     call rk2( derivs, eps, eps+deps, n, y1, y2)
!     !calculate the next step
!
!     eps = eps+deps
!     phi = y2(1)
!     dphideps = y2(2)
!!     deps = max(1.e-6, min(0.01, abs(phi/dphideps)*0.1))
!
!     !Pass out phi, eps, dphideps
!     eps_r = eps
!     dphideps_r = dphideps
!     phi_r = phi
!     write(*,*) 'eps_r is', eps_r
!  end do


  return

end subroutine phi_of_r

end module rk2_polytrope
!
!program test_rk2
!  external :: polytrope
!  double precision :: n, eps1, dphideps1
!
!  !do n = 0.5, 5.0, 0.05
!  do n = 1.4, 3.1, 0.05
!     call polytrope( n, eps1, dphideps1) 
!     write (*,*) n, eps1, -eps1**2*dphideps1
!  end do
!
!  stop
!end program test_rk2
! 
!
