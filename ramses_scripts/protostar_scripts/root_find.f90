module root_find
!------------------------------------------------------------------------------------------
! Find the root, provided the function and derivative.
!------------------------------------------------------------------------------------------
  use constants, only            : r_15
  use protostellar_core
  implicit none
  PUBLIC Beta_root_find
  PUBLIC T_core_root_find
contains

  subroutine Beta_root_find(n, rho_core, P_core, xinit, machine_tolerance, maxiter, result)
    double precision, intent(in)  :: n, rho_core, P_core, xinit, machine_tolerance
    integer, intent(in)           :: maxiter
    double precision, intent(out) :: result
    double precision              :: Beta_func, Beta_func_prime
    double precision              :: x, xnew
    integer                       :: i
    
    result = 0.0
    x = xinit
    do i = 1, maxiter
!       xnew = x - (Beta_func(x) / Beta_func_prime(x))
       xnew = x - (Beta_func(x, rho_core, P_core) / Beta_func_prime(x, rho_core, P_core))
       if (abs(xnew-x) .le.  machine_tolerance) then
          result = xnew
          return
       end if
       x = xnew
    end do
    write(*,*) 'failure.'
  end subroutine Beta_root_find


  subroutine T_core_root_find(n, rho_core, P_core, xinit, machine_tolerance, maxiter, result)
    double precision, intent(in)  :: n, rho_core, P_core, xinit, machine_tolerance
    integer, intent(in)           :: maxiter
    double precision, intent(out) :: result
    double precision              :: T_core_func, T_core_func_prime
    double precision              :: x, xnew
    integer                       :: i
    
    result = 0.0
    x = xinit
    do i = 1, maxiter
!       xnew = x - (T_core_func(x) / T_core_func_prime(x))
       xnew = x - (T_core_func(x, rho_core, P_core) / T_core_func_prime(x, rho_core, P_core))
       if (abs(xnew-x) .le.  machine_tolerance) then
          result = xnew
          return
       end if
       x = xnew
    end do
    write(*,*) 'failure.'
  end subroutine T_core_root_find
end module root_find

function Beta_func(x, rho_core, P_core)
  use constants
  implicit none

  double precision             :: Beta_func
  double precision, intent(in) :: x, rho_core, P_core
  double precision, parameter  :: mu = 0.613 ! Offner et al 2009

!  write(*,*) a_rad, mu, m_p, k_b, rho_core, P_core
  !Beta_core !a_rad, k_b, m_p are defined in constants.f90
  Beta_func = (P_core**3 * a_rad / rho_core**4 / 3. * (mu * m_p / k_b)**4) * x**4 + x - 1.
  !Test function
  !func = x*x - 4.
end function Beta_func

function Beta_func_prime(x, rho_core, P_core)
  use constants
  implicit none

  double precision, intent(in) :: x, rho_core, P_core
  double precision             :: Beta_func_prime
  double precision, parameter  :: mu = 0.613 ! Offner et al 2009

  !Beta_core !a_rad, k_b, m_p are defined in constants.f90
  Beta_func_prime = 4. * (P_core**3 * a_rad / rho_core**4 / 3. * (mu * m_p / k_b)**4) * x**3 + 1.

  !Test function prime
  !func_prime = 2.*x
end function Beta_func_prime


function T_core_func(x, rho_core, P_core)
  use constants
  implicit none

  double precision             :: T_core_func
  double precision, intent(in) :: x, rho_core, P_core
  double precision, parameter  :: mu = 0.613 ! Offner et al 2009

  !T_core of polytrope   !a_rad, k_b, m_p are defined in constants.f90
  T_core_func = 1./3. * a_rad * x**4 + rho_core * k_b * x / (mu * m_p) - P_core
end function T_core_func

function T_core_func_prime(x, rho_core, P_core)
  use constants
  implicit none

  double precision, intent(in) :: x, rho_core, P_core
  double precision             :: T_core_func_prime
  double precision, parameter  :: mu = 0.613 ! Offner et al 2009

  !T_core of polytrope   !a_rad, k_b, m_p are defined in constants.f90
  T_core_func_prime = 4. * a_rad * x**3 / 3. + rho_core * k_b / (mu * m_p)
end function T_core_func_prime


