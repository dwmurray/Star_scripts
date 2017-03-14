module root_find
!------------------------------------------------------------------------------------------
! Find the root, provided the function and derivative.
!------------------------------------------------------------------------------------------
  use constants, only            : r_15
  use protostellar_core
  implicit none
  PUBLIC find_root

contains

  subroutine find_root(n, rho_core, P_core, xinit, machine_tolerance, maxiter, result)
    
    double precision, intent(in)  :: n, rho_core, P_core, xinit, machine_tolerance
    integer, intent(in)           :: maxiter
    double precision, intent(out) :: result
  
    double precision              :: func, func_prime
    double precision              :: x, xnew
    integer                       :: i

!    write(*,*) n, rho_core, P_core, xinit, machine_tolerance, maxiter
    
    result = 0.0
    x = xinit
    do i = 1, maxiter
!       write(*,*) i, x
!       xnew = x - (func(x) / func_prime(x))
       xnew = x - (func(x, rho_core, P_core) / func_prime(x, rho_core, P_core))
       
!       write(*,*) i, x
       if (abs(xnew-x) .le.  machine_tolerance) then
          result = xnew
          return
          !exit
       end if
       x = xnew
    end do
    !  if
    write(*,*) 'failure.'
  
  end subroutine find_root

end module root_find

function func(x, rho_core, P_core)
  use constants

  implicit none

  double precision             :: func
  double precision, intent(in) :: x, rho_core, P_core

  double precision, parameter  :: mu = 0.613 ! Offner et al 2009

!  write(*,*) a_rad, mu, m_p, k_b, rho_core, P_core
  !Beta_core !a_rad, k_b, m_p are defined in constants.f90
  func = (P_core**3 * a_rad / rho_core**4 / 3. * (mu * m_p / k_b)**4) * x**4 + x - 1.

  !T_core of polytrope   !a_rad, k_b, m_p are defined in constants.f90
  !func = 1./3. * a_rad * x**4 + rho_core * k_b * x / (mu * m_p) - P_core
  
  !Test function
  !func = x*x - 4.
  
end function func

function func_prime(x, rho_core, P_core)
  use constants

  implicit none

  double precision, intent(in) :: x, rho_core, P_core
  double precision             :: func_prime

  double precision, parameter  :: mu = 0.613 ! Offner et al 2009
!  write(*,*) a_rad, mu, m_p, k_b

  !Beta_core !a_rad, k_b, m_p are defined in constants.f90
  func_prime = 4. * (P_core**3 * a_rad / rho_core**4 / 3. * (mu * m_p / k_b)**4) * x**3 + 1.
  
  !T_core of polytrope   !a_rad, k_b, m_p are defined in constants.f90
  !func_prime = 4. * a_rad * x**3 / 3. + rho_core * k_b / (mu * m_p)

  !Test function prime
  !func_prime = 2.*x
  
end function func_prime


