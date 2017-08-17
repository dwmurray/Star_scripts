module protostellar_core
!------------------------------------------------------------------------------------------
! Functions that determine the properties of the core of the protostar for a given
! phi and epsilon from the Lane-Emden eqn.
!------------------------------------------------------------------------------------------

  use protostar_constants
  implicit none
  !PUBLIC rho_core, P_core

contains

  function determine_rho_core(M, R, eps1, dphideps1)

    real(dp)             :: determine_rho_core
    real(dp), intent(in) :: M, R, eps1, dphideps1
    real(dp)             :: rho_avg

    rho_avg = 3.*M / (4*pi*R**3)
    determine_rho_core = -(eps1 / 3.) * (1./dphideps1) * rho_avg

    !write(*,*) M, R, rho_avg
    !rho_core = M / (R**3)

  end function determine_rho_core

  function determine_P_core(M, R, n, dphideps1)

    real(dp), intent(in) :: M, R, n, dphideps1
    real(dp)             :: determine_P_core

    determine_P_core = G * M*M / ((n+1)*4*pi*R**4 * dphideps1**2)
    !P_core = 5. ! placeholder
    
  end function Determine_P_core

end module protostellar_core
