
program ps 
  use constants
  implicit none


  integer, parameter :: n_step = 10
  double precision, parameter :: logM_step = 0.1, logR_step = 0.1
  double precision, parameter :: n_min = 1.5, n_max = 3.5
  double precision, parameter :: M_min = 1.0D-2, M_max = 1.0D2 ! In Solar masses.
  double precision, parameter :: R_min = 1.0D-2, R_max = 1.0D3 ! In Solar radii

  ! n_step * (m_max - m_min)/logM_step * (R_max - R_min)/logR_step
  ! Is the minimum size for beta_table_dim.
!  integer, parameter :: beta_table_dim1 = 36000, beta_table_dim2 = 7 !dim2 is n, M, R, rho_core, P_core, Beta_core, Beta_mean
  integer, parameter :: beta_table_dim1 = n_step * (log10(M_max) - log10(M_min))/logM_step * (log10(R_max) - log10(R_min))/logR_step
  !The above technically drops the last line but said line is for n = 3.3 so we do not use it.
  integer, parameter :: beta_table_dim2 = 7!dim2 is n, M, R, rho_core, P_core, Beta_core, Beta_mean
  !this arrays will hold the table to interpolate beta from.
  double precision, dimension(beta_table_dim1, beta_table_dim2) :: Beta_file_table
  double precision            :: Beta, Beta_m_r, dlogBetadlogM
  double precision            :: Beta_core_m_r, dlogBeta_coredlogM

  ! These declarations are used for the driver to follow a protostar through its evolution.
  double precision           :: tage, mdot, deltaT, deltaM, lum, n, M, R, md
  double precision           :: n_in, M_in, R_in
  integer                    :: protostar_state

  !test
  double precision         :: mdot0

  write(*,*) beta_table_dim1
  call create_beta_table_file(n_min, n_max, n_step, M_min, M_max, logM_step, R_min, R_max, logR_step)
  !check reading the values back
  !fill the arrays with the values from the Beta table file.
  call init_beta_table_values(beta_table_dim1, beta_table_dim2, Beta_file_table)


  !as a test hand in a random n, M, R
  n_in = 1.4
  M_in = 1.001*M_sun
  R_in = 1.001*R_sun
  do while( n_in .LE. 3.0)
     ! Now interpolate and obtain Beta(n_in,M_in,R_in)
     call obtain_Beta(n_in, M_in, R_in, beta_table_dim1, beta_table_dim2, Beta_file_table, Beta_m_r, Beta_core_m_r)
     call obtain_dlogBetadlogM(n_in, M_in, R_in, Beta_m_r, Beta_core_m_r, &
          beta_table_dim1, beta_table_dim2, Beta_file_table, dlogBetadlogM, dlogBeta_coredlogM)

     write(*,*) Beta_m_r, Beta_core_m_r, dlogBetadlogM, dlogBeta_coredlogM

     n_in = n_in + 0.1
  end do
  write(*,*) 'Check of Beta & dlogBetadlogM complete.'
! Now proceed for the protostellar evolution.
  tage = 0.
  mdot0 =1e-6*M_sun/yr
  mdot = 1e-6*M_sun/yr
  deltaT = 1.e2*yr!0.01
  deltaM = mdot*deltaT

  m = 0.001 * M_sun
  protostar_state = 0
  do while (tage .le. 1.e6*yr)
     mdot = mdot0 + mdot0 * tage / (1D5*yr)
     deltaM = mdot*deltaT

!     if(tage > 1e6*yr) then
!     end if
     call update_protostar( deltaM, deltaT, m, r, n, md, protostar_state, &
          beta_table_dim1, beta_table_dim2, Beta_file_table)
     m = m + deltaM
     tage = tage + deltaT
  end do

  write(*,*) tage
  write(*,*) "m/m_sun = ", m/m_sun, "r/r_sun = ", r/r_sun
  write(*,*) r," ", lum
  
  stop
end program ps

subroutine update_protostar( deltaM, deltaT, M, R, n, md, protostar_state, &
     beta_table_dim1, beta_table_dim2, Beta_file_table)
  use constants
  implicit none

  double precision, intent( in)    :: deltaM, deltaT, M
  integer, intent( in)             :: beta_table_dim1, beta_table_dim2
  double precision, dimension(beta_table_dim1, beta_table_dim2), intent( in) :: Beta_file_table
  double precision, intent( inout) :: R, n, md
  integer, intent( inout)          :: protostar_state
  double precision, parameter      :: f_rad = 0.33 ! used to update to shell deuterium burning.
  double precision                 :: Beta_m_r, dlogBetadlogM, Beta_core_m_r, dlogBeta_core_dlogM
  double precision                 :: L_int, L_ion, L_D, deltaR
  double precision                 :: T_core
  !ZAMS checks
  double precision                 :: L_ZAMS, R_ZAMS, beta_ZAMS1

103 format(3(A5,(1pe14.5,2x)))
102 format(2(A6,(1pe14.5,2x)))
101 format(1(A5,(1pe14.5,2x)))

  if (protostar_state .EQ. 0) then
     ! if the mass is below 0.01 solar masses then we do not initialize
     if( m .GE. (0.01 * M_sun)) then! if it is now above 0.01 solar masses, it is in the no burning state.
        write(*,*) 'init protostar'
        protostar_state = 1
        call initialize_protostar( deltaM, deltaT, M, R, n, md, T_core) ! sets initial n, r, md, T_core for the particle.
     end if
     return
  end if

  if ( protostar_state .EQ. 1) then !No burning State.
     call obtain_Beta(n, M, R, beta_table_dim1, beta_table_dim2, Beta_file_table, Beta_m_r, Beta_core_m_r) !Returns Beta_m_r
     call obtain_dlogBetadlogM(n, M, R, Beta_m_r, Beta_core_m_r, &
          beta_table_dim1, beta_table_dim2, Beta_file_table, dlogBetadlogM, dlogBeta_core_dlogM) !Return dlogbetadlogM
     call set_L_internal(M, R, L_int) !Set L_int
     call set_L_ionize(deltaM, deltaT, L_ion) !Set L_ion
     call set_L_deuterium(n, M, R, deltaM, deltaT, protostar_state, L_int, L_ion, &
          Beta_m_r, Beta_core_m_r, dlogBetadlogM, dlogBeta_core_dlogM, L_D)! Set L_D
!     call set_L_deuterium(n, M, R, deltaM, deltaT, protostar_state, L_int, L_ion, Beta_m_r, dlogBetadlogM, L_D)! Set L_D
     write(*,102) 'M/Msun', M/M_sun, 'R/Rsun', R/R_sun
     write(*,103) 'L_int', L_int, 'L_ion', L_ion, 'L_D  ', L_D

     call update_radius( deltaM, deltaT, n, M, R, Beta_m_r, dlogBetadlogM, L_int, L_ion, L_D)!returns an updated R

     ! Find Tcore and check to see if we need to update the protostar state.
     call determine_T_core(n, M, R, T_core)
     write(*,101) 'T_c/T', T_core/1.5D6
!     write(*,101) 'MDeut', md
     if (T_core .GE. 1.5D6) then
        write(*,*) 'core burn'
        protostar_state = 2 ! update to core burning at fixed T_core
        !Offner et al. 2009 below B14
        n = 1.5
        !Init deuterium mass to burn.
        !Offner et al. 2009 B3
        md = m

     end if
     return
  end if
  
  if (protostar_state .EQ. 2) then ! Core burning at fixed T_core
     call obtain_Beta(n, M, R, beta_table_dim1, beta_table_dim2, Beta_file_table, Beta_m_r, Beta_core_m_r) !Returns Beta_m_r
     call obtain_dlogBetadlogM(n, M, R, Beta_m_r, Beta_core_m_r, &
          beta_table_dim1, beta_table_dim2, Beta_file_table, dlogBetadlogM, dlogBeta_core_dlogM) !Return dlogbetadlogM
     call set_L_internal(M, R, L_int) !Set L_int
     call set_L_ionize(deltaM, deltaT, L_ion) !Set L_ion
     call set_L_deuterium(n, M, R, deltaM, deltaT, protostar_state, L_int, L_ion, &
          Beta_m_r, Beta_core_m_r, dlogBetadlogM, dlogBeta_core_dlogM, L_D)! Set L_D
     call update_radius( deltaM, deltaT, n, M, R, Beta_m_r, dlogBetadlogM, L_int, L_ion, L_D)!returns an updated R
!     write(*,103) 'L_int', L_int, 'L_ion', L_ion, 'L_D', L_D
     ! Find Tcore
     call determine_T_core(n, M, R, T_core)
     write(*,101) 'T_fix', T_core/1.5D6

     call update_md(deltaM, deltaT, L_D, md) !returns an updated md
!     write(*,101) 'MDeut', md
     if (md .EQ. 0.0) then 
        write(*,*) 'core burn variable'
        protostar_state = 3 ! update to core burning at variable T_core
     end if
     return
  end if

  if (protostar_state .EQ. 3) then ! Core burning at Variable T_core
     call obtain_Beta(n, M, R, beta_table_dim1, beta_table_dim2, Beta_file_table, Beta_m_r, Beta_core_m_r) !Returns Beta_m_r
     call obtain_dlogBetadlogM(n, M, R, Beta_m_r, Beta_core_m_r, &
          beta_table_dim1, beta_table_dim2, Beta_file_table, dlogBetadlogM, dlogBeta_core_dlogM) !Return dlogbetadlogM
     call set_L_internal(M, R, L_int) !Set L_int
     call set_L_ionize(deltaM, deltaT, L_ion) !Set L_ion
     call set_L_deuterium(n, M, R, deltaM, deltaT, protostar_state, L_int, L_ion, &
          Beta_m_r, Beta_core_m_r, dlogBetadlogM, dlogBeta_core_dlogM, L_D)! Set L_D
     call update_radius( deltaM, deltaT, n, M, R, Beta_m_r, dlogBetadlogM, L_int, L_ion, L_D)!returns an updated R
     call determine_T_core(n, M, R, T_core)
!     write(*,101) 'T_var', T_core/1.5D6
!     write(*,101) 'MDeut', md
     !compare L_D to L_ms
     call set_ZAMS( m, L_ZAMS, R_ZAMS, beta_ZAMS1)
!     write
     write(*,103) 'L_D', L_D, 'L_ZAMS', L_ZAMS, 'ratio', L_D/L_ZAMS

     write(*,102) 'M/Msun', M/M_sun, 'R/Rsun', R/R_sun
!     write(*,103) 'L_int', L_int, 'L_ion', L_ion, 'L_D  ', L_D

!     write(*,103) 'L_int', L_int, 'L_ZAMS', L_ZAMS, 'ratio', L_int/L_ZAMS
     !Never seem to validate this condition.
     if( L_D / L_ZAMS .GE. f_rad) then !L_D/L_ZAMS > 0.33
        write(*,*) 'shell burn'
        protostar_state = 4 ! update to shell burning deuterium
        n = 3.
        R = 2.1 * R ! represents a swelling of star due to formation of the radiative barrier
        stop
     end if
     return
  end if

  if (protostar_state .EQ. 4) then
!     write(*,*) 'Shell Burning deuterium'
     call obtain_Beta(n, M, R, beta_table_dim1, beta_table_dim2, Beta_file_table, Beta_m_r, Beta_core_m_r) !Returns Beta_m_r
     call obtain_dlogBetadlogM(n, M, R, Beta_m_r, Beta_core_m_r, &
          beta_table_dim1, beta_table_dim2, Beta_file_table, dlogBetadlogM, dlogBeta_core_dlogM) !Return dlogbetadlogM
     call set_L_internal(M, R, L_int) !Set L_int
     call set_L_ionize(deltaM, deltaT, L_ion) !Set L_ion
     call set_L_deuterium(n, M, R, deltaM, deltaT, protostar_state, L_int, L_ion, &
          Beta_m_r, Beta_core_m_r, dlogBetadlogM, dlogBeta_core_dlogM, L_D)! Set L_D
     call update_radius( deltaM, deltaT, n, M, R, Beta_m_r, dlogBetadlogM, L_int, L_ion, L_D)!returns an updated R
!     write(*,103) 'L_internal', L_int, 'L_ionize', L_ion, 'L_D', L_D
     call determine_T_core(n, M, R, T_core)
!     write(*,101) 'T_c/T', T_core/1.5D6
!     write(*,101) 'MDeut', md


     !compare R to R_ZAMS
     call set_ZAMS( m + deltaM, L_ZAMS, R_ZAMS, beta_ZAMS1)
     write(*,103) 'R/R_sun', R/R_sun, 'RZAMS', R_ZAMS/R_sun, 'ratio', R/R_ZAMS

     if( R .LE. R_ZAMS) then
        protostar_state = 5 ! update to main sequence star.
        R = R_ZAMS !Set R to R_ZAMS and it stays at that value.
        write(*,*) 'At Main Sequence.'
        return
     end if
  end if

  if( protostar_state .EQ. 5) then
     call set_ZAMS( m + deltaM, L_ZAMS, R_ZAMS, beta_ZAMS1)
     R = R_ZAMS
  end if

end subroutine update_protostar


!subroutine protostar( tage, deltaM, deltaT, m, r, lum, n, md, T_core, protostar_state)!, phi_dim1, phi_dim2, phi_eps_table) 
!  use constants
!  implicit none
!  
!  double precision, intent( in) :: tage, deltaM, deltaT, m!, yr
!!  integer, intent(in) :: phi_dim1, phi_dim2
!!  double precision, intent( in), dimension(phi_dim1,phi_dim2) :: phi_eps_table
!
!  double precision, intent( out) :: lum
!  double precision, intent( inout) :: r, n, md, T_core
!  integer, intent( inout) :: protostar_state
!
!  double precision :: deltaR, L_D, L_MS
!
!  call update_radius( deltaM, deltaT, m, r, n, L_D, L_MS, deltaR)!( deltaM, deltaT, m, r)
!
!  return
!end subroutine protostar

subroutine initialize_protostar( deltaM, deltaT, m, r, n, md, T_core)
  use constants
  implicit none

  double precision, intent( in) :: deltaM, deltaT, m
  double precision, intent( out) :: r, n, md, T_core
  double precision :: mdot

  mdot = deltaM / deltaT
  write(*,*) mdot *yr/M_sun
  !init is from Offner et al. 2009 ApJ 703
  r = ( 2.5 * R_sun) * ( deltaM * yr / (deltaT * 1.E-5 * M_sun))**0.2

  n = 5. - 3. / (1.475 + 0.07 * log10(mdot*yr/M_sun))
  write(*,*) n
  if (n .LE. 1.5) then
     n = 1.5
  else if (n .GE. 3) then
     n = 3
  end if

  md = m
  T_core = 0.
  return
end subroutine initialize_protostar

subroutine update_radius( deltaM, deltaT, n, m, r, Beta, dlogBetadlogM, L_int, L_ion, L_D)
  use constants
  implicit none

  double precision, intent( in) :: deltaM, deltaT, n, m, Beta, dlogBetadlogM, L_int, L_ion, L_D
  double precision, intent( out) :: R !deltaR
  double precision :: deltaR
  double precision :: L_ZAMS, R_ZAMS, beta_ZAMS1 !Values from the ZAMS file. !beta_ZAMS1 to avoid conflict with constants beta_ZAMS
  double precision :: alpha_g! alpha_g describes the gravitational binding energy of a poltrope.

  double precision, parameter :: f_k = 0.5 !fraction of the kinetic energy of the infalling material that is radiated away.

  alpha_g = 3./(5. - n)

! Offner et al. deltaR
  deltaR = 2.* (deltaM/m) * (1 - (1-f_k)/(alpha_g*Beta) + &
       (1./2.)*(dlogBetadlogM)) * r - 2. * (deltaT / (alpha_g*Beta)) * (r/G/m/m)*(L_int + L_ion - L_D)*r

  R = R + deltaR
  ! Checking vs ZAMS
  call set_ZAMS( m + deltaM, L_ZAMS, R_ZAMS, beta_ZAMS1)

!103 format(4(A5,(1pe14.5,2x)))
!  write(*,103) 'DeltR', deltaR/R_sun, 'R    ', R/R_sun, 'RZAMS', R_ZAMS/R_sun, 'ratio', R/R_ZAMS

  if( R < R_ZAMS) then ! conditions to ZAMS
     R = R_ZAMS
  end if


  return
end subroutine update_radius


!subroutine update_radius( deltaM, deltaT, m, r)
!  ! Nakano et al. 95 ApJ 450...183N Appendix A equation 27
!  use constants
!
!  implicit none
!
!
!  double precision, intent( in) :: deltaM, deltaT, m
!  double precision, intent( inout) :: r
!
!  double precision :: L_I, L_D
!
!  double precision, parameter :: f_acc=1.0 ! This is a correction factor in Nakano 95. !TODO obtain this value
!  double precision, parameter :: a_e = 3./4. ! 3/4. is for n = 3
!!  double precision, parameter :: L_solar = 3.839D33
!!  double precision, parameter :: M_solar = 1.989D33
!!  double precision, parameter :: f_k = 0.5 !fraction of the kinetic energy of the infalling material that is radiated away.
!!  double precision, parameter :: Temp_H = 3000 ! surface temp of hayashi track star.
!  double precision, parameter :: sigma = 5.67D-5 !stefan-Boltzmann cgs.
!  double precision :: T_eff
!  double precision :: L_star, R_star ! ZAMS values
!  double precision :: mdot, beta
!  double precision :: deltaR
!!  write(*,*) G, pi, a_e, f_acc
!
!!  L_H = 4.*pi*r**2*sigma*Temp_H**4
!  
!! Rough estimate for L_star currently !TODO Update this estimate
!
!  call set_ZAMS( m + deltaM, L_star, R_star, beta)
!!  L_star = set_L_star(m, beta)!  L_star = L_sun * m*m*m / m_sun**3
!  !write(*,*) "called set_ZAMS", m+deltaM, L_star, R_star, beta  
!  mdot = deltaM/deltaT
!  L_D = 15. * L_sun * mdot/( 1.e-5 * M_sun/yr)
!  L_I = 2.5 * L_sun * mdot/( 1.e-5 * M_sun/yr) !dot{M} / 1e-5 M_sun per yr
!
!
!  !r = G*a_e*(m+deltaM)*mdot/L_star*(1. + beta - (1.+f_acc)/a_e*0.5)
!  deltaR = (2. - (1.+f_acc)/(2.*a_e)*r*deltaM/m - r/(a_e*G*m*m)*(L_star + L_I - L_D))*deltaT
!  r = r + deltaR 
!  !r = r + deltaR
!  !deltaR = (2. - (1. + f_acc) / (2.*a_e)) * r * deltaM / m - r*r*deltaT*L_star * deltaM / (a_e * G * deltaM*m*m) &
!  !     - r*r*deltaM / (a_e*G*m*m) * (L_I*deltaT / deltaM - L_D*deltaT / deltaM)
!  
!  if( r < R_star) then ! conditions to ZAMS
!     r = R_star
!  end if
!  write(*,*) "new mass radius", m+deltaM, r/r_sun, G*a_e*(m+deltaM)*mdot/L_star/r_sun, beta, R_star/r_sun
!  write(12,103) "L_* ",L_star , &
!       "L_I ", L_I , &
!       "L_D ", L_D , &
!       "f_acc ", (2. - (1. + f_acc) / (2.*a_e)), &
!       "RHS ", r*deltaT / (a_e * G * deltaM*m), &
!       "LHS ", r*deltaM /m
!103 format( 6(A6,(1pe10.4,2x)))
!
!!  deltaR = ((2. - (1. + f_acc) / (2.*a_e)) - r*deltaT*L_star / (a_e * G * deltaM*m) &
!!       - r / (a_e*G*m) * (L_I*deltaT / deltaM - L_D*deltaT / deltaM)) * (r*deltaM /m)
!
!  !deltaR =  - r*r*deltaT*L_star * deltaM / (a_e * G * deltaM*m*m) &
!  !deltaR = - r*r*deltaM / (a_e*G*m*m) * (L_I*deltaT / deltaM - L_D*deltaT / deltaM)
!  !T_eff = (L_star / (4.*pi*r*r*sigma_b))**0.25
!!  write(11,101) "T_eff", T_eff, "L_star", L_star/L_sun, "M_star", M/M_sun, "R_star", r/r_sun, "deltaR/r_sun", deltaR/r_sun
!!101 format( 4(A6,(1pe10.3,2x)), A12,(1pe10.3,2x))
!
!  return
!end subroutine update_radius

subroutine initialize_ZAMS
  use constants
  implicit none
  character(len=256) :: header
  double precision :: m, r, l, age
  integer :: line, ios
  open( unit=15, file="ZAMS_ezer67.txt", status="old")
  
  ! skip the first line
  read(15,*) header
  line = 0
  ios = 0
  ! read data
  do while( ios .eq. 0) 
     read(15, *, iostat=ios) m, r, l, age
     if( ios .eq. 0) then
        line = line + 1
        lgm_zams(line) = log10(m)
        lgl_zams(line) = log10(l)
        lgr_zams(line) = log10(r)
     end if
  end do
  close(unit=15)

  num_zams = line
  
  ! figure out beta
  do line = 1, num_zams - 1
     beta_zams(line) = (lgl_zams(line+1)-lgl_zams(line))/(lgm_zams(line+1)-lgm_zams(line)) 
     betaR_zams(line) = (lgr_zams(line+1)-lgr_zams(line))/(lgm_zams(line+1)-lgm_zams(line)) 
  end do
  beta_zams(num_zams) = beta_zams(num_zams-1)
  betaR_zams(num_zams) = betaR_zams(num_zams-1)
  initialized_zams = .true.

  return
end subroutine initialize_ZAMS

subroutine set_ZAMS( M_star, L_star, R_star, beta)
  use constants
  implicit none
  double precision, intent(in) :: M_star
  double precision, intent(out) :: L_star, R_star, beta
  double precision :: lgm, lgL, lgR
  integer :: i 
  if( .not. initialized_zams ) then
     call initialize_zams
  end if

  lgm = log10(M_star/m_sun)
  beta = 1.

  if( lgm < lgm_zams(1)) then
     lgL = lgl_zams(1) + beta_zams(1)*(lgm - lgm_zams(1))
     lgR = lgr_zams(1) + betaR_zams(1)*(lgm - lgm_zams(1))
     beta = beta_zams(1)
  elseif( lgm >= lgm_zams(num_zams)) then
     lgL = lgl_zams(num_zams) + beta_zams(num_zams)*(lgm - lgm_zams(num_zams))
     lgR = lgr_zams(num_zams) + betaR_zams(num_zams)*(lgm - lgm_zams(num_zams))
     beta = beta_zams(num_zams)
  else
     do i = 1, num_zams-1 
        if( lgm .ge. lgm_zams(i) .and. lgm < lgm_zams(i+1)) then
           lgL = lgl_zams(i) + beta_zams(i)*(lgm - lgm_zams(i))
           lgR = lgr_zams(i) + betaR_zams(i)*(lgm - lgm_zams(i))
           beta = beta_zams(i)
        end if
     end do
  end if
  
  L_star = 1d1**lgL*L_sun
  R_star = 1d1**lgR*R_sun

  return
end subroutine set_ZAMS

subroutine set_L_internal(M_star, R, L_int)
  use constants
  implicit none

  double precision, intent( in)  :: M_star, R
  double precision, intent( out) :: L_int
  double precision               :: L_Hayashi, L_ZAMS, R_ZAMS, beta
  double precision, parameter    :: T_Hayashi = 3000.

  L_Hayashi = 4.*pi*R**2*a_rad*T_Hayashi**4
  
  call set_ZAMS( M_star, L_ZAMS, R_ZAMS, beta)
  L_int = max(L_ZAMS,L_Hayashi)

end subroutine set_L_internal

subroutine set_L_ionize(deltaM, deltaT, L_ion)
  use constants
  implicit none

  double precision, intent( in)  :: deltaM, deltaT  
  double precision, intent( out) :: L_ion

  L_ion = 2.5 * L_sun * (deltaM/deltaT)/(1d-5 * M_sun/(pi*1D7))

end subroutine set_L_ionize

subroutine set_L_deuterium(n, M, R, deltaM, deltaT, protostar_state, L_int, L_ion, &
     Beta_m_r, Beta_core_m_r, dlogBeta_dlogM, dlogBeta_core_dlogM, L_D)! Set L_D
  use constants
  implicit none
  
  integer,          intent( in)   :: protostar_state
  double precision, intent( in)   :: n, M, R, deltaM, deltaT
  double precision, intent( in)   :: L_int, L_ion
  double precision, intent( in)   :: Beta_m_r, Beta_core_m_r, dlogBeta_dlogM, dlogBeta_core_dlogM
  double precision, intent( out)  :: L_D
  double precision, parameter     :: f_k = 0.5 !fraction of the kinetic energy of the infalling material that is radiated away.
  double precision                :: alpha_g! alpha_g describes the gravitational binding energy of a poltrope.
  double precision                :: mdot

!  !To Do! obtain dlogBeta_core / dlogM
!  dlogBeta_core_dlogM = -4.8911310542155833E-007 !until the subroutine call is put in.

  alpha_g = 3./(5. - n)
  mdot = deltaM / deltaT

  if( protostar_state .EQ. 1) then
     L_D = 0.0
  else if( protostar_state .EQ. 2) then
!     write(*,*) 'core burn'
     !Offner et al. 2009 expression B8
     L_D =  L_int + L_ion + G * M / R * mdot * (1. - f_k - alpha_g*Beta_m_r*0.5 * (1. + dlogBeta_dlogM - dlogBeta_core_dlogM))
  else
!     write(*,*) 'shell burn'
     !Offner et al. 2009 expression B9
     L_D = 15. *L_sun * mdot * yr / 1.e-5 / M_sun
!     L_D = 15. * L_sun * mdot/ 1.e-5 * M_sun/yr


  end if
  return
end subroutine set_L_deuterium

subroutine determine_T_core(n, M, R, T_core)
  use constants
  use rk2_polytrope
  use root_find
  implicit none

  double precision, intent( in)   :: n, M, R
  double precision, intent( out)  :: T_core
  double precision   :: rho_core, P_core, K_value, Beta_core
  double precision   :: Beta_mean
  double precision   :: logMmin, logMmax, logRmin, logRmax
  double precision   :: logM, logR
  double precision   :: T_core_by_betacore, T_core_by_Pcore
  double precision   :: eps1, dphideps1, phi1
  integer, parameter :: find_root_maxiter = 200, phi_eps_iter = 100

  double precision, parameter :: mu = 0.613 !Offner et al. '09 mean molecular weight for fully ionized gas of solar comp.
  double precision, parameter :: machine_tolerance = 1e-14, Beta_guess = 0.5
  logical            :: file_exist

  !To Do make a call to find root to explicitly solve for T_core.
  !We can not call to polytrope because we have CHANGED n across the different protostar states.

  call polytrope(n, eps1, dphideps1, phi1) !Returns eps1, dphideps1 and phi1 all eval. at eps1.
  !Find rho and P at the core for the given n, M, R
  rho_core = determine_rho_core(M, R, eps1, dphideps1)
  P_core = determine_P_core(M, R, n, dphideps1)
  
  !Solve for Beta_core
  call find_root(n, rho_core, P_core, Beta_guess, machine_tolerance, find_root_maxiter, Beta_core)



  ! Checks on the core temperature
  T_core_by_betacore = (k_b * 3. / (mu * m_p * a_rad) *(1 - Beta_core) / Beta_core)**(1./3.) * rho_core**(1./3.)
  T_core_by_Pcore = mu * m_p * Beta_core * P_core / k_b / rho_core
  ! These two methods of calculating the core temperature agree fairly well. 
  T_core = T_core_by_Pcore
!  T_core = (T_core_by_Pcore + T_core_by_betacore) * 0.5
end subroutine determine_T_core

subroutine update_md(deltaM, deltaT, L_D, md)
  use constants
  implicit none
  double precision, intent( in) :: deltaM, deltaT, L_D
  double precision, intent( inout) :: md
  double precision :: deltaMd

!  deltaMd = deltaM - 1.D-5 * M_solar * (L_D / (15.*L_solar)) * deltaT
  deltaMd = deltaM - 1.D-5 * M_sun * L_D / (15. * L_sun) * deltaT / yr

!  write(*,*) 'deltaMd', deltaMd, 'deltaM', deltaM, 'deltaT', deltaT, 'L_D', L_D
!  write(*,*) 1.D-5 * M_sun * L_D / (15. * L_sun) * deltaT / yr

!103 format(2(A5,(1pe14.5,2x)))
!  write(*,103) 'delmd', deltaMd, 'md   ', md


  md = md + deltaMd !deltaMd is negative.

  if (md .LT. 0.0) then
     md = 0.0
  end if
  return
end subroutine update_md


subroutine create_beta_table_file(n_min, n_max, n_step, M_min, M_max, logMstep, R_min, R_max, logRstep)
  use constants
  use rk2_polytrope
  use root_find

  implicit none
  integer, intent( in) :: n_step
  double precision, intent( in) :: n_min, n_max, M_min, M_max, logMstep, R_min, R_max, logRstep
  double precision   :: n, M, R
  double precision   :: eps1, dphideps1, phi1
  double precision   :: rho_core, P_core, K_value, Beta_core
  double precision   :: Beta_mean
  double precision   :: logMmin, logMmax, logRmin, logRmax
  double precision   :: logM, logR
  double precision   :: T_core_by_betacore, T_core_by_Pcore, R_star_check
  integer            :: i
  integer, parameter :: find_root_maxiter = 200, phi_eps_iter = 100

  double precision, parameter :: mu = 0.613 !Offner et al. '09 mean molecular weight for fully ionized gas of solar comp.
  double precision, parameter :: machine_tolerance = 1e-14, Beta_guess = 0.5
  logical            :: file_exist

  !Set initial n, M, R to loop through.
  n = n_min

  logMmin = log10(M_min)
  logMmax = log10(M_max)
!  logMstep = (logMmax - logMmin) / M_step
!  logMstep = 0.1

  logRmin = log10(R_min)
  logRmax = log10(R_max)
!  logRstep = (logRmax - logRmin) / R_step ! -> set to 0.1 for R and M
!  logRstep = 0.1
!  n = n_min + (n_max-n_min)*i/n_step
  write(*,*) n_min, n_max, M_min, M_max, logMstep, R_min, R_max, logRstep

  do i = 1, n_step !stepping in n
     !Solve Lane-Emden equation for the current n
     call polytrope(n, eps1, dphideps1, phi1) !Returns eps1, dphideps1 and phi1 all eval. at eps1.
     logM = logMmin
     do while (logM .LE. logMmax)!Stepping in Mass
        M = 10**logM * M_sun
        logR = logRmin
        do while (logR .LE. logRmax)!Stepping in Radius.
           R = 10**logR * R_sun
           !Find rho and P at the core for the given n, M, R
           rho_core = determine_rho_core(M, R, eps1, dphideps1)
           P_core = determine_P_core(M, R, n, dphideps1)

           !Solve for Beta_core
           call find_root(n, rho_core, P_core, Beta_guess, machine_tolerance, find_root_maxiter, Beta_core)

           !Obtain value of K for the polytrope !P = K rho**gamma
           K_value = P_core / rho_core**(1+1/n)

!           ! Checks on the core temperature
!           T_core_by_betacore = (k_b * 3. / (mu * m_p * a_rad) *(1 - Beta_core) / Beta_core)**(1./3.) * rho_core**(1./3.)
!           T_core_by_Pcore = mu * m_p * Beta_core * P_core / k_b / rho_core
!           write(*,*)'Temp(1D6)', T_core_by_betacore/1.e6, T_core_by_Pcore/1.e6, 'ratio', T_core_by_Pcore/T_core_by_betacore
!
!           ! Check the stellar radius
!           R_star_check = ((n+1)*K_value/4./pi/G)**0.5 *rho_core**((1-n)/2./n) * eps1
!           write(*,*) 'R/R_sun', 'Input: ', R/R_sun , 'Check: ', R_star_check/R_sun, 'ratio: ', R_star_check/R

           call determine_Beta_mean(n, eps1, M, R, rho_core, P_core, K_value, Beta_core, find_root_maxiter, phi_eps_iter, Beta_mean)
!           write(*,*) n, M/M_sun, R/R_sun
!           write(*,*) 'Beta_core', Beta_core, 'Beta_mean', Beta_mean

           inquire(file="Beta_interpolation_file.txt", exist=file_exist)
           if (file_exist) then
              open(12, file="Beta_interpolation_file.txt", status="old", position="append", action="write")
           else
              open(12, file="Beta_interpolation_file.txt", status="new", action="write")
           end if
           write(12, 101) n, m, R, rho_core, P_core, Beta_core, Beta_mean
101           format(8(1pe14.7,2x))
           close(12)

           !Update radius to next radius
           logR = logR + logRstep
        end do
        !Having stepped through all radii at that mass, update to next mass.
        logM = logM + logMstep
     end do
     !Having stepped through all mass and radii at this index, update to next index.
     n = n_min + (n_max-n_min)*i/n_step
  end do
end subroutine create_beta_table_file


!subroutine init_phi_eps(phi_dim1, phi_dim2, phi_eps_table)
!  use rk2_polytrope
!  implicit none
!  integer, intent(in) :: phi_dim1, phi_dim2
!  double precision, intent( inout), dimension(phi_dim1,phi_dim2) :: phi_eps_table
!!  double precision, dimension(lane_emden_dim2, lane_emden_dim1, lane_emden_dim1) :: test_table
!  double precision :: n, eps1, dphideps1, phi1
!  integer :: i
!  double precision, parameter :: nstart=1.4, nend=3.1
!
!!  n = 1.9
!  do i =1 , phi_dim2
!
!     n = nstart + (nend-nstart)*i/phi_dim2 !This is 1001 steps across 3.1-1.4! (3.1-1.4)/1000.!0.01
!!     write(*,*) n
!     call polytrope(n, eps1, dphideps1, phi1)
!!     write (*,*) phi1!n, eps1, dphideps1, phi1, -eps1**2*dphideps1
!     phi_eps_table(1,i) = n
!     phi_eps_table(2,i) = eps1
!     phi_eps_table(3,i) = dphideps1
!     phi_eps_table(4,i) = phi1
!     phi_eps_table(5,i) = -eps1**2*dphideps1
!!     write (*,*) n, eps1, phi_eps_table(5,i), phi1
!  end do
!  write(*,*) 'populated phi_eps_table'
!end subroutine init_phi_eps
!
!subroutine init_P_and_rho_core(phi_dim1, phi_dim2, phi_eps_table, proto_dim1, proto_dim2, proto_core_table)
!  use protostellar_core
!  use constants
!  use root_find
!
!  implicit none
!  
!  integer, intent(in) :: phi_dim1, phi_dim2
!  double precision, intent( in), dimension(phi_dim1,phi_dim2)          :: phi_eps_table
!
!  integer, intent(in) :: proto_dim1, proto_dim2
!  double precision, intent( out), dimension(proto_dim1, proto_dim2)    :: proto_core_table
!
!  double precision    :: rho_core, P_core, K_value, Beta_core
!  double precision    :: T_core_by_betacore, T_core_by_Pcore !check the temperature two ways.
!  double precision, parameter :: mu = 0.613 !Offner et al. '09 mean molecular weight for fully ionized gas of solar comp.
!  double precision    :: R_star_check
!  double precision    :: M_init, R_init, M, R, logR, logM
!  double precision, parameter :: dlogR=0.1, dlogM=0.1
!  double precision    :: n, eps, phi, dphideps, write_reset
!
!  integer             :: i,j,k,table_count
!
!  double precision, parameter    :: machine_tolerance = 1e-14, Beta_guess = 0.5
!  integer, parameter    :: maxiter = 20
!
!  ! for a solar mass then scale linearly.
!  !R = ( 2.5 * R_sun) !* ( deltaM / (deltaT * 1.E-5 * M_solar))**0.2
!  table_count = 1
!  i = 1
!  do i = 1, phi_dim2
!     n = phi_eps_table(1,i)
!     eps = phi_eps_table(2,i)
!     dphideps = phi_eps_table(3,i)
!     phi = phi_eps_table(4,i)
!!     write(*,*) n, eps, dphideps, phi
!     M_init = 1.0*M_sun
!     logM = log10(M_init)
!     do j = 1, 2!300!600
!        write_reset = 0.
!        M = 10**logM
!        R_init = 1.0*R_sun
!        logR = log10(R_init)
!!        R = 0.1*R_sun
!        do k = 1, 1!100
!           R = 10**logR
!
!           !Find rho and P at the core for the given n, M, R
!           rho_core = determine_rho_core(M, R, eps, dphideps)
!           P_core = determine_P_core(M, R, n, dphideps)
!
!           !Solve for Beta_core
!           call find_root(n, rho_core, P_core, Beta_guess, machine_tolerance, maxiter, Beta_core)
!
!           !Obtain value of K for the polytrope !P = K rho**gamma
!           K_value = P_core / rho_core**(1+1/n)
!
!!           ! Checks on the core temperature
!!           T_core_by_betacore = (k_b * 3. / (mu * m_p * a_rad) *(1 - Beta_core) / Beta_core)**(1./3.) * rho_core**(1./3.)
!!           T_core_by_Pcore = mu * m_p * Beta_core * P_core / k_b / rho_core
!!           write(*,*) 'mass', M/M_sun
!!           write(*,*)'Temp', T_core_by_betacore/1.e6, T_core_by_Pcore/1.e6
!!
!!           ! Check the stellar radius
!!           R_star_check = ((n+1)*K_value/4./pi/G)**0.5 *rho_core**((1-n)/2./n) * eps
!!           write(*,*) 'Input: ', R/7.D10 , 'Check: ', R_star_check/7.D10
!!
!           if (write_reset .EQ. 0.) then
!              write(*,102) n, M, R
!102           format(3(1pe14.5,2x))
!              write_reset = 1.
!           end if
!
!           !This if statement is while I do not have log steps, I don't want to seg fault by attempting to write beyond the table.
!           if (table_count .LE. proto_dim2) then
!
!              proto_core_table(1,table_count) = n
!              proto_core_table(2,table_count) = M
!!              proto_core_table(3,table_count) = logM
!              proto_core_table(3,table_count) = R
!!              proto_core_table(5,table_count) = logR
!              proto_core_table(4,table_count) = rho_core
!              proto_core_table(5,table_count) = P_core
!              proto_core_table(6,table_count) = K_value
!              proto_core_table(7,table_count) = Beta_core
!!              write(*,*) proto_core_table(:,table_count)
!           else
!              write(*,*) 'Looping beyond size of table for core values. Missing this data!'
!           end if
!
!           table_count = table_count + 1
!!              !write out to a csv and then plot with python to make sure there are no bumps etc.
!           logR = logR + dlogR
!!           R = R + 0.1*R_sun
!        end do
!        logM = logM + dlogM
!!       M = M + 1.0*M_sun
!     end do
!!     write(*,*) proto_core_table(2,:)
!  end do
!
!end subroutine init_P_and_rho_core


subroutine determine_Beta_mean(n, eps1, M_star, R_star, rho_core, P_core, K_value, Beta_core, &
     find_root_maxiter, phi_eps_iter, Beta_mean)
  ! This subroutine will calculate phi(eps), from that rho(eps) and then Beta(eps) and Beta_mean w.r.t. Pressure.
  ! This subroutine will calculate rho(r) = rho_core * phi(r)**n
  ! and P(r) = K_value * rho(r)**gamma
  ! and then call root_find.f90 to find the beta(r)
  ! which is then used to integrate over to find the mean beta.
  use rk2_polytrope
  use root_find
  use constants
  implicit none

  double precision, intent( in)  :: n, eps1, M_star, R_star, rho_core, P_core, K_value, Beta_core
  integer, intent( in)           :: find_root_maxiter, phi_eps_iter
  double precision, intent( out) :: Beta_mean
  double precision, parameter    :: machine_tolerance = 1e-14, Beta_guess = 0.5
  double precision :: eps_r, dphideps_r, phi_r
  double precision :: eps, deps, phi, dphideps, r_0
  double precision :: rho_r, P_r, Beta_r, r, dr
  double precision :: Pressure_Beta_Integral, Pressure_Integral
  integer          :: i
  !Initialize Integrals at zero
  Pressure_Beta_Integral = 0
  Pressure_Integral = 0

  ! set initial conditions
  eps = 1.0D-5
  deps = (eps1 - eps) / phi_eps_iter
  !Some initial guess for phi(eps) and dphi(eps)deps
  phi = 1. - eps*eps/6. + n/120.*eps**4.
  dphideps = -eps/3. + n/30.*eps**3.
  
  r_0 = ((n+1)*K_value/4./pi/G)**0.5 *rho_core**((1-n)/2./n) !scale for radius
  !set dr for the integral
  dr = r_0 * deps
  
  do i = 1, phi_eps_iter
     call phi_of_r( n, eps, deps, phi, dphideps, eps_r, dphideps_r, phi_r) !rk2 solve from eps=0 to eps
     
     if (phi_r .LT. 0.0) then
        write(*,*) 'phi less than zero'
        exit
     end if
     eps = eps_r
     phi = phi_r
     dphideps = dphideps_r
     
     r = r_0 * eps
     rho_r = rho_core * phi**n
     P_r = K_value * rho_r
     call find_root(n, rho_r, P_r, Beta_guess, machine_tolerance, find_root_maxiter, Beta_r)
     Pressure_Beta_Integral = Pressure_Beta_Integral + 4.*pi*P_r*Beta_r*r**2*dr
     Pressure_Integral = Pressure_Integral + 4.*pi*P_r*r**2*dr
!     write(*,*) 'Beta_r', Beta_r, 'rho_r   ', rho_r, 'P_r   ', P_r, 'r', r/7D10
  end do

  Beta_mean = Pressure_Beta_Integral/Pressure_Integral

!  write(*,*) Pressure_Beta_Integral, Pressure_Integral, Beta_mean
!  write(*,*)'Beta_core', Beta_core, 'Beta_mean', Beta_mean

end subroutine determine_Beta_mean


!subroutine obtain_interp_values_from_beta_table(n_in, M_in, R_in, beta_table_dim)
!
!  use constants
!  implicit none
!  double precision, intent( in)               :: n_in, M_in, R_in, beta_table_dim
!!  double precision, dimension(beta_table_dim) :: n_table, M_table, R_table, rho_c_table, P_c_table, Beta_c_table, Beta_mean_table
!  double precision                            :: n, M, R, rho_core, P_core, Beta_core, Beta_mean
!  double precision                            :: logM_in, logR_in, logBeta_m_r
!  double precision                            :: n_old, M_old, R_old, rho_core_old, P_core_old, Beta_core_old, Beta_mean_old
!  integer                                     :: line, ios
!  integer, parameter                          :: nvars = 2 
!  double precision, dimension(nvars, nvars)   :: n_bracket, M_bracket, R_bracket
!  double precision, dimension(nvars, nvars)   :: rho_core_bracket, P_core_bracket, Beta_core_bracket
!  double precision, dimension(nvars, nvars)   :: Beta_mean_bracket
!  
!  logM_in = log10(M_in)
!  logR_in = log10(R_in)
!
!  !Now obtain the values for X(mi+,ri) and X(mi+1, ri+1)
!  open( unit=16, file="Beta_interpolation_file.txt", status="old")
!  line = 0
!  ios = 0
!  ! read data
!  do while( ios .eq. 0) 
!     read(16, *, iostat=ios) n, M, R, rho_core, P_core, Beta_core, Beta_mean
!     if( ios .eq. 0) then
!        line = line + 1
!        if(line == 1) then
!           write(*,*) n, M, R
!        end if
!        !We are always overestimating in n. the way we get around this is by having many n's to compare to.
!        !thus, if our n_in = 1.6 and our n = 1.61 we will use the correct M_old<M_in<M and R_old<R_in<R
!        ! for this n > n_in.
!        if(n .GE. n_in) then
!           if(M .GE. M_in) then
!              if(R .LT. R_in) then
!                 !We now have the m,r,rho's etc that create the lower right point to bracket the input m_in and r_in
!                 ! on a plot of m vs r.
!                 !convert to log as we will use a linear interpolation in logspace.
!                 !  logM_in = log10(M_in)
!                 n_bracket(2,1) = n
!                 M_bracket(2,1) = log10(M)
!                 R_bracket(2,1) = log10(R)
!                 rho_core_bracket(2,1) = log10(rho_core)
!                 P_core_bracket(2,1) = log10(P_core)
!                 Beta_core_bracket(2,1) = log10(Beta_core)
!                 Beta_mean_bracket(2,1) = log10(Beta_mean)
!              else if(R .GE. R_in) then
!                 n_bracket(2,2) = n
!                 M_bracket(2,2) = log10(M)
!                 R_bracket(2,2) = log10(R)
!                 rho_core_bracket(2,2) = log10(rho_core)
!                 P_core_bracket(2,2) = log10(P_core)
!                 Beta_core_bracket(2,2) = log10(Beta_core)
!                 Beta_mean_bracket(2,2) = log10(Beta_mean)
!                 exit
!              end if
!           end if
!        end if
!     end if !ends if ios
!  end do
!!  close(unit=16)
!  !Now obtain the values for X(mi,ri) and X(mi, ri+1)
!!  open( unit=16, file="Beta_interpolation_file.txt", status="old")
!  rewind(unit=16)
!  line = 0
!  ios = 0
!  ! read data
!  do while( ios .eq. 0) 
!     read(16, *, iostat=ios) n, M, R, rho_core, P_core, Beta_core, Beta_mean
!     if( ios .eq. 0) then
!        line = line + 1
!        if(line == 1) then
!           write(*,*) n, M, R
!        end if
!        !We are always overestimating in n. the way we get around this is by having many n's to compare to.
!        !thus, if our n_in = 1.6 and our n = 1.61 we will use the correct M_old<M_in<M and R_old<R_in<R
!        ! for this n > n_in.
!        if(n .GE. n_in) then
!           if(M .LT. M_in) then! .and. M .GT. M_in) then
!              if(R .LT. R_in) then
!                 n_bracket(1,1) = n
!                 M_bracket(1,1) = log10(M)
!                 R_bracket(1,1) = log10(R)
!                 rho_core_bracket(1,1) = log10(rho_core)
!                 P_core_bracket(1,1) = log10(P_core)
!                 Beta_core_bracket(1,1) = log10(Beta_core)
!                 Beta_mean_bracket(1,1) = log10(Beta_mean)
!              else if(R .GE. R_in) then
!                 n_bracket(1,2) = n
!                 M_bracket(1,2) = log10(M)
!                 R_bracket(1,2) = log10(R)
!                 rho_core_bracket(1,2) = log10(rho_core)
!                 P_core_bracket(1,2) = log10(P_core)
!                 Beta_core_bracket(1,2) = log10(Beta_core)
!                 Beta_mean_bracket(1,2) = log10(Beta_mean)
!                 exit
!              end if
!           end if
!        end if
!     end if !ends if ios
!     n_old = n
!     M_old = M
!     R_old = R
!     rho_core_old = rho_core
!     P_core_old = P_core
!     Beta_core_old = Beta_core
!     Beta_mean_old = Beta_mean
!  end do
!  close(unit=16)
!
!  !Check that we are correctly bracketing what we want.
!103 format(3(1pe14.7,2x))
!104 format(5(1pe14.7,2x))
!108 format(8(1pe14.7,2x))
!
!!  write(*,104) n_mi_ri, n_mi_ri_1, n_mi_1_ri, n_mi_1_ri_1, n_in
!  write(*,104) n_bracket(1,1),n_bracket(1,2), n_in, n_bracket(2,1),n_bracket(2,2)
!  write(*,104) 10**M_bracket(1,1), 10**M_bracket(1,2), M_in, 10**M_bracket(2,1), 10**M_bracket(2,2)
!  write(*,104) 10**R_bracket(1,1), 10**R_bracket(1,2), R_in, 10**R_bracket(2,1), 10**R_bracket(2,2)
!
!
!! and now interpolate (in logspace) to determine what the rho_core, P_core, Beta_core and Beta_mean are
!  ! interpolate over Beta_mean. and return the linear interpolated value of Beta_mean(m_in, r_in)
!  call interpolate(nvars, logM_in, logR_in, M_bracket, R_bracket, Beta_mean_bracket, logBeta_m_r)!feed in X_bracket, M_bracket, R_bracket, nvars
!!  logrho_core_out = logrho_core_old + dlogrho/dlogm_in * logdeltaM
!! logdeltaM = logM_in - logM_old
!!  logrho = logrho_i + (logrho_i+1 - logrho_i)/(logm_i+1 - logm_i) *logdeltaM
! !inter across m with fixed R
!! logrho_core_out = logrho_core_old + (logrho_core - logrho_core_old)/(logM - logM_old) * logdeltaM
! ! now interp across R with fixed M
!
!
!! rho_core_oldi = 10**logrho_core_old
!! rho_core_out = 10**logrho_core_out
!! rho_corei = 10**logrho_core
!! write(*,*) logrho_core_old, rho_core_oldi, rho_core_old
!! write(*,*) logrho_core_out, rho_core_out
!! write(*,*) logrho_core, rho_corei, rho_core
!! 
!
!end subroutine obtain_interp_values_from_beta_table
!
!
subroutine init_beta_table_values(beta_table_dim1, beta_table_dim2, Beta_file_table)
  use constants
  implicit none

  integer, intent( in)                        :: beta_table_dim1, beta_table_dim2
  double precision, dimension(beta_table_dim1, beta_table_dim2) :: Beta_file_table
  double precision                            :: n, M, R, rho_core, P_core, Beta_core, Beta_mean
  integer                                     :: line, ios

  open( unit=16, file="Beta_interpolation_file.txt", status="old")
  line = 0
  ios = 0
  ! read data
  do while( ios .eq. 0) 
     read(16, *, iostat=ios) n, M, R, rho_core, P_core, Beta_core, Beta_mean
     if( ios .eq. 0) then
        line = line + 1
        Beta_file_table(line, 1) = n
        Beta_file_table(line, 2) = M
        Beta_file_table(line, 3) = R
        Beta_file_table(line, 4) = rho_core
        Beta_file_table(line, 5) = P_core
        Beta_file_table(line, 6) = Beta_core
        Beta_file_table(line, 7) = Beta_mean
     end if
  end do
  close(unit=16)
end subroutine init_beta_table_values


subroutine obtain_Beta(n_in, M_in, R_in, beta_table_dim1, beta_table_dim2, Beta_file_table, Beta_m_r, Beta_core_m_r)
  use constants
  implicit none

  double precision, intent( in)                :: n_in, M_in, R_in
  integer, intent( in)                         :: beta_table_dim1, beta_table_dim2
  double precision, intent( in), dimension(beta_table_dim1, beta_table_dim2) :: Beta_file_table
  !Populate these from Beta_file_table.
  double precision, dimension(beta_table_dim1) :: n_table, M_table, R_table, rho_c_table, P_c_table
  double precision, dimension(beta_table_dim1) :: Beta_c_table, Beta_mean_table

  double precision, intent( out)               :: Beta_m_r, Beta_core_m_r
  double precision                             :: n, M, R, rho_core, P_core, Beta_core, Beta_mean
  double precision                             :: logM_in, logR_in, logBeta_m_r, logBeta_core_m_r
  integer                                      :: line, index
  integer, parameter                           :: nvars = 2
  double precision, dimension(nvars, nvars)    :: n_bracket, M_bracket, R_bracket
  double precision, dimension(nvars, nvars)    :: rho_core_bracket, P_core_bracket, Beta_core_bracket
  double precision, dimension(nvars, nvars)    :: Beta_mean_bracket
  
  logM_in = log10(M_in)
  logR_in = log10(R_in)
  line = 0
  do while(line .LE. beta_table_dim1)
     line = line + 1
     n_table(line)         = Beta_file_table(line, 1)
     M_table(line)         = Beta_file_table(line, 2)
     R_table(line)         = Beta_file_table(line, 3)
     rho_c_table(line)     = Beta_file_table(line, 4)
     P_c_table(line)       = Beta_file_table(line, 5)
     Beta_c_table(line)    = Beta_file_table(line, 6)
     Beta_mean_table(line) = Beta_file_table(line, 7)
  end do

  !We are always overestimating in n. the way we get around this is by having many n's to compare to.
  !thus, if our n_in = 1.6 and our n = 1.61 we will use the correct M_old<M_in<M and R_old<R_in<R
  ! for this n (1.61) > n_in (1.6).
  index = 1
  do while (index .LE. beta_table_dim1)
     if(n_table(index) .GE. n_in) then
        if(M_table(index) .LT. M_in .and. M_table(index+1) .GE. M_in) then
        end if
        if(R_table(index) .LT. R_in .and. R_table(index+1) .GE. R_in) then
           if(M_table(index) .LT. M_in) then
              n_bracket(1,1) = n_table(index)
              M_bracket(1,1) = log10(M_table(index))
              R_bracket(1,1) = log10(R_table(index))
              rho_core_bracket(1,1) = log10(rho_c_table(index))
              P_core_bracket(1,1) = log10(P_c_table(index))
              Beta_core_bracket(1,1) = log10(Beta_c_table(index))
              Beta_mean_bracket(1,1) = log10(Beta_mean_table(index))

              n_bracket(1,2) = n_table(index+1)
              M_bracket(1,2) = log10(M_table(index+1))
              R_bracket(1,2) = log10(R_table(index+1))
              rho_core_bracket(1,2) = log10(rho_c_table(index+1))
              P_core_bracket(1,2) = log10(P_c_table(index+1))
              Beta_core_bracket(1,2) = log10(Beta_c_table(index+1))
              Beta_mean_bracket(1,2) = log10(Beta_mean_table(index+1))

           else if(M_table(index) .GE. M_in) then
              n_bracket(2,1) = n_table(index)
              M_bracket(2,1) = log10(M_table(index))
              R_bracket(2,1) = log10(R_table(index))
              rho_core_bracket(2,1) = log10(rho_c_table(index))
              P_core_bracket(2,1) = log10(P_c_table(index))
              Beta_core_bracket(2,1) = log10(Beta_c_table(index))
              Beta_mean_bracket(2,1) = log10(Beta_mean_table(index))

              n_bracket(2,2) = n_table(index+1)
              M_bracket(2,2) = log10(M_table(index+1))
              R_bracket(2,2) = log10(R_table(index+1))
              rho_core_bracket(2,2) = log10(rho_c_table(index+1))
              P_core_bracket(2,2) = log10(P_c_table(index+1))
              Beta_core_bracket(2,2) = log10(Beta_c_table(index+1))
              Beta_mean_bracket(2,2) = log10(Beta_mean_table(index+1))
              exit
           end if
        end if
     end if
     index = index + 1
  end do

104 format(5(1pe14.7,2x))
!  write(*,104) n_bracket(1,1),n_bracket(1,2), n_in, n_bracket(2,1),n_bracket(2,2)
!  write(*,104) 10**M_bracket(1,1), 10**M_bracket(1,2), M_in, 10**M_bracket(2,1), 10**M_bracket(2,2)
!  write(*,104) 10**R_bracket(1,1), 10**R_bracket(1,2), R_in, 10**R_bracket(2,1), 10**R_bracket(2,2)

  call interpolate(nvars, logM_in, logR_in, M_bracket, R_bracket, Beta_mean_bracket, logBeta_m_r)!feed in X_bracket, M_bracket, R_bracket, nvars
!  write(*,104) 10**Beta_mean_bracket(1,1), 10**Beta_mean_bracket(1,2), &
!       10**logBeta_m_r, 10**Beta_mean_bracket(2,1), 10**Beta_mean_bracket(2,2)

  Beta_m_r = 10**logBeta_m_r
!  write(*,*) logBeta_m_r, Beta_m_r

  call interpolate(nvars, logM_in, logR_in, M_bracket, R_bracket, Beta_core_bracket, logBeta_core_m_r)!feed in X_bracket, M_bracket, R_bracket, nvars
  Beta_core_m_r = 10**logBeta_core_m_r

end subroutine obtain_Beta

subroutine obtain_dlogBetadlogM(n_in, M_in, R_in, Beta_m_r, Beta_core_m_r, &
     beta_table_dim1, beta_table_dim2, Beta_file_table, dlogBetadlogM, dlogBeta_coredlogM)
  use constants
  implicit none

  double precision, intent( in)               :: n_in, M_in, R_in, Beta_m_r, Beta_core_m_r
  integer, intent( in)                        :: beta_table_dim1, beta_table_dim2
  double precision, intent( in), dimension(Beta_table_dim1, Beta_table_dim2) :: Beta_file_table
  double precision, intent( out)              :: dlogBetadlogM, dlogBeta_coredlogM
  !for obtaining dlogBeta/dlogM
  double precision, parameter                 :: deriv_step = 1.D-1
  double precision                            :: n_prime, M_prime, R_prime, Beta_prime_m_r, Beta_prime_core_m_r
  double precision                            :: dlogBeta, dlogBeta_core, dlogM

  !We now have Beta(M_in,R_in), so we now pick a location close to Beta(M_in+M_in*deriv_step,R_in+R_in*deriv_step)
  !Interp to find that beta value, and then use that to determine dlogBeta/dlogM
  n_prime = n_in
  M_prime = M_in + M_in * deriv_step
  R_prime = R_in

  call obtain_Beta(n_prime, M_prime, R_prime, Beta_table_dim1, Beta_table_dim2, Beta_file_table, &
       Beta_prime_m_r, Beta_prime_core_m_r)
!  write(*,*) 'Beta & Beta prime', Beta_m_r, Beta_prime_m_r

  dlogBeta = log10(Beta_prime_m_r) - log10(Beta_m_r)
  dlogM    = log10(M_prime) - log10(M_in)
  dlogBetadlogM = dlogBeta/dlogM


  dlogBeta_core = log10(Beta_prime_core_m_r) - log10(Beta_core_m_r)
  dlogM    = log10(M_prime) - log10(M_in)
  dlogBeta_coredlogM = dlogBeta_core/dlogM

end subroutine obtain_dlogBetadlogM

subroutine interpolate(nvars, logM_in, logR_in, M_bracket, R_bracket, X_bracket,  logX_m_r)
  implicit none
  ! This is a linear interpolation in log space. Everything handed in and out should be in log
  integer, intent( in)            :: nvars
  double precision, intent( in), dimension(nvars, nvars) :: M_bracket, R_bracket, X_bracket!these are in log space
  double precision, intent( in)   :: logM_in, logR_in
  double precision, intent( out) :: logX_m_r
  double precision               :: delta_logm, delta_logr
  double precision               :: logX_m_r1, logX_m_r2

  delta_logm = logM_in - M_bracket(1,1)
  delta_logr = logR_in - R_bracket(1,1)

  !mi,ri     = 1,1
  !mi,ri+1   = 1,2
  !mi+1,ri   = 2,1
  !mi+1,ri+1 = 2,2

  !Interpolating over M, then interpolating over R.
  !logX(m,r1) = X(m1,r1) + ((X(m2,r1) - X(m1,r1)) / (M(m2) - M(m1))) - M(m1,r1)
  logX_m_r1 = X_bracket(1,1) + (X_bracket(2,1) - X_bracket(1,1))/ (M_bracket(2,1) - M_bracket(1,1)) * (logM_in - M_bracket(1,1))

  !logX(m,r2) = X(m1,r2) + ((X(m2,r2) - X(m1,r2)) / (M(m2) - M(m1))) - M(m1,r1)
  logX_m_r2 = X_bracket(1,1) + (X_bracket(2,1) - X_bracket(1,1))/ (M_bracket(2,1) - M_bracket(1,1))

  logX_m_r = logX_m_r1 + (logX_m_r2 - logX_m_r1) / (R_bracket(1,2) - R_bracket(1,1)) * (logR_in - R_bracket(1,1))

end subroutine interpolate
