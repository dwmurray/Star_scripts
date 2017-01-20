
program ps 
  use constants
  implicit none
  
  double precision :: tage, deltaM, deltaT, m, r, lum, n, md
  double precision :: T_core
  integer :: protostar_state
  double precision :: mdot

! Now proceed for the protostellar evolution.
  tage = 0.
  mdot = 1e-5*M_sun/yr
  deltaT = 1.e4*yr!0.01
  deltaM = mdot*deltaT

  m = 0.001 * M_sun
  protostar_state = 0

  do while (tage .le. 1.e7*yr)
     if(tage > 1e6*yr) deltaM = 0.
     call update_protostar( deltaM, deltaT, m, r, protostar_state)

!     call protostar(tage, deltaM, deltaT, m, r, lum, protostar_state)
     m = m + deltaM
     tage = tage + deltaT
     
!     write(10,102) tage/3.15e7, r, m
!102  format( 3(1pe10.4,2x))


  end do
  write(*,*) tage
  write(*,*) "m/m_sun = ", m/m_sun, "r/r_sun = ", r/r_sun
  write(*,*) r," ", lum
  
  stop
end program ps


subroutine update_protostar( deltaM, deltaT, m, r, protostar_state)!(deltaM, deltaT, m, r, n, L_D, protostar_state)
  use constants
  implicit none

  double precision, intent( in) :: deltaM, deltaT, m
  double precision, intent( inout) :: r!, n, md, T_core
  integer, intent( inout) :: protostar_state
  double precision :: L_MS, L_D, r_MS

!  double precision, parameter :: f_rad = 0.33 ! used to update to shell deuterium burning.

  if (protostar_state .EQ. 0) then
     if (m .LT. (0.01 * M_sun)) then   ! if the mass is below 0.01 solar masses then we do not initialize
        protostar_state = 0   ! This is considered the pre-collapse state.
        return
     else ! if it is now above 0.01 solar masses, it is in the no burning state.
        protostar_state = 1
        call initialize_protostar( deltaM, deltaT, m, r) ! sets initial r, n, md for the particle.
        return
     end if
  end if

  if (protostar_state .EQ. 1) then
     call update_radius( deltaM, deltaT, m, r)
     return
  end if
end subroutine update_protostar


!subroutine protostar( tage, deltaM, deltaT, m, r, lum, protostar_state)
!  use constants
!  implicit none
!  
!  double precision, intent( in) :: tage, deltaM, deltaT, m
!
!  double precision, intent( out) :: lum
!  double precision, intent( inout) :: r!, n, md, T_core
!  integer, intent( inout) :: protostar_state
!
!  double precision :: deltaR, L_D, L_MS
!
!  call update_radius( deltaM, deltaT, m, r)
!  return
!end subroutine protostar

subroutine initialize_protostar( deltaM, deltaT, m, r)
  use constants
  implicit none

  double precision, intent( in) :: deltaM, deltaT, m
  double precision, intent( out) :: r

  !init is from Offner et al. 2009 ApJ 703
  r = ( 2.5 * R_sun) * ( deltaM * yr / (deltaT * 1.E-5 * M_sun))**0.2

!  n = 5. - 3. * (1.475 + 0.07 * log10( deltaM / deltaT))**(-1.)
!  if (n .LE. 1.5) then
!     n = 1.5
!  else if (n .GE. 3) then
!     n = 3
!  end if
!
!  md = m
!  T_core = 0.

  return
end subroutine initialize_protostar

subroutine update_radius( deltaM, deltaT, m, r)
  use constants
  implicit none

  double precision, intent( in) :: deltaM, deltaT, m
  double precision, intent( inout) :: r

  double precision :: L_I, L_D
  double precision, parameter :: f_acc=1.0 ! This is a correction factor in Nakano 95.
  double precision, parameter :: a_e = 3./4. ! 3/4. is for n = 3
  double precision, parameter :: sigma = 5.67D-5 !stefan-Boltzmann cgs.
  double precision :: T_eff
  double precision :: L_star, R_star ! ZAMS values
  double precision :: mdot, beta
  double precision :: deltaR

  call set_ZAMS( m + deltaM, L_star, R_star, beta)

  !write(*,*) "called set_ZAMS", m+deltaM, L_star, R_star, beta  
  mdot = deltaM/deltaT
  L_D = 15. * L_sun * mdot/( 1.e-5 * M_sun/yr)
  L_I = 2.5 * L_sun * mdot/( 1.e-5 * M_sun/yr) !dot{M} / 1e-5 M_sun per yr

  ! Nakano et al. 95 ApJ 450...183N Appendix A equation 27

  !r = G*a_e*(m+deltaM)*mdot/L_star*(1. + beta - (1.+f_acc)/a_e*0.5)
  deltaR = (2. - (1.+f_acc)/(2.*a_e)*r*deltaM/m - r/(a_e*G*m*m)*(L_star + L_I - L_D))*deltaT
  r = r + deltaR 
  !r = r + deltaR
  !deltaR = (2. - (1. + f_acc) / (2.*a_e)) * r * deltaM / m - r*r*deltaT*L_star * deltaM / (a_e * G * deltaM*m*m) &
  !     - r*r*deltaM / (a_e*G*m*m) * (L_I*deltaT / deltaM - L_D*deltaT / deltaM)
  
  if( r < R_star) then ! conditions to ZAMS
     r = R_star
  end if
!  write(*,*) "new mass radius", m+deltaM, r/r_sun, G*a_e*(m+deltaM)*mdot/L_star/r_sun, beta, R_star/r_sun
!  write(12,103) "L_* ",L_star , &
!       "L_I ", L_I , &
!       "L_D ", L_D , &
!       "f_acc ", (2. - (1. + f_acc) / (2.*a_e)), &
!       "RHS ", r*deltaT / (a_e * G * deltaM*m), &
!       "LHS ", r*deltaM /m
!103 format( 6(A6,(1pe10.4,2x)))

  return
end subroutine update_radius

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
     lgL = lgl_zams(1) - beta_zams(1)*(lgm - lgm_zams(1))
     lgR = lgr_zams(1) - betaR_zams(1)*(lgm - lgm_zams(1))
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
