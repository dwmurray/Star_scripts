
program ps 
  implicit none
  
  double precision :: tage, deltaM, deltaT, m, r, lum, n, md
  double precision, parameter :: M_solar = 1.989D33
  double precision, parameter :: R_solar = 6.96D10
  
  tage = 0.
  deltaM = 0.01
  deltaT = 0.01
  m = 0.02 * M_solar
  n = 5./3.
  md = 0.

  call protostar(tage, deltaM, deltaT, m, r, lum, n, md)

  write(*,*) r," ", lum
  
  stop
end program ps

subroutine protostar( tage, deltaM, deltaT, m, r, lum, n, md) 
  implicit none
  
  double precision, intent( in) :: tage, deltaM, deltaT, m
  double precision, intent( out) :: r, lum
  double precision, intent( inout) :: n, md

  double precision, parameter :: M_solar = 1.989D33
  double precision, parameter :: R_solar = 6.96D10

  ! if the mass is below 0.01 solar masses then we do not initialize
  ! This is considered the pre-collapse state.
  if (m .LT. (0.01 * M_solar)) then
     write (*,*) "m less than 0.01 Msolar"
     r = 0.
     lum = 0.
     n = 0.
     return
  end if
  write (*,*) "m = ", m, 'and r = ', r
  ! First time step where m > 0.01 we set r and n and md
  ! This is the no burning state.
  if ((m .GE. (0.01 * M_solar)) .and. ( r .EQ. 0.)) then
     write (*,*) "initializing radius of protostar."
     call initialize_protostar( tage, deltaM, deltaT, m, r, lum, n, md)
     write(*,*) n, md
     !write(*,*) r," ", lum
  end if

  if ((m .GE. (0.01 * M_solar)) .and. ( r .GT. 0.)) then
     write (*,*) "updating radius of protostar."
     call update_protostar( deltaM, deltaT, m, r, deltaR) !double check this is all we require.
     
  end if
!  r = 0.
!  lum = 0.
  return
end subroutine protostar

subroutine initialize_protostar( tage, deltaM, deltaT, m, r, lum, n, md)
  implicit none

  double precision, intent( in) :: tage, deltaM, deltaT, m
  double precision, intent( out) :: r, lum, md
  double precision, intent( inout) :: n

  double precision, parameter :: M_solar = 1.989D33
  double precision, parameter :: R_solar = 6.96D10

  r = ( 2.5 * R_solar) * ( deltaM / (deltaT * 1.E-5 * M_solar))**0.2

  n = 5. - 3. * (1.475 + 0.07 * log10( deltaM / deltaT))**(-1.)

  md = m

  return
end subroutine initialize_protostar
