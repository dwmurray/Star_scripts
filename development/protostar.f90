
program ps 
  use constants
  implicit none
  
  double precision :: tage, deltaM, deltaT, m, r, lum, n, md
  double precision :: T_core
  integer :: protostar_state
!  double precision, parameter :: M_solar = 1.989D33
!  double precision, parameter :: R_solar = 6.96D10

! Initialize the table of phi and epsilon
  integer, parameter :: phi_dim1=4, phi_dim2=50
  integer, parameter  :: proto_dim1=40, proto_dim2=40
  double precision, dimension(phi_dim1, phi_dim2) :: phi_eps_table
  double precision, dimension(proto_dim1, proto_dim2) :: proto_core_table
  double precision :: mdot
!  call init_phi_eps(phi_dim1, phi_dim2, phi_eps_table)

! Now that we have phi and eps, we need to get rho_core & P_core
! From P_core and rho_core we can root find for T_core and Beta.
!  call init_P_and_rho_core(phi_dim1, phi_dim2, phi_eps_table, proto_dim1, proto_dim2, proto_core_table)

! Now proceed for the protostellar evolution.
  tage = 0.
  mdot = 1e-5*M_sun/yr
  deltaT = 1.e4*yr!0.01
  deltaM = mdot*deltaT

  m = 0.001 * M_sun
  protostar_state = 0
!  n = 5./3.
!  md = 0.

  do while (tage .le. 1.e7*yr)
     if(tage > 1e6*yr) deltaM = 0.
     call update_protostar_state( deltaM, deltaT, m, r, n, md, T_core, yr, protostar_state, phi_dim1, phi_dim2, phi_eps_table)

     call protostar(tage, deltaM, deltaT, m, r, lum, n, md, T_core, protostar_state, phi_dim1, phi_dim2, phi_eps_table)
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


subroutine update_protostar_state( deltaM, deltaT, m, r, n, md, T_core, yr, protostar_state, phi_dim1, phi_dim2, phi_eps_table)!(deltaM, deltaT, m, r, n, L_D, protostar_state)
  implicit none

  double precision, intent( in) :: deltaM, deltaT, m, yr
  integer, intent(in) :: phi_dim1, phi_dim2
  double precision, intent( in), dimension(phi_dim1,phi_dim2) :: phi_eps_table

  double precision, intent( inout) :: r, n, md, T_core
  integer, intent( inout) :: protostar_state
  double precision :: L_MS, L_D, r_MS

  double precision, parameter :: M_solar = 1.989D33
  double precision, parameter :: R_solar = 6.96D10
  double precision, parameter :: f_rad = 0.33 ! used to update to shell deuterium burning.

!  write (*,*) phi_eps_table(:,1)

  if (protostar_state .EQ. 0) then
     if (m .LT. (0.01 * M_solar)) then   ! if the mass is below 0.01 solar masses then we do not initialize
        !write (*,*) "m less than 0.01 Msolar"
        protostar_state = 0   ! This is considered the pre-collapse state.
        return
     else ! if it is now above 0.01 solar masses, it is in the no burning state.
        protostar_state = 1
        call initialize_protostar( deltaM, deltaT, m, r, n, md, T_core, yr) ! sets initial r, n, md for the particle.
        return
     end if
  end if

  if ( protostar_state .EQ. 1) then
     ! Find Tcore
!     call set_T_core()
     if (T_CORE .LT. 1.5D6) then
        protostar_state = 1
        return
     else
        protostar_state = 2 ! update to core burning at fixed T_core
        return
     end if
  end if
  
  if (protostar_state .EQ. 2) then
     call update_md(deltaM, deltaT, L_D, md)
     if (md .GT. 0.0) then
        protostar_state = 2
        return
     else
        protostar_state = 3 ! update to core burning at variable T_core
        return
     end if
  end if

  if (protostar_state .EQ. 3) then
     !compare L_D to L_ms
     call set_L_D()
     call set_L_MS(L_MS)
     if (L_D/L_MS .LT. f_rad) then
        protostar_state = 3
        return
     else
        protostar_state = 4 ! update to shell burning deuterium
        return
     end if
  end if

  if (protostar_state .EQ. 4) then
     ! Compare r to r_ms
     call set_r_MS()
     if (r .GT. (r_ms + (0.05 * r_ms))) then ! Need them to be close, not exactly equal
        protostar_state = 4
        return
     else
        protostar_state = 5 ! update to main sequence star.
        return
     end if
  end if

end subroutine update_protostar_state


subroutine protostar( tage, deltaM, deltaT, m, r, lum, n, md, T_core, protostar_state, phi_dim1, phi_dim2, phi_eps_table) 
  use constants
  implicit none
  
  double precision, intent( in) :: tage, deltaM, deltaT, m!, yr
  integer, intent(in) :: phi_dim1, phi_dim2
  double precision, intent( in), dimension(phi_dim1,phi_dim2) :: phi_eps_table

  double precision, intent( out) :: lum
  double precision, intent( inout) :: r, n, md, T_core
  integer, intent( inout) :: protostar_state

  double precision :: deltaR, L_D, L_MS

  call update_radius( deltaM, deltaT, m, r)
  

  return
end subroutine protostar

subroutine initialize_protostar( deltaM, deltaT, m, r, n, md, T_core, yr)
  implicit none

  double precision, intent( in) :: deltaM, deltaT, m, yr
  double precision, intent( out) :: r, n, md, T_core

  double precision, parameter :: M_solar = 1.989D33
  double precision, parameter :: R_solar = 6.96D10

  !init is from Offner et al. 2009 ApJ 703
  r = ( 2.5 * R_solar) * ( deltaM * yr / (deltaT * 1.E-5 * M_solar))**0.2

  n = 5. - 3. * (1.475 + 0.07 * log10( deltaM / deltaT))**(-1.)
  if (n .LE. 1.5) then
     n = 1.5
  else if (n .GE. 3) then
     n = 3
  end if

  md = m
  T_core = 0.

  return
end subroutine initialize_protostar

!double check this is all we require.
!subroutine update_radius_Offner( deltaM, deltaT, m, r, n, L_D, L_MS, deltaR)
!  implicit none
!
!  double precision, intent( in) :: deltaM, deltaT, m, r, n, L_D, L_MS
!  double precision, intent( out) :: deltaR
!  double precision :: alpha_g! alpha_g describes the gravitational binding energy of a poltrope.
!  double precision :: L_H, L_int, L_ion ! luminosity of hayashi track star of some radius r.
!
!  double precision, parameter :: L_solar = 3.839D33
!  double precision, parameter :: M_solar = 1.989D33
!  double precision, parameter :: f_k = 0.5 !fraction of the kinetic energy of the infalling material that is radiated away.
!  ! Do need to look at how Offner et al. calculate Beta
!  double precision, parameter :: Beta = 1  !Beta is the mean ratio of the gas pressure to total gas plus radiation pressure in the star.
!  double precision, parameter :: Temp_H = 3000 ! surface temp of hayashi track star.
!  double precision, parameter :: sigma = 5.67D-5 !stefan-Boltzmann cgs.
!  double precision, parameter :: pi = 3.1415
!
!
!  alpha_g = 3./(5. - n)
!  write(*,*) alpha_g
!
!
!  L_H = 4.*pi*r**2*sigma*Temp_H**4
!  L_int = MAX(L_MS, L_H)  ! Need L_hayashi & L_ms to find the interior luminosity
!  
!  L_ion = 2.5 * L_solar * (deltaM/deltaT)/(1d-5 * M_solar/(pi*1D7))
!! Offner et al. deltaR
!  deltaR = 2.* (deltaM/m) * (1 - (1-f_k)/(alpha_g*Beta) + (1./2.)*(dlogBeta / dlogm)) * r - 2. * (deltaT / (alpha_g*Beta)) * (r/Gm*m)*(L_int + L_ion - L_D)*r
!
!  return
!end subroutine update_radius_Offner


subroutine update_radius( deltaM, deltaT, m, r)
  ! Nakano et al. 95 ApJ 450...183N Appendix A equation 27
  use constants

  implicit none


  double precision, intent( in) :: deltaM, deltaT, m
  double precision, intent( inout) :: r

  double precision :: L_I, L_D

  double precision, parameter :: f_acc=1.0 ! This is a correction factor in Nakano 95. !TODO obtain this value
  double precision, parameter :: a_e = 3./4. ! 3/4. is for n = 3
!  double precision, parameter :: L_solar = 3.839D33
!  double precision, parameter :: M_solar = 1.989D33
!  double precision, parameter :: f_k = 0.5 !fraction of the kinetic energy of the infalling material that is radiated away.
!  double precision, parameter :: Temp_H = 3000 ! surface temp of hayashi track star.
  double precision, parameter :: sigma = 5.67D-5 !stefan-Boltzmann cgs.
  double precision :: T_eff
  double precision :: L_star, R_star ! ZAMS values
  double precision :: mdot, beta
  double precision :: deltaR
!  write(*,*) G, pi, a_e, f_acc

!  L_H = 4.*pi*r**2*sigma*Temp_H**4
  
! Rough estimate for L_star currently !TODO Update this estimate

  call set_ZAMS( m + deltaM, L_star, R_star, beta)
!  L_star = set_L_star(m, beta)!  L_star = L_sun * m*m*m / m_sun**3
  !write(*,*) "called set_ZAMS", m+deltaM, L_star, R_star, beta  
  mdot = deltaM/deltaT
  L_D = 15. * L_sun * mdot/( 1.e-5 * M_sun/yr)
  L_I = 2.5 * L_sun * mdot/( 1.e-5 * M_sun/yr) !dot{M} / 1e-5 M_sun per yr


  !r = G*a_e*(m+deltaM)*mdot/L_star*(1. + beta - (1.+f_acc)/a_e*0.5)
  deltaR = (2. - (1.+f_acc)/(2.*a_e)*r*deltaM/m - r/(a_e*G*m*m)*(L_star + L_I - L_D))*deltaT
  r = r + deltaR 
  !r = r + deltaR
  !deltaR = (2. - (1. + f_acc) / (2.*a_e)) * r * deltaM / m - r*r*deltaT*L_star * deltaM / (a_e * G * deltaM*m*m) &
  !     - r*r*deltaM / (a_e*G*m*m) * (L_I*deltaT / deltaM - L_D*deltaT / deltaM)
  
  if( r < R_star) then ! conditions to ZAMS
     r = R_star
  end if
  write(*,*) "new mass radius", m+deltaM, r/r_sun, G*a_e*(m+deltaM)*mdot/L_star/r_sun, beta, R_star/r_sun
  write(12,103) "L_* ",L_star , &
       "L_I ", L_I , &
       "L_D ", L_D , &
       "f_acc ", (2. - (1. + f_acc) / (2.*a_e)), &
       "RHS ", r*deltaT / (a_e * G * deltaM*m), &
       "LHS ", r*deltaM /m
103 format( 6(A6,(1pe10.4,2x)))

!  deltaR = ((2. - (1. + f_acc) / (2.*a_e)) - r*deltaT*L_star / (a_e * G * deltaM*m) &
!       - r / (a_e*G*m) * (L_I*deltaT / deltaM - L_D*deltaT / deltaM)) * (r*deltaM /m)

  !deltaR =  - r*r*deltaT*L_star * deltaM / (a_e * G * deltaM*m*m) &
  !deltaR = - r*r*deltaM / (a_e*G*m*m) * (L_I*deltaT / deltaM - L_D*deltaT / deltaM)
  !T_eff = (L_star / (4.*pi*r*r*sigma_b))**0.25
!  write(11,101) "T_eff", T_eff, "L_star", L_star/L_sun, "M_star", M/M_sun, "R_star", r/r_sun, "deltaR/r_sun", deltaR/r_sun
!101 format( 4(A6,(1pe10.3,2x)), A12,(1pe10.3,2x))

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

function set_L_star(m, beta)
  use constants
  implicit none

  double precision :: set_L_star
  double precision, intent( in) :: m, beta ! Beta is power law slope


  set_L_star = L_sun * (m / m_sun)**beta

end function set_L_star

subroutine set_L_MS(L_MS)
  implicit none
  double precision, intent( inout) :: L_MS

  !L_ms = 3.*L_solar !See Trout et al. 1996
end subroutine set_L_MS

subroutine set_L_D()
  implicit none

  return
end subroutine set_L_D

subroutine update_md(deltaM, deltaT, L_D, md)
  implicit none
  double precision, intent( in) :: deltaM, deltaT, L_D
  double precision, intent( inout) :: md
  double precision :: deltaMd

  double precision, parameter :: M_solar = 1.989D33
  double precision, parameter :: L_solar = 3.839D33

  deltaMd = deltaM - (1.D-5 * M_solar * (L_D / (15.*L_solar)) * deltaT)

  md = md - deltaMd

  if (md .LT. 0.0) then
     md = 0.0
  end if
  return
end subroutine update_md

subroutine set_r_MS()
  implicit none
  !See Trout et al. 1996 on how to set r_ZAMS

  return
end subroutine set_r_MS

subroutine init_phi_eps(phi_dim1, phi_dim2, phi_eps_table)
  use rk2_polytrope
  implicit none
  integer, intent(in) :: phi_dim1, phi_dim2
  double precision, intent( inout), dimension(phi_dim1,phi_dim2) :: phi_eps_table
  double precision :: n, eps1, dphideps1
  integer :: i
  double precision, parameter :: nstart=1.4, nend=3.1

!  n = 1.9
!  call polytrope( n, eps1, dphideps1) 
!  do i =1, 1000
!     phi_eps_table(1,i) = n
!     phi_eps_table(2,i) = eps1
!     phi_eps_table(3,i) = dphideps1
!  end do
  do i =1 , phi_dim2

     n = nstart + (nend-nstart)*i/phi_dim2 !This is 1001 steps across 3.1-1.4! (3.1-1.4)/1000.!0.01

     call polytrope( n, eps1, dphideps1) 
!     write (*,*) n, eps1, -eps1**2*dphideps1
     phi_eps_table(1,i) = n
     phi_eps_table(2,i) = eps1
     phi_eps_table(3,i) = dphideps1
     ! trying to include a fourth column gives a segmentation fault.
     phi_eps_table(4,i) = -eps1**2*dphideps1
     !write (*,*) phi_eps_table(:,i)
  end do
!  write (*,*) i
!  write (*,*) phi_eps_table
end subroutine init_phi_eps

subroutine init_P_and_rho_core(phi_dim1, phi_dim2, phi_eps_table, proto_dim1, proto_dim2, proto_core_table)
  use protostellar_core
  use constants

  implicit none
  
  integer, intent(in) :: phi_dim1, phi_dim2
  double precision, intent( in), dimension(phi_dim1,phi_dim2) :: phi_eps_table

  integer, intent(in)            :: proto_dim1, proto_dim2  !proto_core_table !n, m, rho, P, T, Beta
  double precision, intent( out), dimension(proto_dim1, proto_dim2)    :: proto_core_table

  double precision    :: rho_core_value, P_core_value
  double precision    :: M, R
  double precision    :: n, eps1, dphideps1

  integer             :: i,j,k

  ! for a solar mass then scale linearly.
  !R = ( 2.5 * R_sun) !* ( deltaM / (deltaT * 1.E-5 * M_solar))**0.2
  
  i = 1
  do i = 1, phi_dim2
     n = phi_eps_table(1,i)
     eps1 = phi_eps_table(2,i)
     dphideps1 = phi_eps_table(3,i)

     M = 0.01*M_sun
     do j = 1, 600

        R = 0.1*R_sun
        do k = 1, 100
           rho_core_value = determine_rho_core(M, R, eps1, dphideps1)
           P_core_value = determine_P_core(M, R, n, dphideps1)
           if (i .EQ. 5) then
              write(*,101) i, j, n, eps1, dphideps1, M, rho_core_value, P_core_value
101           format(2i8, 6(1pe10.3,2x))

              !write out to a csv and then plot with python to make sure there are no bumps etc.
              !how big is a float? and how many floats per gigabyte?
              
              !Get T and Beta now as well.
           end if
           R = R + 0.1*R_sun
        end do
        M = M + 0.01*M_sun
     end do
  end do


end subroutine init_P_and_rho_core

subroutine init_T_core_Beta(phi_dim1, phi_dim2, phi_eps_table)
  use root_find
  
  integer, intent(in) :: phi_dim1, phi_dim2
  double precision, intent( in), dimension(phi_dim1,phi_dim2) :: phi_eps_table

  double precision    :: rho_core_value, P_core_value
  double precision    :: M, R
  double precision    :: n, eps1, dphideps1

  !call find_root(n, rho_core, P_core, xinit, machine_tolerance, maxiter, result)

end subroutine init_T_core_Beta
