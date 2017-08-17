module protostar_constants
  use amr_parameters, only : dp
  ! From "The Fundamental Physical Constants, E. R. Cohen & B. N. Taylor, Physics Today
  ! p9, August 1990
  implicit none
  
  real (dp),parameter :: pi = 3.14159265358979d0
  real (dp),parameter :: G = 6.67259d-8, c = 2.99792458d10 , e = 4.803206799d-10
  real (dp),parameter :: h = 6.6260755d-27, h_bar = 1.05457266d-27
  real (dp),parameter :: m_p = 1.6726231d-24, m_e = 9.1093897d-28
  real (dp),parameter :: k_b=1.380658d-16
  real (dp),parameter :: eV = 1.60217733d-12
  real (dp),parameter :: nu_lyman = 3.2898419499d15            !Lyman edge
  real (dp),parameter :: f_alpha = 0.415d0, f_beta = 0.119d0   !Lyman alpha, H beta
  real (dp),parameter :: a_edge = 6.3d-18                      !Lyman edge, Osterbrock
  real (dp),parameter :: sigma_b = pi**2 * k_b**4 / (60.d0 * h_bar**3 * c**2)
  real (dp),parameter :: sigma_thomp = 6.65e-25, kappa_es = sigma_thomp/m_p
  real (dp),parameter :: a_rad = 7.5657D-15 !erg cm-3 K-4

  real (dp),parameter :: parsec = 3.08e18, M_sun = 2.d33, L_sun = 3.85d33
  real (dp),parameter :: R_sun = 6.95949d10
  real (dp),parameter :: degrees_to_radians = pi / 180.d0
  real (dp),parameter :: AU = 1.5d13, R_J = 7.d9, M_J = 1.91d30!, year = 3.15e7
  real (dp),parameter :: C_Spitzer = 1.2d-6   !Spitzer conductivity
  real (dp),parameter :: yr = 24.*3600.*365.24 

end module protostar_constants
