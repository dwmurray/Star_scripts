module constants
  ! From "The Fundamental Physical Constants, E. R. Cohen & B. N. Taylor, Physics Today
  ! p9, August 1990
  implicit none
  
  integer, parameter    :: r_15 = selected_real_kind(15,307)
  real (r_15),parameter :: pi = 3.14159265358979d0
  real (r_15),parameter :: G = 6.67259d-8, c = 2.99792458d10 , e = 4.803206799d-10
  real (r_15),parameter :: h = 6.6260755d-27, h_bar = 1.05457266d-27
  real (r_15),parameter :: m_p = 1.6726231d-24, m_e = 9.1093897d-28
  real (r_15),parameter :: k_b=1.380658d-16
  real (r_15),parameter :: eV = 1.60217733d-12
  real (r_15),parameter :: nu_lyman = 3.2898419499d15            !Lyman edge
  real (r_15),parameter :: f_alpha = 0.415d0, f_beta = 0.119d0   !Lyman alpha, H beta
  real (r_15),parameter :: a_edge = 6.3d-18                      !Lyman edge, Osterbrock
  real (r_15),parameter :: sigma_b = pi**2 * k_b**4 / (60.d0 * h_bar**3 * c**2)
  real (r_15),parameter :: sigma_thomp = 6.65e-25, kappa_es = sigma_thomp/m_p
  real (r_15),parameter :: a_rad = 7.5657D-15 !erg cm-3 K-4

  real (r_15),parameter :: parsec = 3.08e18, M_sun = 2.d33, L_sun = 3.85d33
  real (r_15),parameter :: R_sun = 6.95949d10
  real (r_15),parameter :: degrees_to_radians = pi / 180.d0
  real (r_15),parameter :: AU = 1.5d13, R_J = 7.d9, M_J = 1.91d30, year = 3.15e7
  real (r_15),parameter :: C_Spitzer = 1.2d-6   !Spitzer conductivity
  real (r_15),parameter :: yr = 24.*3600.*365.24 
  real (r_15), dimension(100) :: lgm_zams, lgl_zams, lgr_zams, beta_zams, betaR_zams
  integer :: num_zams
  logical :: initialized_zams = .false.
  
end module constants
