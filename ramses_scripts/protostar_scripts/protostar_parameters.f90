module protostar_parameters

  implicit none
  !For determining Zero Age Main Sequence values.
  integer, parameter    :: r_15 = selected_real_kind(15,307)
  real (r_15), dimension(100) :: lgm_zams, lgl_zams, lgr_zams, beta_zams, betaR_zams
  integer :: num_zams
  logical :: initialized_zams = .false.

  !Creating the table of Beta values to interpolate across for given N, R and M.
  integer, parameter :: n_step = 10
  double precision, parameter :: logM_step = 0.1, logR_step = 0.1
  double precision, parameter :: n_min = 1.5, n_max = 3.3
  double precision, parameter :: M_min = 1.0D-2, M_max = 1.0D2 ! In Solar masses.
  double precision, parameter :: R_min = 1.0D-2, R_max = 1.0D3 ! In Solar radii

  !Used in determining DeltaR in the Offner prescription (udpate_radius call)
  !As well as in setting the deuterium luminosity
  double precision, parameter :: f_k = 0.5 !fraction of the kinetic energy of the infalling material that is radiated away.

  !Used in set_L_internal Hayashi track temperature.
  double precision, parameter :: T_Hayashi = 3000.

  !number of interations before admitting defeat in root finding.
  integer, parameter :: find_root_maxiter = 200, phi_eps_iter = 100

  !Used in determine_T_core, and create_beta_table_file
  double precision, parameter :: mu = 0.613 !Offner et al. '09 mean molecular weight for fully ionized gas of solar comp.
  double precision, parameter :: machine_tolerance = 1e-14, Beta_guess = 0.5

  !Used in obtain_dlogBetadlogM
  double precision, parameter :: deriv_step = 1.D-1


end module protostar_parameters
