module protostar_parameters
  use amr_parameters, only : dp

  implicit none
  !For determining Zero Age Main Sequence values.
  real(dp), dimension(100) :: lgm_zams, lgl_zams, lgr_zams, beta_zams, betaR_zams
  integer :: num_zams
  logical :: initialized_zams = .false.
  logical :: initialized_protostar = .false.
  !Creating the table of Beta values to interpolate across for given N, R and M.
  integer, parameter :: n_step = 10
  real(dp), parameter :: logM_step = 0.1, logR_step = 0.1
  real(dp), parameter :: n_min = 1.5, n_max = 3.3
  real(dp), parameter :: M_min = 1.0D-2, M_max = 1.0D2 ! In Solar masses.
  real(dp), parameter :: R_min = 1.0D-2, R_max = 1.0D3 ! In Solar radii

  !Used in determining DeltaR in the Offner prescription (udpate_radius call)
  !As well as in setting the deuterium luminosity
  real(dp), parameter :: f_k = 0.5 !fraction of the kinetic energy of the infalling material that is radiated away.

  !Used in set_L_internal Hayashi track temperature.
  real(dp), parameter :: T_Hayashi = 3000.

  !number of interations before admitting defeat in root finding.
  integer, parameter :: find_root_maxiter = 200, phi_eps_iter = 100

  !Used in determine_T_core, and create_beta_table_file
  real(dp), parameter :: mu = 0.613 !Offner et al. '09 mean molecular weight for fully ionized gas of solar comp.
  real(dp), parameter :: machine_tolerance = 1e-14, Beta_guess = 0.5

  !Used in obtain_dlogBetadlogM
  real(dp), parameter :: deriv_step = 1.D-1

  integer  :: beta_table_dim1
  integer, parameter :: beta_table_dim2 = 7!dim2 is n, M, R, rho_core, P_core, Beta_core, Beta_mean

  !this array will hold the table to interpolate beta from.
  real(dp),allocatable,dimension(:,:) :: Beta_file_table

end module protostar_parameters
