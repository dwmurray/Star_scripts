program test_machine_epsilon

  use constants, only      : r_15
  use machine_epsilon

  implicit none

  real(r_15)        :: smallest_number

  smallest_number = determine_machine_epsilon()
  write(*,*)'machine epsilon = ', smallest_number


end program test_machine_epsilon
