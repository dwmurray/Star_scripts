program test_newton

  use root_find
  
  implicit none
  
  double precision      :: x! = -3.
  double precision      :: xinit! = 1.D-7
  double precision      :: n
  double precision      :: root
  double precision      :: func, func_prime
  integer               :: maxiter = 100
  integer               :: j
  double precision      :: machine_tolerance = 1D-15
  
  !  write(*,*) 'test, should find 0'
  ! Writing out to find out what the bracketed region is.
  x = -3.
  do j=1, 100
     write(9,*) x, func(x), func_prime(x)
     x = x + 6./100.
  end do
  
  xinit = 1.D-7

  !for now, just providing n, rho_core and P core. future will have func calls
  n = 3


  call find_root(n, rho_core, P_core, xinit, machine_tolerance, maxiter, root)
  write(*,*) root
  
end program test_newton
