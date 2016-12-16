module machine_epsilon
!--------------------------------------------------------------------------
! Identify the smallest difference the current machine can find. Subtracting
! two number that differ by less than this amount will generate noise
!--------------------------------------------------------------------------
  use constants, only      : r_15
  implicit none

contains

  function determine_machine_epsilon()
    real(r_15)              :: determine_machine_epsilon
    real(r_15)              :: u, den

    integer                 :: i, large_integer = 100

!    write(*,*)r_15
  
    u = 1.d0
    do i =1, large_integer
       u = u / 2.d0
       den = 1.0d0 + u
!       write(*,*)'i=', i, 'den=', den
       if(den .le. 1.d0) exit
    end do

    if(i .gt. large_integer) then
       write(*,*)'cannot determine machine epsilon, quitting'
       stop
    end if

    determine_machine_epsilon = 10.d0 * u
!    write(*,*)'machine epsilon =', determine_machine_epsilon
  end function determine_machine_epsilon
  
end module machine_epsilon
