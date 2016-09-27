
subroutine moment_of_inertia( r,x,y,z,cellMass,parsec,radiusMin,lgradiusMin,lgradiusSph,bins,Ixxbin, Iyybin, & 
     Izzbin, Ixybin, Ixzbin, Izybin, size)
  implicit none
  integer, intent(in) :: bins, size
  double precision, intent(in) :: parsec,radiusMin,lgradiusMin,lgradiusSph
  double precision, dimension(bins), intent(out) :: Ixxbin, Iyybin, Izzbin, Ixybin, Ixzbin, Izybin
  double precision, dimension(size), intent(in) :: r, x, y, z, cellMass 
!f2py intent(in) bins, size
!f2py intent(in) parsec,radiusMin,lgradiusMin,lgradiusSph
!f2py intent(in) r, x, y, z, cellMass 
!f2py intent(out) Ixxbin, Iyybin, Izzbin, Ixybin, Ixzbin, Izybin
!f2py depend(size) r

  integer :: i, index 
  double precision :: Ixx, Iyy, Izz, Ixy, Ixz, Izy

  do i = 1, size 
     if( r(i)/parsec < radiusMin) continue
     
     index = int((log10(r(i)/parsec)-lgradiusMin)*bins/(lgradiusSph - lgradiusMin))

     if(index >= 0 .and. index < bins) then
        Ixx = (y(i)*y(i) + z(i)*z(i)) * cellMass(i)
        Iyy = (x(i)*x(i) + z(i)*z(i)) * cellMass(i)
        Izz = (x(i)*x(i) + y(i)*y(i)) * cellMass(i)
        Ixy = -1d0*x(i)*y(i) * cellMass(i)
        Ixz = -1d0*x(i)*z(i) * cellMass(i)
        Izy = -1d0*y(i)*z(i) * cellMass(i)
        
        Ixxbin(index) = Ixxbin(index) + Ixx
        Iyybin(index) = Iyybin(index) + Iyy
        Izzbin(index) = Izzbin(index) + Izz
        Ixybin(index) = Ixybin(index) + Ixy
        Ixzbin(index) = Ixzbin(index) + Ixz
        Izybin(index) = Izybin(index) + Izy
        
     end if
  end do

  return 
end subroutine moment_of_inertia
