module grid

  implicit none

  private
  save

  public :: setup_grid

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!m layers from 1 to m, 1 and m are ghost points at distance dz/2 from boundaries (inside the domain).
!Surface is near k = 1, bottom of M.L. near k = m
!surface at "depth" = 0
!total depth = (m-1)*dz
!first point dz/2
!last point (m-1/2)*dz

!-------------------------------------------------------------------------------------------------------------------------!

subroutine setup_grid(m,dz,h,z)
implicit none
integer, intent(in) :: m
double precision, intent(in) :: dz
double precision, intent(out) ::h,z(:)
integer :: k

h = dz*float(m-1)

do k = 1,m
  z(k) = (float(k)-0.5)*dz
end do

end subroutine setup_grid

!-------------------------------------------------------------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module grid
