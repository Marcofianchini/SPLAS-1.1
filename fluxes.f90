module fluxes


implicit none

private
save

integer, parameter :: dim=70
integer, parameter :: dn=3
integer, parameter :: m=2

public :: p_diff,p_sink_input,nut_mix,n_load

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-------------------------------------------------------------------------------------------------------------------------!

subroutine p_diff(D,dz,p,p_diff_flux) !phytoplankton eddy diffusion (vertical)
implicit none
double precision, intent(in) :: D(m),dz,p(dim,m)
double precision, intent(out) :: p_diff_flux(dim,m)
double precision :: d2
integer :: k

d2 = (1.0/(dz*dz))

if (m == 2) then
  p_diff_flux(:,1) = D(1)*(p(:,m) - p(:,1))
  p_diff_flux(:,m) = D(m)*(p(:,1) - p(:,m))
else
  do k = 2,m-1
    p_diff_flux(:,k) = D(k)*d2*(p(:,k+1) -2*p(:,k) + p(:,k-1))
  end do

  !p(:,1) = p(:,2)
  !p(:,m) = p(:,m-1)

endif

end subroutine p_diff

!-------------------------------------------------------------------------------------------------------------------------!

subroutine p_sink_input(sig,p,p_sink_flux) !phytoplankton sinking input from the layer above
implicit none
double precision, intent(in) :: sig(dim),p(dim,m)
double precision, intent(out) :: p_sink_flux(dim,m)
integer :: k

if (m == 2) then
  p_sink_flux(:,1) = -sig(:)*p(:,1)
  p_sink_flux(:,m) = sig(:)*(p(:,1)-p(:,m))
else
  do k = 2,m-1
    p_sink_flux(:,k) = sig(:)*(-p(:,k)+p(:,k-1))
  end do

  !p(:,1) = p(:,2)
  !p(:,m) = p(:,m-1)

endif

end subroutine p_sink_input

!-------------------------------------------------------------------------------------------------------------------------!

subroutine nut_mix(kmix,S,n,nut_flux)
implicit none
double precision, intent(in) :: kmix,S,n(m)
double precision, intent(out) :: nut_flux
integer :: k

if (m == 2) then
  nut_flux(1) = kmix*(n(m)-n(1))
  nut_flux(m) = S
else
  do k = 2,m-1
  nut_flux(k) = 0.0
  end do

  nut_flux(m) = S

  !n(1) = n(2)
  !n(m) = n(m-1)

endif

end subroutine nut_mix

!-------------------------------------------------------------------------------------------------------------------------!

subroutine n_load()
return
end subroutine n_load

!-------------------------------------------------------------------------------------------------------------------------!


!-------------------------------------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------------!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end module fluxes
