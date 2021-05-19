module euler_modified


  implicit none
  public
  save

  integer, parameter :: dim=70
  integer, parameter :: dn=3

  public :: nfirst_step,nsecond_step,pfirst_step,psecond_step

contains

  !---------------------------------------------------------------------------------------------------------------------!

  subroutine nfirst_step(pred,old,dt,fo)

  implicit none

  double precision, dimension(dn), intent(in) :: old
  double precision, dimension(dn), intent(in) :: fo
  double precision, intent(in) :: dt
  double precision, dimension(dn), intent(out) :: pred

  pred = old + 0.5*dt*fo

  !FIRST STEP : EXPLICIT EULER with half time step

end subroutine nfirst_step

  !---------------------------------------------------------------------------------------------------------------------!


  subroutine nsecond_step(pred,old,dt,fp,new)

  implicit none

  double precision, dimension(dn), intent(in) :: pred,old
  double precision, dimension(dn), intent(in) :: fp
  double precision, intent(in) :: dt
  double precision, dimension(dn), intent(out) :: new

  new = old + dt*fp

  !SECOND STEP : IMPLICIT TRAPEZOIDAL RULE with full time step

end subroutine nsecond_step

  !---------------------------------------------------------------------------------------------------------------------!
  !---------------------------------------------------------------------------------------------------------------------!
  !---------------------------------------------------------------------------------------------------------------------!
  !---------------------------------------------------------------------------------------------------------------------!
  !---------------------------------------------------------------------------------------------------------------------!


  !---------------------------------------------------------------------------------------------------------------------!

  subroutine pfirst_step(pred,old,dt,fo)

  implicit none

  double precision, dimension(dim), intent(in) :: old
  double precision, dimension(dim), intent(in) :: fo
  double precision, intent(in) :: dt
  double precision, dimension(dim), intent(out) :: pred

  pred = old + 0.5*dt*fo

  !FIRST STEP : EXPLICIT EULER with half time step

end subroutine pfirst_step

  !---------------------------------------------------------------------------------------------------------------------!


  subroutine psecond_step(pred,old,dt,fp,new)

  implicit none

  double precision, dimension(dim), intent(in) :: pred,old
  double precision, dimension(dim), intent(in) :: fp
  double precision, intent(in) :: dt
  double precision, dimension(dim), intent(out) :: new

  new = old + dt*fp

  !SECOND STEP : IMPLICIT TRAPEZOIDAL RULE with full time step

end subroutine psecond_step

  !---------------------------------------------------------------------------------------------------------------------!



end module euler_modified
