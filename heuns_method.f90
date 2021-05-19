module heuns_method


  implicit none
  public
  save

  integer, parameter :: dim=70
  integer, parameter :: dn=3

  public :: ppredictor,pcorrector,npredictor,ncorrector
  public :: detpredictor,detcorrector

contains

!---------------------------------------------------------------------------------------------------------------------!

subroutine ppredictor(pred,old,dt,fo)

implicit none

double precision, dimension(dim), intent(in) :: old
double precision, dimension(dim), intent(in) :: fo
double precision, intent(in) :: dt
double precision, dimension(dim), intent(out) :: pred

!PREDICTOR : EXPLICIT EULER
pred = old + dt*fo

end subroutine ppredictor

!---------------------------------------------------------------------------------------------------------------------!


subroutine pcorrector(pred,old,dt,fp,new)

implicit none

double precision, dimension(dim), intent(in) :: pred,old
double precision, dimension(dim), intent(in) :: fp
double precision, intent(in) :: dt
double precision, dimension(dim), intent(out) :: new

!CORRECTOR : IMPLICIT TRAPEZOIDAL RULE
new = old + 0.5*(pred - old) + 0.5*dt*fp

end subroutine pcorrector

!---------------------------------------------------------------------------------------------------------------------!

!---------------------------------------------------------------------------------------------------------------------!
!---------------------------------------------------------------------------------------------------------------------!
!TO AVOID ALLOCATION PROBLEMS 4 SUBROUTINES ARE DEFINED, THAT TAKE AS INPUT DIFFERENT ARRAYS (plankton above, and nutrients below)
!---------------------------------------------------------------------------------------------------------------------!
!---------------------------------------------------------------------------------------------------------------------!

!---------------------------------------------------------------------------------------------------------------------!

subroutine npredictor(pred,old,dt,fo)

implicit none

double precision, dimension(dn), intent(in) :: old
double precision, dimension(dn), intent(in) :: fo
double precision, intent(in) :: dt
double precision, dimension(dn), intent(out) :: pred

!PREDICTOR : EXPLICIT EULER
pred = old + dt*fo

end subroutine npredictor

!---------------------------------------------------------------------------------------------------------------------!

subroutine ncorrector(pred,old,dt,fp,new)

implicit none

double precision, dimension(dn), intent(in) :: pred,old
double precision, dimension(dn), intent(in) :: fp
double precision, intent(in) :: dt
double precision, dimension(dn), intent(out) :: new

!CORRECTOR : IMPLICIT TRAPEZOIDAL RULE
new = old + 0.5*(pred - old) + 0.5*dt*fp

end subroutine ncorrector

!---------------------------------------------------------------------------------------------------------------------!

!---------------------------------------------------------------------------------------------------------------------!
!---------------------------------------------------------------------------------------------------------------------!
!---------------------------------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! NEW DETRITUS COMPARTMENT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------------------------------------------------------!
!---------------------------------------------------------------------------------------------------------------------!

!---------------------------------------------------------------------------------------------------------------------!

subroutine detpredictor(pred,old,dt,fo)

implicit none

double precision, intent(in) :: old
double precision, intent(in) :: fo
double precision, intent(in) :: dt
double precision, intent(out) :: pred

!PREDICTOR : EXPLICIT EULER
pred = old + dt*fo

end subroutine detpredictor

!---------------------------------------------------------------------------------------------------------------------!

subroutine detcorrector(pred,old,dt,fp,new)

implicit none

double precision, intent(in) :: pred,old
double precision, intent(in) :: fp
double precision, intent(in) :: dt
double precision, intent(out) :: new

!CORRECTOR : IMPLICIT TRAPEZOIDAL RULE
new = old + 0.5*(pred - old) + 0.5*dt*fp

end subroutine detcorrector

!---------------------------------------------------------------------------------------------------------------------!

!---------------------------------------------------------------------------------------------------------------------!
!---------------------------------------------------------------------------------------------------------------------!
!---------------------------------------------------------------------------------------------------------------------!

end module heuns_method
