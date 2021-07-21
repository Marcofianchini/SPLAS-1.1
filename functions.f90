  module functions


  implicit none

  private
  save

  integer, parameter :: dim=70
  integer, parameter :: dn=3

  public :: linspace_dp,logspace_dp,grazing,calcul_phi,fiperp
  public :: allometric,monod,p_growth
  public :: nut_metrics,graz_metrics,time_mean_biomass
  public :: shannon_evenness,initialise,eff_grazing,small_biomass
  public :: supply_per,supply_season, supply_day
  public :: mu_tang
  public :: phyto_sink
  public :: phyto_carbon,phyto_c_export,zoo_c_export
  public :: det_c_export
  public :: irradiance,p_light_lim,ss_irradiance

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-------------------------------------------------------------------------------------------------------------------------!

subroutine linspace_dp(xmin,xmax,x)
    implicit none
    double precision,intent(in) :: xmin,xmax
    double precision,intent(out) :: x(:)
    integer :: i,n
    n = size(x)
    if (n == 1) then
       if(xmin /= xmax) then
          write(0,'("ERROR: Cannot call linspace with n=1 and xmin /= xmax")')
          stop
       else
          x = xmin
       end if
    else
       do i=1,n
          x(i) = (xmax-xmin) * real(i-1) / real(n-1) + xmin
       end do
    end if
  end subroutine linspace_dp

!-------------------------------------------------------------------------------------------------------------------------!

  subroutine logspace_dp(xmin,xmax,x)
    implicit none
    double precision,intent(in) :: xmin,xmax
    double precision,intent(out) :: x(:)
    if (size(x) == 1 .and. xmin /= xmax) then
       write(0,'("ERROR: Cannot call logspace with n=1 and xmin /= xmax")')
       stop
    end if
    call linspace_dp(log10(xmin),log10(xmax),x)
    x = 10.**x
  end subroutine logspace_dp

!-------------------------------------------------------------------------------------------------------------------------!

subroutine grazing(p,z,phi,ir,k,totp,gr) !zooplankton grazing on phytoplankton, Holling Type II function
implicit none
double precision, intent(in) :: p(dim),z(dim),phi(dim,dim),ir(dim),k,totp(dim)
double precision, intent(out) :: gr(:,:)
integer :: i,j
double precision :: temp(dim),asd(dim)

temp = 1.0/(k + totp)      !totp calculated as call fiperp(phi,p,totp)

asd = ir*z*temp

  do j = 1,dim
    do i = 1,dim

      gr(i,j) = asd(j)*phi(i,j)*p(i)

     end do
   end do

end subroutine grazing

!-------------------------------------------------------------------------------------------------------------------------!

subroutine calcul_phi(x,xo,dx,phi) !calculate relative preference of Zj for Pi
implicit none
double precision, intent(in) :: x(dim),xo(dim),dx
double precision, intent(out) :: phi(:,:)
integer :: i,j

   do j = 1,dim
     do i = 1,dim

      phi(i,j) = exp(-(((log10(x(i)) - log10(xo(j)))/dx)**2))

      end do
    end do

end subroutine calcul_phi

!-------------------------------------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------------!

subroutine fiperp(a,c,b)   !matrix multiplication for a vector

implicit none

#ifdef MKLSPLAS
include 'mkl.fi'
#else
external :: dgemv
#endif

double precision, intent(in) :: a(dim,dim),c(dim)
double precision, intent(out) :: b(dim)
!!!!!!INTERNAL PARAMETERS
double precision, parameter :: alpha = 1.0
double precision, parameter :: beta = 1.0
integer, parameter :: incrx=1
integer, parameter :: incry=1

  !!!!!!!!!!!!!!!!!!! CALLING BLAS SUBROUTINE THAT PERFOMES MATRIX X VECTOR
  !!!!!!! WE WANT TO DO b = a*c + b with a matrix, b and c vectors
  !( alpha = 1.0 and beta = 1.0 )

  b=0.
  call dgemv('n',dim,dim,alpha,a,dim,c,incrx,beta,b,incry)

! or one can use fortran intrinsic command for matrix multiplication
  !b = matmul(a,c)

end subroutine fiperp

!-------------------------------------------------------------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-------------------------------------------------------------------------------------------------------------------------!

subroutine allometric(xp,xz,xo,pvol,vmax,pquota,mu,hs,ir,s,ss)  !calculates allometric relationships
implicit none
double precision, intent(in) :: xp(dim),xz(dim),s  !arrays of phyto and zoo sizes + P sinking
double precision, intent(out) :: xo(dim),vmax(dim),pquota(dim),mu(dim),hs(dim),ir(dim),ss(dim) !parameters
double precision, intent(inout) :: pvol(dim)
double precision :: scaling,sphere
integer :: i

scaling = 1.00000000E-12
sphere = 0.52359877559

                            ! phytoplankton growth rate - Tang (1995)
 call mu_tang(xp,mu)
 !mu = 3.49*(pvol**(-0.15))

!non isometric model by Wirtz, 2012 that accounts for non-size related features of predators
xo = 0.16*xz*exp(-0.02*log(xz)*log(xz))

!!!!!!!!!!!!!!!!!!!!!!!!!

pvol = sphere*(xp**3) ! Volume of phytoplankton cell calculated from ESD (in micro-meters^3)

vmax = 9.10*scaling*(pvol**0.67)    ! maximum uptake rate of phytoplankton cell - Litchman et al., 2007
pquota = 1.36*scaling*(pvol**0.77)  ! minimum internal quota of phytoplankton cell - Litchman et al., 2007

hs = 0.1*xp                 ! phytoplankton half saturation for nutrients - Eppley et al., 1969
!hs = 0.17*(10**(-3))*(pvol**0.27) ! half saturation by Litchman et al., 2007

ir = 26.0*(xz**(-0.4))      ! used by Banas; zooplankton ingestion rate - Hansen et al., 1997

call phyto_sink(s,xp,ss)    ! phytoplankton size-dependent sinking, Smayda (1970)

end subroutine allometric

!-------------------------------------------------------------------------------------------------------------------------!

subroutine monod(n,hs,up) !phytoplankton uptake of nutrients (Monod)
implicit none
double precision, intent(in) :: n,hs(dim)
double precision, intent(out) :: up(dim)
double precision :: temp(dim)
double precision :: x

x = max(n,0.0)

  temp = 1.0/(x+hs)

  up = x*temp

end subroutine monod

!-------------------------------------------------------------------------------------------------------------------------!

subroutine p_growth(mu,qmin,q,g) !phytoplankton growth according to quota model (Droop)
implicit none
double precision, intent(in) :: mu(dim),qmin(dim),q(dim)
double precision, intent(out) :: g(dim)
double precision :: temp(dim)
integer :: i

do i = 1,dim
  temp(i) = 1.0/(max(q(i),qmin(i)))
end do

  g = mu*(1.0 - qmin*temp)

end subroutine p_growth

!-------------------------------------------------------------------------------------------------------------------------!

subroutine nut_metrics(hs,n,nl)
implicit none
double precision, intent(in) :: n,hs(dim)
double precision, intent(out) :: nl(dim)
double precision :: temp(dim)
double precision :: nut

nut = max(n,0.0)

temp = 1.0/(hs+nut)

nl = hs*temp

end subroutine nut_metrics

!-------------------------------------------------------------------------------------------------------------------------!

subroutine graz_metrics(mu,p,gr,grl)    !may have some numerical problems if p is very small
implicit none
double precision, intent(in) :: mu(dim),p(dim),gr(dim,dim)
double precision, intent(out) :: grl(dim)
double precision :: tmp(dim)
integer :: j

tmp = 1.0/(mu*p)

grl = tmp*sum(gr,2)

end subroutine graz_metrics

!-------------------------------------------------------------------------------------------------------------------------!

subroutine eff_grazing(mu,z,phi,ir,k,tot,gm) !zooplankton EFFECTIVE grazing on phytoplankton (per unit mass of phytoplankton)*1/mu
implicit none
double precision, intent(in) :: mu(dim),z(dim),phi(dim,dim),ir(dim),k,tot(dim)
double precision :: gr(dim,dim)
double precision, intent(out) :: gm(dim)
integer :: i,j
double precision :: temp(dim),umu(dim),asd(dim)

temp = 1.0/(k + tot)
umu = 1.0/mu
asd = ir*z*temp

  do j = 1,dim
    do i = 1,dim

      gr(i,j) = asd(j)*phi(i,j)

     end do
   end do

   gm = umu*sum(gr,2)

end subroutine eff_grazing

!-------------------------------------------------------------------------------------------------------------------------!

subroutine time_mean_biomass(plkt,period,dt,ave)
implicit none

double precision, intent(in) :: plkt(:,:)
double precision, intent(in) :: period,dt
double precision, intent(out) :: ave(dim)
double precision :: freq

freq = 1.0/period

ave = freq*sum(dt*plkt,2)

end subroutine time_mean_biomass

!-------------------------------------------------------------------------------------------------------------------------!

subroutine shannon_evenness(pl,num,se)
implicit none

#ifdef MKLSPLAS
include 'mkl.fi'
#else
double precision :: dasum
#endif

double precision, intent(in) :: pl(dim),num
double precision, intent(out) :: se
double precision :: rec,norm,den(dim)
integer :: i

rec = 1.0/log(num)                  !should have a (-) sign but we use dasum (sums abs values)
norm = 1.0/(dasum(dim,pl,1))

!den = log(pl*norm)                  !is negative

do i = 1,dim

   if(pl(i)*norm == 0.0) then
      den(i) = 0.0
   else
      den(i) = log(pl(i)*norm)
   end if

end do

!print*,rec*dasum(dim,norm*pl*den,1), -rec*sum(norm*pl*den) !rec should be negative (-), but we use dasum
                                                            !so se is already positive

se = rec*dasum(dim,(pl*norm)*den,1)

end subroutine shannon_evenness

!-------------------------------------------------------------------------------------------------------------------------!

subroutine small_biomass(plkt,xp,sma)
implicit none

double precision, intent(in) :: plkt(dim),xp(dim)
double precision, intent(out) :: sma

sma = sum(plkt,mask=xp<2.0)

end subroutine small_biomass

!-------------------------------------------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------------------------------------------------!

subroutine mu_tang(x,mu) !valid in the range 1-100 micro-meters
implicit none

double precision,intent(in) :: x(dim)
double precision,intent(out) :: mu(dim)

mu = 2.6*(x**(-0.45))

end subroutine mu_tang

!-------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------------!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-------------------------------------------------------------------------------------------------------------------------!

subroutine initialise(p,z,n,quota,ni,pi,zi,iquota) !initialises everything
implicit none
double precision, intent(inout) :: p(dim),z(dim),n(dn),quota(dim)
double precision, intent(in) :: ni,pi,zi,iquota
!integer :: i

!no cycle for now..
!do i =1,dn

!N
 n(1) = ni
!P
 n(2) = ni/16.0
!C
 n(3) = (106.0*ni)/16.0

!end do

 p = pi
 z = zi

!some non homogeneous initial conditions for P and Z (winners and complementary to winners for S = 5)
  !p = 0.0
  !z = 0.0
  !p(1) = 0.0 !pi
  !p(11) = 0.0 !pi
  !p(22) = 0.0 !pi
  !p(30) = 0.0 !pi
  !p(31) = 0.0 !pi
  !p(32) = 0.0 !pi
  !p(41) = 0.0 !pi
  !p(42) = 0.0 !pi
  !p(43) = 0.0 !pi
  !p(52) = 0.0 !pi
  !p(57) = 0.0 !pi
  !p(58) = 0.0 !pi
  !p(59) = 0.0 !pi
  !p(60) = 0.0 !pi
  !z(15) = 0.0 !zi
  !z(25) = 0.0 !zi
  !z(26) = 0.0 !zi
  !z(27) = 0.0 !zi
  !z(35) = 0.0 !zi
  !z(36) = 0.0 !zi
  !z(37) = 0.0 !zi
  !z(44) = 0.0 !zi
  !z(45) = 0.0 !zi
  !z(46) = 0.0 !zi
  !z(58) = 0.0 !zi
  !z(59) = 0.0 !zi
  !z(60) = 0.0 !zi
  !z(69) = 0.0 !zi

 !low cellular quota
 quota = iquota*1.0E-11
 !very high cellular quota
 !quota = iquota
!leave a very low quota, for stability of export...order 10^-11

end subroutine initialise

!-------------------------------------------------------------------------------------------------------------------------!

double precision function supply_per(dt,i) result(r) !nutrient supply..periodic (time-varying)
implicit none
double precision, intent(in) :: dt
integer, intent(in) :: i
double precision, parameter :: pie = 3.1415926535
double precision, parameter :: OM = pie/365.0
double precision, parameter :: A = 1.0

r = A*(sin(OM*i*dt))**2    

end function supply_per

double precision function supply_season(dt,i) result(r) !Light supply..seasonal(time-varying) Fianchini (2021)
implicit none
double precision, intent(in) :: dt
integer, intent(in) :: i
double precision, parameter :: pie = 3.1415926535
double precision, parameter :: OM = pie/(365.0/2)
double precision, parameter :: A = 1.0

r = A*(sin(OM*i*dt))**2

end function supply_season

double precision function supply_day(dt,i) result(r) !Light supply..daily (time-varying)  Fianchini (2021)
implicit none									     ! BE CAREFUL! not implemented in solver.f90
double precision, intent(in) :: dt
integer, intent(in) :: i
double precision, parameter :: pie = 3.1415926535
double precision, parameter :: OM = pie
double precision, parameter :: A = 1.0

r = A*(sin(OM*i*dt))**2

end function supply_day

!-------------------------------------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------------!

subroutine phyto_sink(a,xp,sink) !Sinking of phytoplankton (size-dependent)
implicit none
double precision,intent(in) :: a,xp(dim)
double precision,intent(out) :: sink(dim)
double precision :: c    !exponent -> Stoke's law suggests 2
                         !         -> Smayda, 1970 poses 1.17
c = 1.17

sink = a*(xp**c)

end subroutine phyto_sink

!-------------------------------------------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------------!
!!!!!!!!!!!!CARBON EXPORT!!!!!!!!!!!!!!!!!CARBON EXPORT!!!!!!!!!!!!!!!!!!!!!!!!CARBON EXPORT!!!!!!!!!!!!!!!!!!!!!!!!!!CARBON EXPORT!!!!!!!!!!!!!!!!!!!!!CARBON EXPORT!!!!!!!!!!!!!!
!-------------------------------------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------------------------------------------------!

subroutine phyto_carbon(vol,carbon)
implicit none

double precision,intent(in) :: vol(dim)
double precision,intent(out) :: carbon(dim)

carbon = 0.0109*((vol)**0.991) !Montagnes et al., 1994


end subroutine phyto_carbon

!-------------------------------------------------------------------------------------------------------------------------!
!NOW MUST CONSIDER DETRITUS DYNAMICS TO CALCULATE CARBON EXPORT
!-------------------------------------------------------------------------------------------------------------------------!

subroutine phyto_c_export(sink,ncells,c,pexpc,totpexpc)
implicit none

double precision,intent(in) :: sink(dim),ncells(dim),c(dim)
double precision,intent(out) :: pexpc(dim),totpexpc
double precision :: scaling

scaling = 1.00000000E-09

pexpc = sink*ncells*c*scaling   !mgC m-2 day-1

totpexpc = sum(pexpc)

!uncomment lines below only for debug
!print*,'CELLS ',sum(ncells)
!print*,'P EXPORT ',totpexpc

end subroutine phyto_c_export

!-------------------------------------------------------------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!zoo_c_export() is not used in 1-D version, because detritus dynamics is included
!and zooplankton egestion (export) is a term in the detritus equation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zoo_c_export(ege,graz1,zexpc)
implicit none

double precision,intent(in) :: ege,graz1(:,:)
double precision,intent(out) :: zexpc
double precision :: scaling
double precision :: a
double precision :: fraction

!a = 0.01
a = 6.625  !Redfield ratio C:N = 106:16
!a = 1.0

!units = milligrams
scaling = 1.00000000

!fraction of feg, if there is another egestion term in dN/dt then fraction < 1 (see solver.f90)
fraction = 1.0

zexpc = fraction*ege*a*sum(graz1)*14.0*scaling  !mgC m-2 day-1

end subroutine zoo_c_export
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-------------------------------------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-------------------------------------------------------------------------------------------------------------------------!

subroutine det_c_export(w,det,dexpc)
implicit none

double precision,intent(in) :: w,det
double precision,intent(out) :: dexpc
double precision :: a,b

a = 6.625  !Redfield ratio C:N = 106:16
b = 14.0   !Molar weight of Nitrogen

dexpc = w*det*a*b !mgC m-2 day-1

end subroutine det_c_export

!-------------------------------------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-------------------------------------------------------------------------------------------------------------------------!
!IRRADIANCE AND P LIGHT LIMITATION TO GROWTH RATE
!-------------------------------------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-------------------------------------------------------------------------------------------------------------------------!

double precision function irradiance(k,dz,i0) result(irr)
implicit none

integer,intent(in) :: k
double precision,intent(in) :: dz
double precision :: z
double precision,intent(in) :: i0
double precision,parameter :: kd = 0.062

!default = float lovbio015c during 18.03.2015 (late winter/spring), central western Mediterranean
!A challenge will be to use SPLAS with more argo floats, in different periods of the year

z = (float(k)-0.5)*dz

irr = i0*exp(-kd*z)

end function irradiance

!-------------------------------------------------------------------------------------------------------------------------!

double precision function ss_irradiance(k,dz,p,i0) result(irr)
implicit none

integer,intent(in) :: k
double precision,intent(in) :: dz
double precision,intent(in) :: p(:)
double precision :: z
double precision,intent(in) :: i0
double precision,parameter :: kw = 0.03   !(Oguz et al., 2013)
double precision,parameter :: kp = 0.01   !(Oguz et al., 2013)
integer :: i
double precision :: shadow

!default = float lovbio015c during 18.03.2015 (late winter/spring), central western Mediterranean
!A challenge will be to use SPLAS with more argo floats, in different periods of the year

z = (float(k)-0.5)*dz

!introduce self shading of phytoplankton
shadow = 0.0

if (k == 2) then
  irr = i0*exp(-kw*z)
else
  do i=2,k
    shadow = shadow + p(i)
  end do
  irr = i0*exp(-kw*z - kp*shadow*dz)
end if

end function ss_irradiance

!-------------------------------------------------------------------------------------------------------------------------!

subroutine p_light_lim(irr,mu0,ll)
implicit none

double precision,intent(in) :: irr,mu0(dim)
double precision,intent(out) :: ll(:)
double precision :: invmu(dim)
double precision :: th,alfa
!
integer :: j

th = 0.02        !Chl a to C ratio.... A first approximation for the Chl:C range may be derived from Figure 1 in Cloern et al. (1995)
                 !to be order 0.06 (mg Chl a) / (mg C)...for details see file in the home of my MAC /Desktop/SPLAS
!alfa = 3.83E-07 initial slope of P-I curve, but this is for energy in micro Einstein for m2, found in Ward & Follows, 2016
!EQF(micro Einstein) = E(W/m2)*wavelength*0.836E-02
alfa = 0.32      !this is for energy in W/m2, from Oguz et al., 2013

invmu = 1.0/mu0

do j = 1,dim
  ll(j) = 1-exp(-irr*th*alfa*invmu(j))
end do

end subroutine p_light_lim

!-------------------------------------------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end module functions
