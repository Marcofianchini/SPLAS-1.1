module solver

use functions
use heuns_method

implicit none

  public :: update_time_varying_stuff,eq_solver

contains

!-------------------------------------------------------------------------------------------------------------------------!

subroutine update_time_varying_stuff(nut,phs,phim,phyto,zoo,zir,kz,up,totpredation,grazij,mu0,qmin,q,growth)

implicit none
double precision, dimension(dim,dim), intent(in) :: phim
double precision, dimension(dn),intent(in) :: nut
double precision, dimension(dim), intent(in) :: phyto,zoo
double precision, dimension(dim), intent(in) :: phs,zir
double precision, dimension(dim), intent(in) :: mu0,qmin,q
double precision, intent(in) :: kz
double precision, dimension(dim), intent(out) :: up,growth
double precision, dimension(dim), intent(inout) :: totpredation
double precision, dimension(dim,dim), intent(out) :: grazij

call monod(nut(1),phs,up)  !phytoplankton uptake function (without uptake rate inside)

call p_growth(mu0,qmin,q,growth)   !phytoplankton growth function (growth rate calculated with Droop's internal quota model)

call fiperp(phim,phyto,totpredation)   !total predation zoo on phyto
call grazing(phyto,zoo,phim,zir,kz,totpredation,grazij) !zoo on phyto

end subroutine update_time_varying_stuff

!-------------------------------------------------------------------------------------------------------------------------!

subroutine eq_solver(tv,n,p,z,uptake,tot,graz,phim,phs,zir,kz,supply,mu_0,eps,feg,pm,zhp,&
                     delta_t,nt,pstep,nmet,gmet,pmean,zmean,psi,zsi,totalp,totalz,smallp,xp,&
                     totalmeanphytomass,totalmeanzoomass,pmeansi,zmeansi,&
                     vmax,qmin,psizesink,pq,pgrowth,&
                     p_c,p_c_export,z_c_export,p_c_total_export,&
                     D,dz,m,detritus,irradiance,ivalue)

implicit none

#ifdef MKLSPLAS
include 'mkl.fi'
#else
double precision :: dasum
external :: dcopy
#endif

!common variables,parameter and input/output
logical,intent(in) :: tv
double precision, dimension(dim,dim), intent(in) :: phim
double precision, dimension(dim), intent(in) :: phs,zir,mu_0,xp,vmax,qmin,psizesink
double precision, intent(in) :: kz,supply,eps,feg,pm,zhp
double precision, intent(in) :: delta_t  !time-step
integer, intent(in) :: m  !number of vertical layers
double precision, intent(in) :: dz !step of the space grid
double precision, dimension(m), intent(in) :: D !diffusion coefficient (dimension m, allowed to vary with depth)
integer, intent(in) :: nt !number of time steps in the simulation
integer, intent(in) :: pstep

double precision, dimension(dn,int(nt/pstep)+1,m),intent(inout) :: n
double precision, dimension(dim,int(nt/pstep)+1,m), intent(inout) :: p,z,pq
double precision, dimension(dim,m), intent(inout) :: uptake,pgrowth,tot
double precision, dimension(dim,dim,m), intent(inout) :: graz
!some important output quantities
double precision, dimension(dim,int(nt/pstep)+1,m), intent(out) :: nmet,gmet
double precision, dimension(dim,m), intent(out) :: pmean,zmean
double precision, dimension(int(nt/pstep)+1,m), intent(out) :: psi,zsi,totalp,totalz
double precision, dimension(m), intent(out) :: smallp,totalmeanphytomass,totalmeanzoomass
double precision, dimension(m), intent(out) :: pmeansi,zmeansi

!internal variables, private (not known outside the subroutine)
double precision, dimension(dim,m) :: pb,zb,pqb !buffer quantities calculated at every step!
double precision, dimension(dn,m) :: nb     !buffer quantities calculated at every step!
double precision, dimension(dn,m) :: nold,n_pred
double precision, dimension(dim,m) :: pold,zold,uptake_old,tot_old
double precision, dimension(dim,dim,m) :: graz_old
double precision, dimension(dim,m) :: p_pred,z_pred,uptake_pred,tot_pred
double precision, dimension(dim,dim,m) :: graz_pred
double precision, dimension(dn,m) :: nfold,nfpred
double precision, dimension(dim,m) :: pfold,zfold,pfpred,zfpred
double precision :: n_species = real(dim)
double precision :: supply_mod
double precision, dimension(dn,m) :: ntemp
double precision, dimension(dim,m) :: ptemp,ztemp
double precision, dimension(dim,m) :: pq_old,pq_pred,pq_fold,pq_fpred,pq_temp
double precision, dimension(dim,m) :: pgrowth_old,pgrowth_pred

!carbon variables (export is only from the last layer)!!
double precision, dimension(dim,int(nt/pstep)+1), intent(inout) :: p_c_export
double precision, dimension(dim), intent(inout) :: p_c
double precision, dimension(int(nt/pstep)+1), intent(inout) :: z_c_export,p_c_total_export

!DETRITUS variables
double precision, dimension(m) :: detfold,detfpred,det_old,det_pred,detb,det_temp
double precision, dimension(int(nt/pstep)+1,m), intent(inout) :: detritus

!LIGHT VARIABLES
double precision, dimension(dim) :: pll
double precision :: irr
double precision, dimension(m) :: irr_temp,irrb
double precision, dimension(int(nt/pstep)+1,m), intent(inout) :: irradiance
double precision, dimension(nt+1), intent(in) :: ivalue
double precision :: i0
!
integer :: k,j,i
double precision :: A = 0.10

!detritus parameters (per il momento non ho voglia/tempo di metterli come parametri del json, poi vedremo)
double precision :: remin = 0.40
double precision :: wdet = 0.0
double precision :: gamma = 0.8

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!I.C.
!print*,det_old
do i = 1,m

  call dcopy(dn,n(:,1,i),1,nb(:,i),1)
  call dcopy(dim,p(:,1,i),1,pb(:,i),1)
  call dcopy(dim,z(:,1,i),1,zb(:,i),1)
  call dcopy(dim,pq(:,1,i),1,pqb(:,i),1)
  !detritus
  detb(i) = detritus(1,i)
  !irradiance
  irrb(i) = ivalue(1)

  call update_time_varying_stuff(nb(:,i),phs,phim,pb(:,i),zb(:,i),zir,kz,uptake(:,i),tot(:,i),graz(:,:,i),mu_0,qmin,pqb(:,i),pgrowth(:,i))

  call nut_metrics(phs,n(1,1,i),nmet(:,1,i))
  call graz_metrics(mu_0,p(:,1,i),graz,gmet(:,1,i))

  call shannon_evenness(p(:,1,i),n_species,psi(1,i))
  call shannon_evenness(z(:,1,i),n_species,zsi(1,i))

  totalp(1,i) = dasum(dim,pb(:,i),1)
  totalz(1,i) = dasum(dim,zb(:,i),1)

  call dcopy(dn,nb(:,i),1,ntemp(:,i),1)
  call dcopy(dim,pb(:,i),1,ptemp(:,i),1)
  call dcopy(dim,zb(:,i),1,ztemp(:,i),1)
  call dcopy(dim,pqb(:,i),1,pq_temp(:,i),1)
  !detritus
  det_temp(i) = detb(i)
  !irradiance
  irr_temp(i) = irrb(i)
  
end do

!!!!!!!! must recalculate export based on detritus
call phyto_c_export(psizesink,p(:,1,m)/pq(:,1,m),p_c,p_c_export(:,1),p_c_total_export(1))
call det_c_export(wdet,detritus(1,m),z_c_export(1))
!!!!!!!!

!B.C. already in (all equal initial distribution)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!TIME LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do k =1,nt
i0 = ivalue(k+1)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!UPDATE OF THE VARIABLES...OUT OF SPACE LOOP!!!!!!!!!!
  do i = 1,m

    call update_time_varying_stuff(nb(:,i),phs,phim,pb(:,i),zb(:,i),zir,kz,uptake(:,i),tot(:,i),graz(:,:,i),mu_0,qmin,pqb(:,i),pgrowth(:,i))

    !UPDATE OF THE VARIABLES
    call dcopy(dn,nb(:,i),1,nold(:,i),1)
    call dcopy(dim,pb(:,i),1,pold(:,i),1)
    call dcopy(dim,zb(:,i),1,zold(:,i),1)
    call dcopy(dim,pqb(:,i),1,pq_old(:,i),1)

    call dcopy(dim,uptake(:,i),1,uptake_old(:,i),1)
    call dcopy(dim,pgrowth(:,i),1,pgrowth_old(:,i),1)
    call dcopy(dim,tot(:,i),1,tot_old(:,i),1)
    call dcopy(dim,det_temp(i),1,det_old(i),1)
    do j = 1,dim
      call dcopy(dim,graz(:,j,i),1,graz_old(:,j,i),1)
    end do
  !print*,"det_old_pre",det_old
  end do
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!UPDATE OF THE VARIABLES...OUT OF SPACE LOOP!!!!!!!!!!

  if(tv) then
    supply_mod = supply_per(delta_t,k)
  else
    supply_mod = 1.0
  end if

  !!!!!!!!!!!!!SPACE LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do i = 2,m-1
    
    !LIGHT LIMITATION, DEPENDING ON DEPTH
    !call p_light_lim(irradiance(i,dz),mu_0,pll)
    !
    !LIGHT LIMITATION, DEPENDING ON DEPTH (with self shading)
    irr_temp(i) = ss_irradiance(i,dz,sum(pold,1),i0)
    !print*,"irradiance",i,irr_temp(i)
    call p_light_lim(irr_temp(i),mu_0,pll)

    !INTEGRATION OF THE EQUATIONS
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    nfold(:,i) = ( i*A*supply*supply_mod - dasum(dim,vmax(:)*uptake_old(:,i)*(pold(:,i)/pq_old(:,i)),1) &
                  + (1-eps-feg)*sum(graz_old(:,:,i)) + remin*det_old(i) & !+ (pm/5.0)*dasum(dim,pold(:,i),1) + (feg/10.0)*sum(graz_old(:,:,i)) &
                  + (D(i)/(dz*dz))*(nold(:,i+1) - 2*nold(:,i) + nold(:,i-1)) )
    pfold(:,i) = ( pgrowth_old(:,i)*pll(:)*pold(:,i) - sum(graz_old(:,:,i),2) - pm*pold(:,i) - (psizesink/dz)*(pold(:,i) - pold(:,i-1)) &
                  + (D(i)/(dz*dz))*(pold(:,i+1) - 2*pold(:,i) + pold(:,i-1)) )
    zfold(:,i) = ( eps*sum(graz_old(:,:,i),1) - zhp*zold(:,i)*dasum(dim,zold(:,i),1) &
                  + (D(i)/(dz*dz))*(zold(:,i+1) - 2*zold(:,i) + zold(:,i-1)) )
    pq_fold(:,i) = ( vmax*uptake_old(:,i) - pgrowth_old(:,i)*pq_old(:,i) )
    !new compartment: detritus
    detfold(i) = ( -remin*det_old(i) + pm*dasum(dim,pold(:,i),1) + gamma*zhp*dasum(dim,zold(:,i),1)*dasum(dim,zold(:,i),1) &
                   + feg*sum(graz_old(:,:,i)) - (wdet/dz)*(det_old(i) - det_old(i-1)) &
                   + (D(i)/(dz*dz))*(det_old(i+1) - 2*det_old(i) + det_old(i-1)) )

    call npredictor(n_pred(:,i),nold(:,i),delta_t,nfold(:,i))
    call ppredictor(p_pred(:,i),pold(:,i),delta_t,pfold(:,i))
    call ppredictor(z_pred(:,i),zold(:,i),delta_t,zfold(:,i))
    call ppredictor(pq_pred(:,i),pq_old(:,i),delta_t,pq_fold(:,i))
    call detpredictor(det_pred(i),det_old(i),delta_t,detfold(i))

    call update_time_varying_stuff(n_pred(:,i),phs,phim,p_pred(:,i),z_pred(:,i),zir,kz,&
                                   uptake_pred(:,i),tot_pred(:,i),graz_pred(:,:,i),&
                                   mu_0,qmin,pq_pred(:,i),pgrowth_pred(:,i))

    !LIGHT LIMITATION, DEPENDING ON DEPTH (with self shading)
  
    irr_temp(i) = ss_irradiance(i,dz,sum(p_pred,1),i0)  
    !print*,"irradiance_post",i,irr_temp(i)
    call p_light_lim(irr_temp(i),mu_0,pll)

    nfpred(:,i) = ( i*A*supply*supply_mod - dasum(dim,vmax(:)*uptake_pred(:,i)*(p_pred(:,i)/pq_pred(:,i)),1) &
                   + (1-eps-feg)*sum(graz_pred(:,:,i)) + remin*det_pred(i) & !+ (pm/5.0)*dasum(dim,p_pred(:,i),1) + (feg/10.0)*sum(graz_pred(:,:,i)) &
                   + (D(i)/(dz*dz))*(n_pred(:,i+1) - 2*n_pred(:,i) + n_pred(:,i-1)) )
    pfpred(:,i) = ( pgrowth_pred(:,i)*pll(:)*p_pred(:,i) - sum(graz_pred(:,:,i),2) - pm*p_pred(:,i) - (psizesink/dz)*(p_pred(:,i) - p_pred(:,i-1)) &
                   + (D(i)/(dz*dz))*(p_pred(:,i+1) - 2*p_pred(:,i) + p_pred(:,i-1)) )
    zfpred(:,i) = ( eps*sum(graz_pred(:,:,i),1) - zhp*z_pred(:,i)*dasum(dim,z_pred(:,i),1) &
                   + (D(i)/(dz*dz))*(z_pred(:,i+1) - 2*z_pred(:,i) + z_pred(:,i-1)) )
    pq_fpred(:,i) = ( vmax*uptake_pred(:,i) - pgrowth_pred(:,i)*pq_pred(:,i) )
    !new compartment: detritus
    detfpred(i) = ( -remin*det_pred(i) + pm*dasum(dim,p_pred(:,i),1) + gamma*zhp*dasum(dim,z_pred(:,i),1)*dasum(dim,z_pred(:,i),1) &
                   + feg*sum(graz_pred(:,:,i)) - (wdet/dz)*(det_pred(i) - det_pred(i-1)) &
                   + (D(i)/(dz*dz))*(det_pred(i+1) - 2*det_pred(i) + det_pred(i-1)) )

    call ncorrector(n_pred(:,i),nold(:,i),delta_t,nfpred(:,i),nb(:,i))
    call pcorrector(p_pred(:,i),pold(:,i),delta_t,pfpred(:,i),pb(:,i))
    call pcorrector(z_pred(:,i),zold(:,i),delta_t,zfpred(:,i),zb(:,i))
    call pcorrector(pq_pred(:,i),pq_old(:,i),delta_t,pq_fpred(:,i),pqb(:,i))
    call detcorrector(det_pred(i),det_old(i),delta_t,detfpred(i),detb(i))

    do j=1,dn
      nb(j,i) = max(nb(j,i),0.0)
    end do

    ntemp(:,i) = ntemp(:,i) + nb(:,i)
    ptemp(:,i) = ptemp(:,i) + pb(:,i)
    ztemp(:,i) = ztemp(:,i) + zb(:,i)
    pq_temp(:,i) = pq_temp(:,i) + pqb(:,i)
    !detritus
    det_temp(i) = det_temp(i) + detb(i)
    det_old(i) = detb(i)  
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  end do
  !!!!!!!!!END OF SPACE LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !BOUNDARY CONDITIONS (Neumann = zero -or prescribed- derivatives at the top and bottom layers)
  !Pay attention to sinking flux from m-1 to m and what should go out from m-th layer
                         !zero flux of anything at the free surface
  nb(:,1) = nb(:,2)
  pb(:,1) = pb(:,2)
  zb(:,1) = zb(:,2)
  pqb(:,1) = pqb(:,2)
  !detritus
  detb(1) = detb(2)
  !irradiance
  irr_temp(1)= i0
                        !fluxes at the bottom layer
  !nb(:,m) = nb(:,m-1)                                   !zero flux of nutrients
  nb(:,m) = nb(:,m-1) + supply*supply_mod*delta_t      !prescribed supply of nutrients
  !nb(:,m) = 1.0                                        !chemostat (Dirichelet boundary, fixed concentration)
  pb(:,m) = pb(:,m-1) + psizesink*(delta_t/dz)*pb(:,m-1)     !no exchange other than sinking
  zb(:,m) = zb(:,m-1)
  pqb(:,m) = pqb(:,m-1) + psizesink*(delta_t/dz)*pqb(:,m-1)  !no exchange other than sinking
  !detritus
  detb(m) = detb(m-1) + wdet*(delta_t/dz)*detb(m-1)        !no exchange other than sinking

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!average calculation for the boundaries

  ntemp(:,1) = ntemp(:,1) + nb(:,1)
  ptemp(:,1) = ptemp(:,1) + pb(:,1)
  ztemp(:,1) = ztemp(:,1) + zb(:,1)
  pq_temp(:,1) = pq_temp(:,1) + pqb(:,1)
  !detritus
  det_temp(1) = det_temp(1) + detb(1)

  ntemp(:,m) = ntemp(:,m) + nb(:,m)
  ptemp(:,m) = ptemp(:,m) + pb(:,m)
  ztemp(:,m) = ztemp(:,m) + zb(:,m)
  pq_temp(:,m) = pq_temp(:,m) + pqb(:,m)
  !detritus
  det_temp(m) = det_temp(m) + detb(m)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!PRINT
  if( mod(k,pstep) == 0 ) then

    do i = 1,m
      call dcopy(dn,ntemp(:,i)/pstep,1,n(:,int(k/pstep)+1,i),1)
      call dcopy(dim,ptemp(:,i)/pstep,1,p(:,int(k/pstep)+1,i),1)
      call dcopy(dim,ztemp(:,i)/pstep,1,z(:,int(k/pstep)+1,i),1)
      call dcopy(dim,pq_temp(:,i)/pstep,1,pq(:,int(k/pstep)+1,i),1)
      detritus(int(k/pstep)+1,i) = det_temp(i)/pstep
      irradiance(int(k/pstep)+1,i) = irr_temp(i)
      !print*,"detritus,detritus(int(k/pstep)+1,10)","d",det_old(10)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call nut_metrics(phs,n(1,int(k/pstep)+1,i),nmet(:,int(k/pstep)+1,i))
      call eff_grazing(mu_0,z(:,int(k/pstep)+1,i),phim,zir,kz,tot_pred(:,i),gmet(:,int(k/pstep)+1,i))
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call shannon_evenness(p(:,int(k/pstep)+1,i),n_species,psi(int(k/pstep)+1,i))
      call shannon_evenness(z(:,int(k/pstep)+1,i),n_species,zsi(int(k/pstep)+1,i))

      totalp(int(k/pstep)+1,i) = dasum(dim,p(:,int(k/pstep)+1,i),1)
      totalz(int(k/pstep)+1,i) = dasum(dim,z(:,int(k/pstep)+1,i),1)

   !reset temporary array to 0, to calculate average on next time printed step
      ntemp(:,i) = 0.0
      ptemp(:,i) = 0.0
      ztemp(:,i) = 0.0
      pq_temp(:,i) = 0.0
      !detritus
      det_temp(i) = 0.0
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!CARBON EXPORT NOW IS ONLY FROM THE LAST LAYER
      call phyto_c_export(psizesink,p(:,int(k/pstep)+1,m)/pq(:,int(k/pstep)+1,m),p_c,p_c_export(:,int(k/pstep)+1),p_c_total_export(int(k/pstep)+1))
      call det_c_export(wdet,detritus(int(k/pstep)+1,m),z_c_export(int(k/pstep)+1))
    !!!!!!now second component of export is calculated from sinking detritus
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!PRINT

end do
!!!!!!!!!END OF TIME LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

do i = 1,m

  call time_mean_biomass(p(:,:,i),nt*delta_t,pstep*delta_t,pmean(:,i))
  call time_mean_biomass(z(:,:,i),nt*delta_t,pstep*delta_t,zmean(:,i))

  totalmeanphytomass(i) = dasum(dim,pmean(:,i),1)
  totalmeanzoomass(i) = dasum(dim,zmean(:,i),1)

  call small_biomass(pmean(:,i),xp,smallp(i))

  pmeansi(i) = (1.0/(nt*delta_t))*dasum(int(nt/pstep)+1,pstep*delta_t*psi(:,i),1)
  zmeansi(i) = (1.0/(nt*delta_t))*dasum(int(nt/pstep)+1,pstep*delta_t*zsi(:,i),1)
  
end do
   print*,"detritus",detritus(10000,10)
end subroutine eq_solver

!-------------------------------------------------------------------------------------------------------------------------!


!-------------------------------------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------------!

end module solver
