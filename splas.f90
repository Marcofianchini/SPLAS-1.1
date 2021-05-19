module splas_class

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  SPLAS module
!
!  Sea PLAnkton Simulator  (SPLAS)
!
!  description
!
!
!  @author : mdepasquale
!  @year : 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use functions
use solver

implicit none

  type splas


    integer :: phyto_class_size,zoo_class_size,number_nutrients
    integer :: nt,print_step
    logical :: splas_time_variation
    double precision :: input_supply_nut,nut_value,initial_p,initial_z,initial_pquota
    double precision :: phyto_x_min,phyto_x_max,phyto_mort,phyto_sinking
    double precision :: zoo_x_min,zoo_x_max
    double precision :: zoo_growth_eff,zoo_eg_frac,zoo_mort,zoo_half_sat
    double precision :: zoo_prey_size_tolerance
    double precision :: sim_length,delta_t
    double precision :: delta_z                                    !SPLAS 1-D VERTICAL SPACING
    double precision :: initial_diff              !initial test value for diffusion coefficient
    integer :: number_layers                                       !SPLAS 1-D NUMBER OF LAYERS
    double precision,allocatable,dimension(:) :: diffusion_coeff   !VERTICAL DIFFUSION COEFFICIENT (allowed to vary among different layers)
    double precision,allocatable,dimension(:,:,:) :: splas_1_phyto_class        !3rd dimension is space (1-D splas)
    double precision,allocatable,dimension(:,:,:) :: splas_1_zoo_class          !3rd dimension is space (1-D splas)
    double precision,allocatable,dimension(:,:,:) :: splas_1_nutrients          !3rd dimension is space (1-D splas)
    double precision,allocatable,dimension(:,:,:) :: splas_1_p_quota            !3rd dimension is space (1-D splas)

    !new output variables needed for calculations done by the module functions.f90
    double precision,allocatable,dimension(:) :: splas_phyto_small_biomass           !one for each layer (1-D splas)
    double precision,allocatable,dimension(:) :: splas_phyto_total_mean_biomass      !one for each layer (1-D splas)
    double precision,allocatable,dimension(:) :: splas_zoo_total_mean_biomass        !one for each layer (1-D splas)
    double precision,allocatable,dimension(:) :: splas_mean_shannon_p,splas_mean_shannon_z   !one for each layer (1-D splas)
    double precision,allocatable,dimension(:) :: xp,xz,xoptz,mu_0,irz,ks,p_size_sink  !!!!!!!!!!!!!!!!!!!!!!!!!!same!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    double precision,allocatable,dimension(:,:) :: nutrient_uptake,p_growth_rate      !2nd dimension is space (1-D splas)
    double precision,allocatable,dimension(:,:) :: total_p_predation_z                !2nd dimension is space (1-D splas)
    double precision,allocatable,dimension(:,:) :: splas_phyto_mean_biomass           !2nd dimension is space (1-D splas)
    double precision,allocatable,dimension(:,:) :: splas_zoo_mean_biomass             !2nd dimension is space (1-D splas)
    double precision,allocatable,dimension(:,:) :: splas_phyto_diversity              !2nd dimension is space (1-D splas)
    double precision,allocatable,dimension(:,:) :: splas_zoo_diversity                !2nd dimension is space (1-D splas)
    double precision,allocatable,dimension(:,:) :: splas_phyto_total_biomass          !2nd dimension is space (1-D splas)
    double precision,allocatable,dimension(:,:) :: splas_zoo_total_biomass            !2nd dimension is space (1-D splas)
    double precision,allocatable,dimension(:,:) :: phizp                              !!!!!!!!!!!!!!!!!!!!!!!!!!same!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    double precision,allocatable,dimension(:,:,:) :: splas_graz_limit                 !3rd dimension is space (1-D splas)
    double precision,allocatable,dimension(:,:,:) :: splas_nut_limit                  !3rd dimension is space (1-D splas)
    double precision,allocatable,dimension(:,:,:) :: z_grazing_p                      !3rd dimension is space (1-D splas)
    double precision,allocatable,dimension(:) :: p_cell_volume,p_v_uptake_rate,p_min_quota

    !carbon variables
    double precision,allocatable,dimension(:,:) :: splas_p_c_export
    double precision,allocatable,dimension(:) :: splas_p_carbon,splas_z_c_export,splas_p_c_tot_export

    !DETRITUS
    double precision,allocatable,dimension(:,:) :: splas_detritus
    
    !IRRADIANCE
    double precision,allocatable,dimension(:,:) :: splas_irradiance
    double precision,allocatable,dimension(:) :: ivalue 
  contains
    procedure :: init
    procedure :: integration
    procedure :: calculation

    procedure :: set_splas_1_nutrients
    procedure :: set_splas_1_phyto_class
    procedure :: set_splas_1_zoo_class
    procedure :: set_splas_1_p_quota
    procedure :: set_splas_1_detritus
    procedure :: set_splas_nut_supply
    procedure :: set_splas_number_nutrients
    procedure :: set_splas_phyto_dim
    procedure :: set_splas_phyto_xmin
    procedure :: set_splas_phyto_xmax
    procedure :: set_splas_phyto_mort
    procedure :: set_splas_phyto_sink
    procedure :: set_splas_zoo_dim
    procedure :: set_splas_zoo_xmin
    procedure :: set_splas_zoo_xmax
    procedure :: set_splas_zoo_growth
    procedure :: set_splas_zoo_eg
    procedure :: set_splas_zoo_mort
    procedure :: set_splas_zoo_hs
    procedure :: set_splas_z_preysize_tolerance
    procedure :: set_splas_nut_value
    procedure :: set_splas_initial_p
    procedure :: set_splas_initial_z
    procedure :: set_splas_initial_pquota
    procedure :: set_splas_sim_length
    procedure :: set_splas_delta_t
    procedure :: set_splas_nt
    procedure :: set_splas_print_step
    procedure :: set_splas_nut_metrics
    procedure :: set_splas_graz_metrics
    procedure :: set_splas_time_mean_mass_p
    procedure :: set_splas_time_mean_mass_z
    procedure :: set_splas_shannon_index_p
    procedure :: set_splas_shannon_index_z
    procedure :: set_splas_total_p_biomass
    procedure :: set_splas_total_z_biomass
    procedure :: set_splas_time_variation

    procedure :: set_splas_phyto_carbon
    procedure :: set_splas_phyto_c_export
    procedure :: set_splas_phyto_c_total_export
    procedure :: set_splas_zoo_c_export

    procedure :: set_splas_delta_z
    procedure :: set_splas_number_layers
    procedure :: set_splas_initial_diff
    procedure :: set_splas_diffusion_coefficient

    procedure :: get_splas_1_nutrients
    procedure :: get_splas_1_phyto_class
    procedure :: get_splas_1_zoo_class
    procedure :: get_splas_1_p_quota
    procedure :: get_splas_1_detritus
    procedure :: get_splas_nut_supply
    procedure :: get_splas_number_nutrients
    procedure :: get_splas_phyto_dim
    procedure :: get_splas_phyto_xmin
    procedure :: get_splas_phyto_xmax
    procedure :: get_splas_phyto_mort
    procedure :: get_splas_phyto_sink
    procedure :: get_splas_zoo_dim
    procedure :: get_splas_zoo_xmin
    procedure :: get_splas_zoo_xmax
    procedure :: get_splas_zoo_growth
    procedure :: get_splas_zoo_eg
    procedure :: get_splas_zoo_mort
    procedure :: get_splas_zoo_hs
    procedure :: get_splas_z_preysize_tolerance
    procedure :: get_splas_sim_length
    procedure :: get_splas_delta_t
    procedure :: get_splas_nt
    procedure :: get_splas_print_step
    procedure :: get_splas_nut_metrics
    procedure :: get_splas_graz_metrics
    procedure :: get_splas_time_mean_mass_p
    procedure :: get_splas_time_mean_mass_z
    procedure :: get_splas_shannon_index_p
    procedure :: get_splas_shannon_index_z
    procedure :: get_splas_total_p_biomass
    procedure :: get_splas_total_z_biomass
    procedure :: get_splas_phyto_small_biomass
    procedure :: get_splas_phyto_total_mean_biomass
    procedure :: get_splas_zoo_total_mean_biomass
    procedure :: get_splas_mean_shannon_p
    procedure :: get_splas_mean_shannon_z

    procedure :: get_splas_phyto_xp
    procedure :: get_splas_zoo_xz
    procedure :: get_splas_zp_xoptz
    procedure :: get_splas_phyto_mu_0
    procedure :: get_splas_zoo_irz
    procedure :: get_splas_phyto_ks
    procedure :: get_splas_zp_preference_phi

    procedure :: get_splas_nutrient_uptake
    procedure :: get_splas_p_growth_rate
    procedure :: get_splas_grazing_z_p
    procedure :: get_splas_total_p_predation_z
    procedure :: get_splas_p_cell_volume
    procedure :: get_splas_p_v_uptake
    procedure :: get_splas_p_min_cell_quota

    procedure :: get_splas_phyto_carbon
    procedure :: get_splas_phyto_c_export
    procedure :: get_splas_phyto_c_total_export
    procedure :: get_splas_zoo_c_export

    procedure :: get_splas_delta_z
    procedure :: get_splas_number_layers
    procedure :: get_splas_diffusion_coefficient
   
    procedure :: set_splas_irr_profile 
    procedure :: get_splas_irr_profile
    procedure :: set_ivalue    
  end type splas

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  CONSTRUCTOR
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------------!

subroutine init(this)
  CLASS(splas), intent(out) :: this
  print *,"Create SPLAS SIMULATION"

end subroutine init

!------------------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  SIMULATION SUBROUTINE
!
!  subroutine that integrates the differential equations
!  calling procedures from the module solver.f90 that in turn
!  calls each integration procedure from heuns_method.f90 module
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------------!

subroutine integration(this)

 CLASS(splas), intent(inout) :: this

 call eq_solver(this%splas_time_variation,this%splas_1_nutrients,this%splas_1_phyto_class,this%splas_1_zoo_class,this%nutrient_uptake,this%total_p_predation_z, &
                this%z_grazing_p,this%phizp,this%ks,this%irz,this%zoo_half_sat,this%input_supply_nut,this%mu_0,this%zoo_growth_eff, &
                this%zoo_eg_frac,this%phyto_mort,this%zoo_mort,this%delta_t,this%nt,this%print_step, &
                this%splas_nut_limit,this%splas_graz_limit,this%splas_phyto_mean_biomass,this%splas_zoo_mean_biomass, &
                this%splas_phyto_diversity,this%splas_zoo_diversity,this%splas_phyto_total_biomass,this%splas_zoo_total_biomass, &
                this%splas_phyto_small_biomass,this%xp,this%splas_phyto_total_mean_biomass,this%splas_zoo_total_mean_biomass,&
                this%splas_mean_shannon_p,this%splas_mean_shannon_z,&
                this%p_v_uptake_rate,this%p_min_quota,this%p_size_sink,this%splas_1_p_quota,this%p_growth_rate,&
                !CARBON
                this%splas_p_carbon,this%splas_p_c_export,this%splas_z_c_export,this%splas_p_c_tot_export,&
                !SPLAS 1-D
                this%diffusion_coeff,this%delta_z,this%number_layers,&
                !DETRITUS
                this%splas_detritus,&
                !IRRADIANCE 
                this%splas_irradiance,this%ivalue)

end subroutine integration

!------------------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
! CALCULATION OF OTHER QUANTITIES NEEDED LATER IN THE EQUATIONS
!
! calculation is performed using procedures taken by functions.f90
!
! EX: subdivision of size classes for P,Z,MESOZ (log-spaced classes),
!     calculation of optimum prey size for all zooplankton, p uptake rate,
!     p growth rate, z and mesoz ingestion rate, nutrient half saturation for p,
!     p internal minimum quota, p size dependent sinking, p cell volume,
!     preference matrix phi(i,j) for all interactions (z-p, mesoz-p, mesoz-z),
!     initial distribution of p,z,n,mesoz,pquota.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------------!

subroutine calculation(this)
  implicit none
  CLASS(splas), intent(inout) :: this
  integer :: error,i

  allocate(this%xp(1:this%phyto_class_size),this%xz(1:this%zoo_class_size),this%xoptz(1:this%zoo_class_size),&
  this%mu_0(1:this%phyto_class_size),this%irz(1:this%zoo_class_size),&
  this%ks(1:this%phyto_class_size),this%phizp(1:this%phyto_class_size,1:this%zoo_class_size),&
  this%p_cell_volume(1:this%phyto_class_size),this%p_v_uptake_rate(1:this%phyto_class_size),&
  this%p_min_quota(1:this%phyto_class_size),this%p_size_sink(1:this%phyto_class_size),STAT=error)

  if(error/=0) then
    print*, "no memory"
  else
    print*, "correct allocation for constant parameters"
  end if

  allocate(this%nutrient_uptake(1:this%phyto_class_size,1:this%number_layers),this%p_growth_rate(1:this%phyto_class_size,1:this%number_layers),&
           this%total_p_predation_z(1:this%zoo_class_size,1:this%number_layers),&
           this%z_grazing_p(1:this%phyto_class_size,1:this%zoo_class_size,1:this%number_layers),STAT=error)

  if(error/=0) then
    print*, "no memory"
  else
    print*, "correct allocation for time dependent quantities"
  end if

  allocate(this%splas_p_carbon(1:this%phyto_class_size),STAT=error)

  if(error/=0) then
    print*, "no memory"
  else
    print*, "correct allocation for carbon quantities"
  end if

  allocate(this%splas_phyto_small_biomass(1:this%number_layers),this%splas_phyto_total_mean_biomass(1:this%number_layers),&
           this%splas_zoo_total_mean_biomass(1:this%number_layers),this%splas_mean_shannon_p(1:this%number_layers),&
           this%splas_mean_shannon_z(1:this%number_layers),STAT=error)

  if(error/=0) then
    print*, "no memory"
  else
    print*, "correct allocation for 1-D quantities"
  end if

  call logspace_dp(this%phyto_x_min,this%phyto_x_max,this%xp)
  call logspace_dp(this%zoo_x_min,this%zoo_x_max,this%xz)

  call allometric(this%xp,this%xz,this%xoptz,this%p_cell_volume,this%p_v_uptake_rate,this%p_min_quota,&
                  this%mu_0,this%ks,this%irz,this%phyto_sinking,this%p_size_sink)

  call phyto_carbon(this%p_cell_volume,this%splas_p_carbon)

  call calcul_phi(this%xp,this%xoptz,this%zoo_prey_size_tolerance,this%phizp)

  !initialisation of variables
  do i = 1,this%number_layers
    call initialise(this%splas_1_phyto_class(:,1,i),this%splas_1_zoo_class(:,1,i),this%splas_1_nutrients(:,1,i),&
                    this%splas_1_p_quota(:,1,i),i*this%nut_value,this%initial_p,this%initial_z,this%initial_pquota)
    this%diffusion_coeff(i) = this%initial_diff
    this%splas_detritus(1,i) = 0.0E0
    this%splas_irradiance(1,i) = 0.0E0
  end do

end subroutine calculation

!------------------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   SET SUBROUTINES
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------------!

subroutine set_splas_1_nutrients(this,input_array,dim1,dim2,dim3)
! Set nutrients array subroutine
!
! @param input_array : array of nutrients
! @param dim1,dim2,dim3 : size of nutrients array
!
  CLASS(splas), intent(inout) :: this
  integer,intent(in) :: dim1,dim2,dim3
  double precision,intent(in),dimension(dim1,dim2,dim3) :: input_array

  allocate(this%splas_1_nutrients(dim1,dim2,dim3))
  this%splas_1_nutrients = input_array

end subroutine set_splas_1_nutrients

!------------------------------------------------------------------------------!

subroutine set_splas_1_phyto_class(this,input_array,dim1,dim2,dim3)
  ! Set phytoplankton array subroutine
  !
  ! @param input_array : array of phytoplankton (index = size class)
  ! @param dim1,dim2,dim3 : size of phytoplankton array
  !
  CLASS(splas), intent(inout) :: this
  integer,intent(in) :: dim1,dim2,dim3
  double precision,intent(in),dimension(dim1,dim2,dim3) :: input_array

  allocate(this%splas_1_phyto_class(dim1,dim2,dim3))
  this%splas_1_phyto_class = input_array

end subroutine set_splas_1_phyto_class

!------------------------------------------------------------------------------!

subroutine set_splas_1_zoo_class(this,input_array,dim1,dim2,dim3)
  ! Set nano and micro-zooplankton array subroutine
  !
  ! @param input_array : array of zooplankton (index = size class)
  ! @param dim1,dim2,dim3 : size of zooplankton array
  !
  CLASS(splas), intent(inout) :: this
  integer,intent(in) :: dim1,dim2,dim3
  double precision,intent(in),dimension(dim1,dim2,dim3) :: input_array

  allocate(this%splas_1_zoo_class(dim1,dim2,dim3))
  this%splas_1_zoo_class = input_array

end subroutine set_splas_1_zoo_class

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

subroutine set_splas_1_p_quota(this,input_array,dim1,dim2,dim3)
  ! Set internal phytoplankton quota array subroutine
  !
  ! @param input_array : array of phytoplankton cell quota (index = size class)
  ! @param dim1,dim2,dim3 : size of phytoplankton quota array
  !
  CLASS(splas), intent(inout) :: this
  integer,intent(in) :: dim1,dim2,dim3
  double precision,intent(in),dimension(dim1,dim2,dim3) :: input_array

  allocate(this%splas_1_p_quota(dim1,dim2,dim3))
  this%splas_1_p_quota = input_array

end subroutine set_splas_1_p_quota

!------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DETRITUS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------!

subroutine set_splas_1_detritus(this,input_array,dim2,dim3) !for now only time and depth dimensions
  ! Set detritus
  !
  ! @param input_array : array of detritus (time x depth) maybe size will be included in the future (size x time x depth)
  ! @param dim2,dim3 : dimensions of detritus (in the future to include size include it as THE FIRST DIMENSION)
  !
  CLASS(splas), intent(inout) :: this
  integer,intent(in) :: dim2,dim3
  double precision,intent(in),dimension(dim2,dim3) :: input_array

  allocate(this%splas_detritus(dim2,dim3))
  this%splas_detritus = input_array

end subroutine set_splas_1_detritus

!------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_splas_irr_profile(this,input_array,dim2,dim3)


  CLASS(splas), intent(inout) :: this
  integer,intent(in) :: dim2,dim3
  double precision,intent(in),dimension(dim2,dim3) :: input_array
  
  allocate(this%splas_irradiance(dim2,dim3))
  this%splas_irradiance = input_array 

end subroutine set_splas_irr_profile

subroutine set_ivalue(this,input,dim)

  CLASS(splas),intent(inout) :: this
  integer,intent(in) :: dim
  double precision,intent(in),dimension(dim) :: input
  
  allocate(this%ivalue(dim))
  this%ivalue = input

end subroutine set_ivalue

!------------------------------------------------------------------------------!

subroutine set_splas_nut_supply(this,input)
  ! Set nutrients supply subroutine
  !
  ! @param input : nutrient supply rate (uMN/day)
  !
  CLASS(splas), intent(inout) :: this
  double precision,intent(in) :: input

  this%input_supply_nut = input

end subroutine set_splas_nut_supply

!------------------------------------------------------------------------------!

subroutine set_splas_time_variation(this,input)
  ! Set flag of time varying nutrient supply
  !
  ! @param input : flag
  !
  CLASS(splas), intent(inout) :: this
  logical,intent(in) :: input

  this%splas_time_variation = input

end subroutine set_splas_time_variation

!------------------------------------------------------------------------------!

subroutine set_splas_number_nutrients(this,input)
  ! Set nutrients number (dimension of array of nutrients) subroutine
  !
  ! @param input : nutrients number
  !
  CLASS(splas), intent(inout) :: this
  integer,intent(in) :: input

  this%number_nutrients = input

end subroutine set_splas_number_nutrients

!------------------------------------------------------------------------------!
subroutine set_splas_phyto_dim(this,input)
  ! Set phytoplankton array dimension subroutine
  !
  ! @param input : dimension of the array of phytoplankton
  !
  CLASS(splas), intent(inout) :: this
  integer,intent(in) :: input

  this%phyto_class_size = input

end subroutine set_splas_phyto_dim

!------------------------------------------------------------------------------!

subroutine set_splas_phyto_xmin(this,input)
  ! Set minimum size of phytoplankton subroutine
  !
  ! @param input : minimum size of phytoplankton
  !
  CLASS(splas), intent(inout) :: this
  double precision,intent(in) :: input

  this%phyto_x_min = input

end subroutine set_splas_phyto_xmin

!------------------------------------------------------------------------------!

subroutine set_splas_phyto_xmax(this,input)
  ! Set maximum size of phytoplankton subroutine
  !
  ! @param input : maximum size of phytoplankton
  !
  CLASS(splas), intent(inout) :: this
  double precision,intent(in) :: input

  this%phyto_x_max = input

end subroutine set_splas_phyto_xmax

!------------------------------------------------------------------------------!

subroutine set_splas_phyto_mort(this,input)
  ! Set mortality fraction of phytoplankton subroutine (no units)
  !
  ! @param input : mortality fraction of phytoplankton
  !
  CLASS(splas), intent(inout) :: this
  double precision,intent(in) :: input

  this%phyto_mort = input

end subroutine set_splas_phyto_mort

!------------------------------------------------------------------------------!

subroutine set_splas_phyto_sink(this,input)
  ! Set mortality fraction of phytoplankton subroutine due to sinking
  !
  ! @param input : mortality fraction of phytoplankton by sinking
  !
  CLASS(splas), intent(inout) :: this
  double precision,intent(in) :: input

  this%phyto_sinking = input

end subroutine set_splas_phyto_sink

!------------------------------------------------------------------------------!

subroutine set_splas_zoo_dim(this,input)
  ! Set nano/micro-zooplankton array dimension subroutine
  !
  ! @param input : dimension of the array of zooplankton
  !
  CLASS(splas), intent(inout) :: this
  integer,intent(in) :: input

  this%zoo_class_size = input

end subroutine set_splas_zoo_dim

!------------------------------------------------------------------------------!

subroutine set_splas_zoo_xmin(this,input)
  ! Set minimum size of nano/micro-zooplankton subroutine
  !
  ! @param input : minimum size of zooplankton
  !
  CLASS(splas), intent(inout) :: this
  double precision,intent(in) :: input

  this%zoo_x_min = input

end subroutine set_splas_zoo_xmin

!------------------------------------------------------------------------------!

subroutine set_splas_zoo_xmax(this,input)
  ! Set maximum size of nano/micro-zooplankton subroutine
  !
  ! @param input : maximum size of zooplankton
  !
  CLASS(splas), intent(inout) :: this
  double precision,intent(in) :: input

  this%zoo_x_max = input

end subroutine set_splas_zoo_xmax

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

subroutine set_splas_zoo_growth(this,input)
  ! Set zooplankton growth efficiency subroutine (no units)
  !
  ! @param input : growth efficiency of zooplankton
  !
  CLASS(splas), intent(inout) :: this
  double precision,intent(in) :: input

  this%zoo_growth_eff = input

end subroutine set_splas_zoo_growth

!------------------------------------------------------------------------------!

subroutine set_splas_zoo_eg(this,input)
  ! Set zooplankton egestion fraction subroutine (no units)
  !
  ! @param input : egestion fraction of zooplankton
  !
  CLASS(splas), intent(inout) :: this
  double precision,intent(in) :: input

  this%zoo_eg_frac = input

end subroutine set_splas_zoo_eg

!------------------------------------------------------------------------------!

subroutine set_splas_zoo_mort(this,input)
  ! Set nano/micro-zooplankton mortality subroutine (quadratic: 1/(mass x day))
  !
  ! @param input : mortality of zooplankton
  !
  CLASS(splas), intent(inout) :: this
  double precision,intent(in) :: input

  this%zoo_mort = input

end subroutine set_splas_zoo_mort

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

subroutine set_splas_zoo_hs(this,input)
  ! Set nano/micro-zooplankton prey half saturation level subroutine (uMN)
  !
  ! @param input : prey half saturation of zooplankton
  !
  CLASS(splas), intent(inout) :: this
  double precision,intent(in) :: input

  this%zoo_half_sat = input

end subroutine set_splas_zoo_hs

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

subroutine set_splas_z_preysize_tolerance(this,input)
  ! Set nano/micro-zooplankton prey size tolerance for zooplankton subroutine (no units)
  !
  ! @param input : prey size tolerance of zooplankton
  !
  CLASS(splas), intent(inout) :: this
  double precision,intent(in) :: input

  this%zoo_prey_size_tolerance = input

end subroutine set_splas_z_preysize_tolerance

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

subroutine set_splas_nut_value(this,input)
  ! Set initial concentration of N subroutine
  !
  ! @param input : initial N concentration
  !
  CLASS(splas), intent(inout) :: this
  double precision,intent(in) :: input

  this%nut_value = input

end subroutine set_splas_nut_value

!------------------------------------------------------------------------------!

subroutine set_splas_initial_p(this,input)
  ! Set phytoplankton initial uniform abundance subroutine
  !
  ! @param input : initial phytoplankton concentration (if uniform)
  !
  CLASS(splas), intent(inout) :: this
  double precision,intent(in) :: input

  this%initial_p = input

end subroutine set_splas_initial_p

!------------------------------------------------------------------------------!

subroutine set_splas_initial_z(this,input)
  ! Set nano/micro-zooplankton initial uniform abundance subroutine
  !
  ! @param input : initial zooplankton concentration (if uniform)
  !
  CLASS(splas), intent(inout) :: this
  double precision,intent(in) :: input

  this%initial_z = input

end subroutine set_splas_initial_z

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

subroutine set_splas_initial_pquota(this,input)
  ! Set phyto quota initial uniform abundance subroutine
  !
  ! @param input : initial phytoplankton cellular quota (if uniform)
  !
  CLASS(splas), intent(inout) :: this
  double precision,intent(in) :: input

  this%initial_pquota = input

end subroutine set_splas_initial_pquota

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SIMULATION PARAMETERS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------------!

subroutine set_splas_sim_length(this,input)
  ! Set simulation length subroutine (usually the unit is days)
  !
  ! @param input : simulation length
  !
  CLASS(splas), intent(inout) :: this
  double precision,intent(in) :: input

  this%sim_length = input

end subroutine set_splas_sim_length

!------------------------------------------------------------------------------!

subroutine set_splas_delta_t(this,input)
  ! Set time step subroutine (usually the unit is days)
  !
  ! @param input : time step
  !
  CLASS(splas), intent(inout) :: this
  double precision,intent(in) :: input

  this%delta_t = input

end subroutine set_splas_delta_t

!------------------------------------------------------------------------------!

subroutine set_splas_nt(this)
  ! Set the number of time steps in the simulation
  !
  !
  !
  CLASS(splas), intent(inout) :: this

  this%nt = int(this%sim_length/this%delta_t)

end subroutine set_splas_nt

!------------------------------------------------------------------------------!

subroutine set_splas_print_step(this,input)
  ! Set the printing step (each print_step steps the quantity is printed to file/plot/screen)
  !
  ! @param input : printing step
  !
  CLASS(splas), intent(inout) :: this
  integer,intent(in) :: input

  this%print_step = input

end subroutine set_splas_print_step

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!! SPLAS 1-D !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

subroutine set_splas_delta_z(this,input)
  ! Set the step of the vertical space-grid
  !
  ! @param input : delta z
  !
  CLASS(splas), intent(inout) :: this
  double precision,intent(in) :: input

  this%delta_z = input

end subroutine set_splas_delta_z

!------------------------------------------------------------------------------!

subroutine set_splas_number_layers(this,input)
  ! Set the number of vertical layers
  !
  ! @param input : number of layers
  !
  CLASS(splas), intent(inout) :: this
  integer,intent(in) :: input

  this%number_layers = input

end subroutine set_splas_number_layers

!------------------------------------------------------------------------------!

subroutine set_splas_initial_diff(this,input)
  ! Set initial value (uniform) diffusion
  !
  ! @param input : initial diffusion coefficient (if uniform)
  !
  CLASS(splas), intent(inout) :: this
  double precision,intent(in) :: input

  this%initial_diff = input

end subroutine set_splas_initial_diff

!------------------------------------------------------------------------------!

subroutine set_splas_diffusion_coefficient(this,input_array,dim1)
  ! Set diffusion coefficient for phytoplankton (allowed to vary with depth)
  !
  ! @param input_array : array of phytoplankton diffusion coefficient
  ! @param dim1 : size of diffusion coefficient array (dimension = number of layers)
  !
  CLASS(splas), intent(inout) :: this
  integer,intent(in) :: dim1
  double precision,intent(in),dimension(dim1) :: input_array

  allocate(this%diffusion_coeff(dim1))
  this%diffusion_coeff = input_array

end subroutine set_splas_diffusion_coefficient

!------------------------------------------------------------------------------!

!-w-----------------------------------------------------------------------------!!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  OTHER VARIABLES THAT WILL NEED TO BE PRINTED OUT...
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------------!

subroutine set_splas_nut_metrics(this,input_array,dim1,dim2,dim3)
! Set nutrients metrics array subroutine
!
! @param input_array : array of nutrients limitation metrics
! @param dim1,dim2,dim3 : size of nutrients metrics array
!
  CLASS(splas), intent(inout) :: this
  integer,intent(in) :: dim1,dim2,dim3
  double precision,intent(in),dimension(dim1,dim2,dim3) :: input_array

  allocate(this%splas_nut_limit(dim1,dim2,dim3))
  this%splas_nut_limit = input_array

end subroutine set_splas_nut_metrics

!------------------------------------------------------------------------------!

subroutine set_splas_graz_metrics(this,input_array,dim1,dim2,dim3)
! Set grazing metrics array subroutine
!
! @param input_array : array of grazing limitation metrics
! @param dim1,dim2,dim3 : size of grazing metrics array
!
  CLASS(splas), intent(inout) :: this
  integer,intent(in) :: dim1,dim2,dim3
  double precision,intent(in),dimension(dim1,dim2,dim3) :: input_array

  allocate(this%splas_graz_limit(dim1,dim2,dim3))
  this%splas_graz_limit = input_array

end subroutine set_splas_graz_metrics

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

subroutine set_splas_time_mean_mass_p(this,input_array,dim1,dim2)
! Set phytoplankton mean biomass array subroutine
!
! @param input_array : array of mean phytoplankton biomass
! @param dim1,dim2 : size of phytoplankton mean biomass array
!
  CLASS(splas), intent(inout) :: this
  integer,intent(in) :: dim1,dim2
  double precision,intent(in),dimension(dim1,dim2) :: input_array

  allocate(this%splas_phyto_mean_biomass(dim1,dim2))
  this%splas_phyto_mean_biomass = input_array

end subroutine set_splas_time_mean_mass_p

!------------------------------------------------------------------------------!

subroutine set_splas_time_mean_mass_z(this,input_array,dim1,dim2)
! Set nano/micro-zooplankton mean biomass array subroutine
!
! @param input_array : array of mean zooplankton biomass
! @param dim1,dim2 : size of zooplankton mean biomass array
!
  CLASS(splas), intent(inout) :: this
  integer,intent(in) :: dim1,dim2
  double precision,intent(in),dimension(dim1,dim2) :: input_array

  allocate(this%splas_zoo_mean_biomass(dim1,dim2))
  this%splas_zoo_mean_biomass = input_array

end subroutine set_splas_time_mean_mass_z

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

subroutine set_splas_shannon_index_p(this,input_array,dim1,dim2)
! Set phytoplankton diversity subroutine
!
! @param input_array : array of phytoplankton diversity
! @param dim1,dim2 : size of phytoplankton diversity array
!
  CLASS(splas), intent(inout) :: this
  integer,intent(in) :: dim1,dim2
  double precision,intent(in),dimension(dim1,dim2) :: input_array

  allocate(this%splas_phyto_diversity(dim1,dim2))
  this%splas_phyto_diversity = input_array

end subroutine set_splas_shannon_index_p

!------------------------------------------------------------------------------!

subroutine set_splas_shannon_index_z(this,input_array,dim1,dim2)
! Set nano/micro-zooplankton diversity subroutine
!
! @param input_array : array of zooplankton diversity
! @param dim1,dim2 : size of zooplankton diversity array
!
  CLASS(splas), intent(inout) :: this
  integer,intent(in) :: dim1,dim2
  double precision,intent(in),dimension(dim1,dim2) :: input_array

  allocate(this%splas_zoo_diversity(dim1,dim2))
  this%splas_zoo_diversity = input_array

end subroutine set_splas_shannon_index_z

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

subroutine set_splas_total_p_biomass(this,input_array,dim1,dim2)
! Set phytoplankton total biomass subroutine
!
! @param input_array : array of phytoplankton total mass in time
! @param dim1,dim2 : size of phytoplankton total biomass array
!
  CLASS(splas), intent(inout) :: this
  integer,intent(in) :: dim1,dim2
  double precision,intent(in),dimension(dim1,dim2) :: input_array

  allocate(this%splas_phyto_total_biomass(dim1,dim2))
  this%splas_phyto_total_biomass = input_array

end subroutine set_splas_total_p_biomass

!------------------------------------------------------------------------------!

subroutine set_splas_total_z_biomass(this,input_array,dim1,dim2)
! Set nano/micro-zooplankton total biomass subroutine
!
! @param input_array : array of zooplankton total mass in time
! @param dim1,dim2 : size of zooplankton total biomass array
!
  CLASS(splas), intent(inout) :: this
  integer,intent(in) :: dim1,dim2
  double precision,intent(in),dimension(dim1,dim2) :: input_array

  allocate(this%splas_zoo_total_biomass(dim1,dim2))
  this%splas_zoo_total_biomass = input_array

end subroutine set_splas_total_z_biomass

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  GET FUNCTIONS
!
!  in each function the result is get = the corresponding variable
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------------!

function get_splas_1_nutrients(this,idx1,idx2,idx3) result(get)

  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx1,idx2,idx3
  double precision :: get
  get = this%splas_1_nutrients(idx1,idx2,idx3)

end function get_splas_1_nutrients

!------------------------------------------------------------------------------!

function get_splas_1_phyto_class(this,idx1,idx2,idx3) result(get)

  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx1,idx2,idx3
  double precision :: get
  get = this%splas_1_phyto_class(idx1,idx2,idx3)

end function get_splas_1_phyto_class

!------------------------------------------------------------------------------!

function get_splas_1_zoo_class(this,idx1,idx2,idx3) result(get)

  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx1,idx2,idx3
  double precision :: get
  get = this%splas_1_zoo_class(idx1,idx2,idx3)

end function get_splas_1_zoo_class

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

function get_splas_1_p_quota(this,idx1,idx2,idx3) result(get)

  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx1,idx2,idx3
  double precision :: get
  get = this%splas_1_p_quota(idx1,idx2,idx3)

end function get_splas_1_p_quota

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DETRITUS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------!

function get_splas_1_detritus(this,idx2,idx3) result(get)

  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx2,idx3
  double precision :: get
  get = this%splas_detritus(idx2,idx3)

end function get_splas_1_detritus

!------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function get_splas_irr_profile(this,idx2,idx3) result(get)

  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx2,idx3
  double precision :: get
  get = this%splas_irradiance(idx2,idx3)

end function get_splas_irr_profile
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

function get_splas_nut_supply(this) result(get)

  CLASS(splas), intent(in) :: this
  double precision :: get

  get = this%input_supply_nut

end function get_splas_nut_supply

!------------------------------------------------------------------------------!

function get_splas_number_nutrients(this) result(get)

  CLASS(splas), intent(in) :: this
  integer :: get

  get = this%number_nutrients

end function get_splas_number_nutrients

!------------------------------------------------------------------------------!

function get_splas_phyto_dim(this) result(get)

  CLASS(splas), intent(in) :: this
  integer :: get

  get = this%phyto_class_size

end function get_splas_phyto_dim

!------------------------------------------------------------------------------!

function get_splas_phyto_xmin(this) result(get)

  CLASS(splas), intent(in) :: this
  double precision :: get

  get = this%phyto_x_min

end function get_splas_phyto_xmin

!------------------------------------------------------------------------------!

function get_splas_phyto_xmax(this) result(get)

  CLASS(splas), intent(in) :: this
  double precision :: get

  get = this%phyto_x_max

end function get_splas_phyto_xmax

!------------------------------------------------------------------------------!

function get_splas_phyto_mort(this) result(get)

  CLASS(splas), intent(in) :: this
  double precision :: get

  get = this%phyto_mort

end function get_splas_phyto_mort

!------------------------------------------------------------------------------!

function get_splas_phyto_sink(this) result(get)

  CLASS(splas), intent(in) :: this
  double precision :: get

  get = this%phyto_sinking

end function get_splas_phyto_sink

!------------------------------------------------------------------------------!

function get_splas_zoo_dim(this) result(get)

  CLASS(splas), intent(in) :: this
  integer :: get

  get = this%zoo_class_size

end function get_splas_zoo_dim

!------------------------------------------------------------------------------!

function get_splas_zoo_xmin(this) result(get)

  CLASS(splas), intent(in) :: this
  double precision :: get

  get = this%zoo_x_min

end function get_splas_zoo_xmin

!------------------------------------------------------------------------------!

function get_splas_zoo_xmax(this) result(get)

  CLASS(splas), intent(in) :: this
  double precision :: get

  get = this%zoo_x_max

end function get_splas_zoo_xmax

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

function get_splas_zoo_growth(this) result(get)

  CLASS(splas), intent(in) :: this
  double precision :: get

  get = this%zoo_growth_eff

end function get_splas_zoo_growth

!------------------------------------------------------------------------------!

function get_splas_zoo_eg(this) result(get)

  CLASS(splas), intent(in) :: this
  double precision :: get

  get = this%zoo_eg_frac

end function get_splas_zoo_eg

!------------------------------------------------------------------------------!

function get_splas_zoo_mort(this) result(get)

  CLASS(splas), intent(in) :: this
  double precision :: get

  get = this%zoo_mort

end function get_splas_zoo_mort

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

function get_splas_zoo_hs(this) result(get)

  CLASS(splas), intent(in) :: this
  double precision :: get

  get = this%zoo_half_sat

end function get_splas_zoo_hs

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

function get_splas_z_preysize_tolerance(this) result(get)

  CLASS(splas), intent(in) :: this
  double precision :: get

  get = this%zoo_prey_size_tolerance

end function get_splas_z_preysize_tolerance

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

function get_splas_sim_length(this) result(get)

  CLASS(splas), intent(in) :: this
  double precision :: get

  get = this%sim_length

end function get_splas_sim_length

!------------------------------------------------------------------------------!

function get_splas_delta_t(this) result(get)

  CLASS(splas), intent(in) :: this
  double precision :: get

  get = this%delta_t

end function get_splas_delta_t

!------------------------------------------------------------------------------!

function get_splas_nt(this) result(get)

  CLASS(splas), intent(in) :: this
  integer :: get

  get = this%nt

end function get_splas_nt

!------------------------------------------------------------------------------!

function get_splas_print_step(this) result(get)

  CLASS(splas), intent(in) :: this
  integer :: get

  get = this%print_step

end function get_splas_print_step

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!! SPLAS 1-D !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

function get_splas_delta_z(this) result(get)

  CLASS(splas), intent(in) :: this
  integer :: get

  get = this%delta_z

end function get_splas_delta_z

!------------------------------------------------------------------------------!

function get_splas_number_layers(this) result(get)

  CLASS(splas), intent(in) :: this
  integer :: get

  get = this%number_layers

end function get_splas_number_layers

!------------------------------------------------------------------------------!

function get_splas_diffusion_coefficient(this,idx) result(get)
!varying with depth
  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx
  double precision :: get

  get = this%diffusion_coeff(idx)

end function get_splas_diffusion_coefficient

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  OTHER VARIABLES THAT WILL NEED TO BE PRINTED OUT...
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------------!

function get_splas_nut_metrics(this,idx1,idx2,idx3) result(get)
!varying with depth
  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx1,idx2,idx3
  double precision :: get
  get = this%splas_nut_limit(idx1,idx2,idx3)

end function get_splas_nut_metrics

!------------------------------------------------------------------------------!

function get_splas_graz_metrics(this,idx1,idx2,idx3) result(get)
!varying with depth
  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx1,idx2,idx3
  double precision :: get
  get = this%splas_graz_limit(idx1,idx2,idx3)

end function get_splas_graz_metrics

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

function get_splas_time_mean_mass_p(this,idx1,idx2) result(get)
!varying with depth
  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx1,idx2
  double precision :: get
  get = this%splas_phyto_mean_biomass(idx1,idx2)

end function get_splas_time_mean_mass_p

!------------------------------------------------------------------------------!

function get_splas_time_mean_mass_z(this,idx1,idx2) result(get)
!varying with depth
  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx1,idx2
  double precision :: get
  get = this%splas_zoo_mean_biomass(idx1,idx2)

end function get_splas_time_mean_mass_z

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

function get_splas_shannon_index_p(this,idx1,idx2) result(get)
!varying with depth
  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx1,idx2
  double precision :: get
  get = this%splas_phyto_diversity(idx1,idx2)

end function get_splas_shannon_index_p

!------------------------------------------------------------------------------!

function get_splas_shannon_index_z(this,idx1,idx2) result(get)
!varying with depth
  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx1,idx2
  double precision :: get
  get = this%splas_zoo_diversity(idx1,idx2)

end function get_splas_shannon_index_z

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
function get_splas_total_p_biomass(this,idx1,idx2) result(get)
!varying with depth
  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx1,idx2
  double precision :: get
  get = this%splas_phyto_total_biomass(idx1,idx2)

end function get_splas_total_p_biomass

!------------------------------------------------------------------------------!

function get_splas_total_z_biomass(this,idx1,idx2) result(get)
!varying with depth
  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx1,idx2
  double precision :: get
  get = this%splas_zoo_total_biomass(idx1,idx2)

end function get_splas_total_z_biomass

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

function get_splas_phyto_small_biomass(this,idx) result(get)
!varying with depth (1-D)
  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx
  double precision :: get

  get = this%splas_phyto_small_biomass(idx)

end function get_splas_phyto_small_biomass

!------------------------------------------------------------------------------!

function get_splas_phyto_total_mean_biomass(this,idx) result(get)
!varying with depth (1-D)
  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx
  double precision :: get

  get = this%splas_phyto_total_mean_biomass(idx)

end function get_splas_phyto_total_mean_biomass

!------------------------------------------------------------------------------!

function get_splas_zoo_total_mean_biomass(this,idx) result(get)
!varying with depth (1-D)
  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx
  double precision :: get

  get = this%splas_zoo_total_mean_biomass(idx)

end function get_splas_zoo_total_mean_biomass

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

function get_splas_mean_shannon_p(this,idx) result(get)
!varying with depth (1-D)
  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx
  double precision :: get

  get = this%splas_mean_shannon_p(idx)

end function get_splas_mean_shannon_p

!------------------------------------------------------------------------------!

function get_splas_mean_shannon_z(this,idx) result(get)
!varying with depth (1-D)
  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx
  double precision :: get

  get = this%splas_mean_shannon_z(idx)

end function get_splas_mean_shannon_z

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!GET FUNCTIONS FOR THE QUANTITIES COMPUTED IN THE ABOVE SUBROUTINE calculation()
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------------!

function get_splas_phyto_xp(this,idx) result(get)

  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx
  double precision :: get

  get = this%xp(idx)

end function get_splas_phyto_xp

!------------------------------------------------------------------------------!

function get_splas_zoo_xz(this,idx) result(get)

  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx
  double precision :: get

  get = this%xz(idx)

end function get_splas_zoo_xz

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

function get_splas_zp_xoptz(this,idx) result(get)

  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx
  double precision :: get

  get = this%xoptz(idx)

end function get_splas_zp_xoptz

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

function get_splas_phyto_mu_0(this,idx) result(get)

  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx
  double precision :: get

  get = this%mu_0(idx)

end function get_splas_phyto_mu_0

!------------------------------------------------------------------------------!

function get_splas_zoo_irz(this,idx) result(get)

  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx
  double precision :: get

  get = this%irz(idx)

end function get_splas_zoo_irz

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

function get_splas_phyto_ks(this,idx) result(get)

  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx
  double precision :: get

  get = this%ks(idx)

end function get_splas_phyto_ks

!------------------------------------------------------------------------------!

function get_splas_zp_preference_phi(this,idx1,idx2) result(get)

  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx1,idx2
  double precision :: get

  get = this%phizp(idx1,idx2)

end function get_splas_zp_preference_phi

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

function get_splas_p_cell_volume(this,idx) result(get)

  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx
  double precision :: get

  get = this%p_cell_volume(idx)

end function get_splas_p_cell_volume

!------------------------------------------------------------------------------!

function get_splas_p_v_uptake(this,idx) result(get)

  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx
  double precision :: get

  get = this%p_v_uptake_rate(idx)

end function get_splas_p_v_uptake

!------------------------------------------------------------------------------!

function get_splas_p_min_cell_quota(this,idx) result(get)

  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx
  double precision :: get

  get = this%p_min_quota(idx)

end function get_splas_p_min_cell_quota

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! GET FUNCTIONS FOR TIME-DEPENDENT QUANTITIES
! THAT VARY STEP BY STEP DURING THE INTEGRATION, NEEDED IN THE EQUATIONS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------------!

function get_splas_nutrient_uptake(this,idx1,idx2) result(get)
!varying with depth
  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx1,idx2
  double precision :: get

  get = this%nutrient_uptake(idx1,idx2)

end function get_splas_nutrient_uptake

!------------------------------------------------------------------------------!

function get_splas_p_growth_rate(this,idx1,idx2) result(get)
!varying with depth
  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx1,idx2
  double precision :: get

  get = this%p_growth_rate(idx1,idx2)

end function get_splas_p_growth_rate

!------------------------------------------------------------------------------!

function get_splas_grazing_z_p(this,idx1,idx2,idx3) result(get)
!varying with depth
  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx1,idx2,idx3
  double precision :: get

  get = this%z_grazing_p(idx1,idx2,idx3)

end function get_splas_grazing_z_p

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

function get_splas_total_p_predation_z(this,idx1,idx2) result(get)
!varying with depth
  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx1,idx2
  double precision :: get

  get = this%total_p_predation_z(idx1,idx2)

end function get_splas_total_p_predation_z

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SET AND GET FUNCTIONS FOR CARBON RELATED QUANTITIES
! NEEDED TO CONNECT EVERYTHING TO THE CARBON PUMP...
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------------!
!SET
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

subroutine set_splas_phyto_carbon(this,input_array,dim)
  ! Set phytoplankton carbon array subroutine
  !
  ! @param input_array : array of phytoplankton carbon (index = size)
  ! @param dim : size of phytoplankton carbon
  !
  CLASS(splas), intent(inout) :: this
  integer,intent(in) :: dim
  double precision,intent(in),dimension(dim) :: input_array

  allocate(this%splas_p_carbon(dim))
  this%splas_p_carbon = input_array

end subroutine set_splas_phyto_carbon

!------------------------------------------------------------------------------!

subroutine set_splas_phyto_c_export(this,input_array,dim1,dim2)
  ! Set phytoplankton carbon export array subroutine
  !
  ! @param input_array : array of phytoplankton carbon export (index1 = size class, index2 = time)
  ! @param dim1,dim2 : size of phytoplankton export array
  !
  CLASS(splas), intent(inout) :: this
  integer,intent(in) :: dim1,dim2
  double precision,intent(in),dimension(dim1,dim2) :: input_array

  allocate(this%splas_p_c_export(dim1,dim2))
  this%splas_p_c_export = input_array

end subroutine set_splas_phyto_c_export

!------------------------------------------------------------------------------!

subroutine set_splas_phyto_c_total_export(this,input_array,dim)
  ! Set phytoplankton total carbon export array subroutine
  !
  ! @param input_array : array of phytoplankton total carbon export (index = time)
  ! @param dim : size of phytoplankton total export array
  !
  CLASS(splas), intent(inout) :: this
  integer,intent(in) :: dim
  double precision,intent(in),dimension(dim) :: input_array

  allocate(this%splas_p_c_tot_export(dim))
  this%splas_p_c_tot_export = input_array

end subroutine set_splas_phyto_c_total_export

!------------------------------------------------------------------------------!

subroutine set_splas_zoo_c_export(this,input_array,dim)
  ! Set zooplankton total carbon export array subroutine
  !
  ! @param input_array : array of zooplankton total carbon export (index = time)
  ! @param dim : size of zooplankton total export array
  !
  CLASS(splas), intent(inout) :: this
  integer,intent(in) :: dim
  double precision,intent(in),dimension(dim) :: input_array

  allocate(this%splas_z_c_export(dim))
  this%splas_z_c_export = input_array

end subroutine set_splas_zoo_c_export

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!GET
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

function get_splas_phyto_carbon(this,idx) result(get)

  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx
  double precision :: get

  get = this%splas_p_carbon(idx)

end function get_splas_phyto_carbon

!------------------------------------------------------------------------------!

function get_splas_phyto_c_export(this,idx1,idx2) result(get)

  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx1,idx2
  double precision :: get
  get = this%splas_p_c_export(idx1,idx2)

end function get_splas_phyto_c_export

!------------------------------------------------------------------------------!

function get_splas_phyto_c_total_export(this,idx) result(get)

  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx
  double precision :: get

  get = this%splas_p_c_tot_export(idx)

end function get_splas_phyto_c_total_export

!------------------------------------------------------------------------------!

function get_splas_zoo_c_export(this,idx) result(get)

  CLASS(splas), intent(in) :: this
  integer,intent(in) :: idx
  double precision :: get

  get = this%splas_z_c_export(idx)

end function get_splas_zoo_c_export

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

end module splas_class
