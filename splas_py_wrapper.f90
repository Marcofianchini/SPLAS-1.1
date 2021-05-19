module wrapper_py

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  SPLAS python wrapper
!
!  python wrapper for fortran code of Sea PLAnkton Simulator (SPLAS)
!  includes all calls to fortran subroutines
!
!
!  @author : mdepasquale, epascolo
!  @year : 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

USE splas_class

implicit none

private
save
  type(splas) :: simulation

public :: set_nutrient_number,get_nutrient_number
public :: init_wrapper,set_nutrients,get_nutrients_x
public :: set_phytoplankton,get_phytoplankton_x,set_zooplankton,get_zooplankton_x
public :: set_nutrient_supply,get_nutrient_supply
public :: set_phytoplankton_dim,get_phytoplankton_dim,set_zooplankton_dim,get_zooplankton_dim
public :: set_phytoplankton_xmin,get_phytoplankton_xmin,set_phytoplankton_xmax,get_phytoplankton_xmax
public :: set_zooplankton_xmin,get_zooplankton_xmin,set_zooplankton_xmax,get_zooplankton_xmax
public :: set_phytoplankton_mort,get_phytoplankton_mort,set_zoooplankton_mort,get_zooplankton_mort
public :: set_zooplankton_growth,get_zooplankton_growth,set_zooplankton_egestion,get_zooplankton_egestion
public :: set_zoooplankton_half_sat,get_zooplankton_half_sat,set_zoooplankton_size_tol,get_zooplankton_size_tol
public :: set_simulation_length,get_simulation_length,set_time_step,get_time_step
public :: calculate_q
public :: get_phytoplankton_size_xp,get_zooplankton_size_xz,get_zp_optimum_prey_size,get_phytoplankton_growth_rate
public :: get_zooplankton_ingestion_rate,get_phytoplankton_half_sat,get_zp_pref_matrix_elem
public :: get_phytoplankton_size_xp_array,get_zooplankton_size_xz_array,get_zp_optimum_prey_size_array,get_phytoplankton_growth_rate_array
public :: get_zooplankton_ingestion_rate_array,get_phytoplankton_half_sat_array,get_zp_preference_matrix
public :: integrate_eqs
public :: get_nutrient_uptake_p, get_nutrient_uptake_p_array,get_grazing_matrix_elem,get_grazing_matrix,get_total_p_predation,get_total_p_predation_array
public :: get_whole_nutrients_array,get_whole_phyto_array,get_whole_zoo_array,set_nt,get_nt
public :: set_print_step,get_print_step
public :: get_1_nut_all_times
public :: set_nutrient_metrics,get_nut_metrics_x,get_whole_nut_metrics
public :: set_grazing_metrics,get_graz_metrics_x,get_whole_graz_metrics
public :: set_phyto_time_mean,get_phyto_time_mean,get_phyto_time_mean_array
public :: set_zoo_time_mean,get_zoo_time_mean,get_zoo_time_mean_array
public :: set_phyto_shannon_index,get_phyto_shannon_index,get_phyto_sh_idx_array
public :: set_zoo_shannon_index,get_zoo_shannon_index,get_zoo_sh_idx_array
public :: set_phyto_biomass,get_phyto_biomass,get_phyto_biomass_array
public :: set_zoo_biomass,get_zoo_biomass,get_zoo_biomass_array
public :: set_initial_n,set_initial_p,set_initial_z
public :: get_phyto_small_biomass,get_phyto_total_mean_biomass,get_zoo_total_mean_biomass
public :: get_phyto_mean_shannon,get_zoo_mean_shannon
public :: get_phyto_small_biomass_array,get_phyto_total_mean_biomass_array,get_zoo_total_mean_biomass_array
public :: get_phyto_mean_shannon_array,get_zoo_mean_shannon_array
public :: set_time_variation
public :: set_delta_z,set_number_layers,set_diffusion,set_initial_diff
public :: get_delta_z,get_number_layers,get_diffusion,get_diffusion_array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
public :: set_phyto_sinking,get_phyto_sinking
public :: get_pcell_volume,get_pcell_volume_array
public :: get_pv_uptake,get_pv_uptake_array
public :: get_pmin_cell_quota,get_pmin_cell_quota_array
public :: set_pquota,get_pquota_x,get_whole_pquota
public :: set_initial_pquota
public :: get_p_growth_rate,get_p_growth_rate_array
!carbon functions and subroutines
public :: set_p_cexport,set_p_ctotalexport
public :: set_z_cexport
public :: get_p_carbon,get_p_carbon_array
public :: get_p_cexport,get_p_cexport_array,get_p_cexport_time,get_whole_p_cexport
public :: get_p_ctotalexport,get_p_ctotalexport_array
public :: get_z_cexport,get_z_cexport_array
!detritus functions and subroutines
public :: set_detritus,get_detritus_x,get_whole_detritus
!IRRADIANCE fil
public :: set_irradiance_profile,get_irradiance_x,get_irradiance_profile,set_i_val

contains


!------------------------------------------------------------------------------!

subroutine init_wrapper()
! Call SPLAS constructor

  implicit none

  call simulation%init()

end subroutine init_wrapper

!------------------------------------------------------------------------------!

subroutine set_nutrients(input_a,dim1,dim2,dim3)
! Wrapper for set_splas_1_nutrients
!
! @param input_a : array of nutrients
! @param dim1,dim2,dim3 : size of nutrients array
!
  implicit none

  integer,intent(in) :: dim1,dim2,dim3
  double precision,intent(inout),dimension(dim1,dim2,dim3) :: input_a

  call simulation%set_splas_1_nutrients(input_a,dim1,dim2,dim3)

end subroutine set_nutrients

!------------------------------------------------------------------------------!

function get_nutrients_x(idx1,idx2,idx3) result(get)
! Wrapper for get_splas_1_nutrients
!
! @param idx1,idx2,idx3  : indeces of array
! @return get : value of nutrients at idx1 at time idx2 in layer idx3
!
  integer,intent(in) :: idx1,idx2,idx3
  double precision :: get
  get = simulation%get_splas_1_nutrients(idx1,idx2,idx3)

end function get_nutrients_x

!------------------------------------------------------------------------------!

function get_1_nut_all_times(idx1,dm2,dm3) result(get)
! Wrapper for get_splas_1_nutrients that returns one nutrient at all times at all layers
!
! @param dm2,dm3 : dimension of array
! @param idx1 : index of nutrient----> tipically NITROGEN (1)
! @return get : array with the values of nutrient idx1 at all times and in all layers
!
  integer :: idx,idxx
  integer, intent(in) :: dm3,dm2,idx1
  double precision, dimension(dm2,dm3) :: get

  do idx = 1,dm2
    do idxx = 1,dm3
      get(idx,idxx) = simulation%get_splas_1_nutrients(idx1,idx,idxx)
    end do
  end do

end function get_1_nut_all_times

!------------------------------------------------------------------------------!

function get_whole_nutrients_array(dm1,dm2,dm3) result(get)
! Wrapper for get_splas_1_nutrients that returns the whole array
!
! @param dm1,dm2,dm3  : dimension of array
! @return get : array with the values of nutrients at all times and in all layers
!
  integer :: idx1,idx2,idx3
  integer, intent(in) :: dm1,dm2,dm3
  double precision, dimension(dm1,dm2,dm3) :: get

  do idx3 = 1,dm3
    do idx2 = 1,dm2
      do idx1 = 1,dm1
        get(idx1,idx2,idx3) = simulation%get_splas_1_nutrients(idx1,idx2,idx3)
      end do
    end do
  end do

end function get_whole_nutrients_array

!------------------------------------------------------------------------------!

subroutine set_phytoplankton(input_a,dim1,dim2,dim3)
! Wrapper for set_splas_1_phyto_class
!
! @param input_a : array of phytoplankton
! @param dim1,dim2,dim3 : size of phytoplankton array
!
  implicit none

  integer,intent(in) :: dim1,dim2,dim3
  double precision,intent(inout),dimension(dim1,dim2,dim3) :: input_a

  call simulation%set_splas_1_phyto_class(input_a,dim1,dim2,dim3)

end subroutine set_phytoplankton

!------------------------------------------------------------------------------!

function get_phytoplankton_x(idx1,idx2,idx3) result(get)
! Wrapper for get_splas_1_phyto_class
!
! @param idx1,idx2,idx3  : index of array
! @return get : value of phytoplankton at idx1 at time idx2 in layer idx3
!
  integer,intent(in) :: idx1,idx2,idx3
  double precision :: get
  get = simulation%get_splas_1_phyto_class(idx1,idx2,idx3)

end function get_phytoplankton_x

!------------------------------------------------------------------------------!

function get_whole_phyto_array(dm1,dm2,dm3) result(get)
! Wrapper for get_splas_1_phyto_class that returns the whole array
!
! @param dm1,dm2,dm3  : dimension of array
! @return get : array with the values of phytoplankton at all times and in all layers
!
  integer :: idx1,idx2,idx3
  integer, intent(in) :: dm1,dm2,dm3
  double precision, dimension(dm1,dm2,dm3) :: get

  do idx3 = 1,dm3
    do idx2 = 1,dm2
      do idx1 = 1,dm1
        get(idx1,idx2,idx3) = simulation%get_splas_1_phyto_class(idx1,idx2,idx3)
      end do
    end do
  end do

end function get_whole_phyto_array

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

subroutine set_zooplankton(input_a,dim1,dim2,dim3)
! Wrapper for set_splas_1_zoo_class
!
! @param input_a : array of zooplankton
! @param dim1,dim2,dim3 : size of zooplankton array
!
  implicit none

  integer,intent(in) :: dim1,dim2,dim3
  double precision,intent(inout),dimension(dim1,dim2,dim3) :: input_a

  call simulation%set_splas_1_zoo_class(input_a,dim1,dim2,dim3)

end subroutine set_zooplankton

!------------------------------------------------------------------------------!

function get_zooplankton_x(idx1,idx2,idx3) result(get)
! Wrapper for get_splas_1_zoo_class
!
! @param idx1,idx2,idx3  : index of array
! @return get : value of zooplankton at idx1 at time idx2 in layer idx3
!
  integer,intent(in) :: idx1,idx2,idx3
  double precision :: get
  get = simulation%get_splas_1_zoo_class(idx1,idx2,idx3)

end function get_zooplankton_x

!------------------------------------------------------------------------------!

function get_whole_zoo_array(dm1,dm2,dm3) result(get)
! Wrapper for get_splas_1_zoo_class that returns the whole array
!
! @param dm1,dm2,dm3  : dimension of array
! @return get : array with the values of zooplankton at all times and in all layers
!
  integer :: idx1,idx2,idx3
  integer, intent(in) :: dm1,dm2,dm3
  double precision, dimension(dm1,dm2,dm3) :: get

  do idx3 = 1,dm3
    do idx2 = 1,dm2
      do idx1 = 1,dm1
        get(idx1,idx2,idx3) = simulation%get_splas_1_zoo_class(idx1,idx2,idx3)
      end do
    end do
  end do

end function get_whole_zoo_array

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

subroutine set_pquota(input_a,dim1,dim2,dim3)
! Wrapper for set_splas_1_p_quota
!
! @param input_a : array of internal nutrient quota in phytoplankton
! @param dim1,dim2,dim3 : size of cellular quota array
!
  implicit none

  integer,intent(in) :: dim1,dim2,dim3
  double precision,intent(inout),dimension(dim1,dim2,dim3) :: input_a

  call simulation%set_splas_1_p_quota(input_a,dim1,dim2,dim3)

end subroutine set_pquota

!------------------------------------------------------------------------------!

function get_pquota_x(idx1,idx2,idx3) result(get)
! Wrapper for get_splas_1_p_quota
!
! @param idx1,idx2,idx3 : indeces of array
! @return get : value of internal quota of phyto idx1 at time idx2 in layer idx3
!
  integer,intent(in) :: idx1,idx2,idx3
  double precision :: get
  get = simulation%get_splas_1_p_quota(idx1,idx2,idx3)

end function get_pquota_x

!------------------------------------------------------------------------------!

function get_whole_pquota(dm1,dm2,dm3) result(get)
! Wrapper for get_splas_1_p_quota that returns the whole array
!
! @param dm1,dm2,dm3  : dimension of array
! @return get : array with the values of p cell quota at all times and in all layers
!
  integer :: idx1,idx2,idx3
  integer, intent(in) :: dm1,dm2,dm3
  double precision, dimension(dm1,dm2,dm3) :: get

  do idx3 = 1,dm3
    do idx2 = 1,dm2
      do idx1 = 1,dm1
        get(idx1,idx2,idx3) = simulation%get_splas_1_p_quota(idx1,idx2,idx3)
      end do
    end do
  end do

end function get_whole_pquota

!------------------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DETRITUS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------!

subroutine set_detritus(input_a,dim2,dim3)
! Wrapper for set_splas_1_detritus
!
! @param input_a : array of detritus
! @param dim2,dim3 : size of detritus array
!
  implicit none

  integer,intent(in) :: dim2,dim3
  double precision,intent(inout),dimension(dim2,dim3) :: input_a

  call simulation%set_splas_1_detritus(input_a,dim2,dim3)

end subroutine set_detritus

!------------------------------------------------------------------------------!

function get_detritus_x(idx2,idx3) result(get)
! Wrapper for get_splas_1_detritus
!
! @param idx2,idx3  : index of array
! @return get : value of detritus at time idx2 in layer idx3
!
  integer,intent(in) :: idx2,idx3
  double precision :: get
  get = simulation%get_splas_1_detritus(idx2,idx3)

end function get_detritus_x

!------------------------------------------------------------------------------!

function get_whole_detritus(dm2,dm3) result(get)
! Wrapper for get_splas_1_detritus that returns the whole array
!
! @param dm2,dm3  : dimension of array
! @return get : array with the values of detritus at all times and in all layers
!
  integer :: idx2,idx3
  integer, intent(in) :: dm2,dm3
  double precision, dimension(dm2,dm3) :: get

  do idx3 = 1,dm3
    do idx2 = 1,dm2
      !do idx1 = 1,dm1
      get(idx2,idx3) = simulation%get_splas_1_detritus(idx2,idx3)
      !end do
    end do
  end do

end function get_whole_detritus

!------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------------!

subroutine set_nutrient_supply(input)
! Wrapper for set_splas_nut_supply
!
! @param input: nutrient supply rate (micro moles of N per day usually)
!
  implicit none

  double precision,intent(inout) :: input

  call simulation%set_splas_nut_supply(input)

end subroutine set_nutrient_supply

!------------------------------------------------------------------------------!

function get_nutrient_supply() result(get)
! Wrapper for get_splas_nut_supply
!
! @return get : value of nutrient supply rate
!
  double precision :: get
  get = simulation%get_splas_nut_supply()

end function get_nutrient_supply

!------------------------------------------------------------------------------!

subroutine set_time_variation(input)
! Wrapper for set_splas_time_variation
!
! @param input: flag for time variation of nutrient supply
!
  implicit none

  logical,intent(inout) :: input

  call simulation%set_splas_time_variation(input)

end subroutine set_time_variation

!------------------------------------------------------------------------------!

subroutine set_nutrient_number(input)
! Wrapper for set_splas_number_nutrients
!
! @param input: number of nutrient elements
!
  implicit none

  integer,intent(inout) :: input

  call simulation%set_splas_number_nutrients(input)

end subroutine set_nutrient_number

!------------------------------------------------------------------------------!

function get_nutrient_number() result(get)
! Wrapper for get_splas_number_nutrients
!
! @return get : value of nutrient number
!
  integer :: get
  get = simulation%get_splas_number_nutrients()

end function get_nutrient_number

!------------------------------------------------------------------------------!

subroutine set_phytoplankton_dim(input)
! Wrapper for set_splas_phyto_dim
!
! @param input: number of size classes for phytoplankton (dimension of array)
!
  implicit none

  integer,intent(inout) :: input

  call simulation%set_splas_phyto_dim(input)

end subroutine set_phytoplankton_dim

!------------------------------------------------------------------------------!

function get_phytoplankton_dim() result(get)
! Wrapper for get_splas_phyto_dim
!
! @return get : number of size classes of phytoplankton
!
  integer :: get
  get = simulation%get_splas_phyto_dim()

end function get_phytoplankton_dim

!------------------------------------------------------------------------------!

subroutine set_phytoplankton_xmin(input)
! Wrapper for set_splas_phyto_xmin
!
! @param input: minimum size for phytoplankton
!
  implicit none

  double precision,intent(inout) :: input

  call simulation%set_splas_phyto_xmin(input)

end subroutine set_phytoplankton_xmin

!------------------------------------------------------------------------------!

function get_phytoplankton_xmin() result(get)
! Wrapper for get_splas_phyto_xmin
!
! @return get : minimum size for phytoplankton
!
  double precision :: get
  get = simulation%get_splas_phyto_xmin()

end function get_phytoplankton_xmin

!------------------------------------------------------------------------------!

subroutine set_phytoplankton_xmax(input)
! Wrapper for set_splas_phyto_xmax
!
! @param input: maximum size for phytoplankton
!
  implicit none

  double precision,intent(inout) :: input

  call simulation%set_splas_phyto_xmax(input)

end subroutine set_phytoplankton_xmax

!------------------------------------------------------------------------------!

function get_phytoplankton_xmax() result(get)
! Wrapper for get_splas_phyto_xmax
!
! @return get : maximum size for phytoplankton
!
  double precision :: get
  get = simulation%get_splas_phyto_xmax()

end function get_phytoplankton_xmax

!------------------------------------------------------------------------------!

subroutine set_phytoplankton_mort(input)
! Wrapper for set_splas_phyto_mort
!
! @param input: mortality fraction for phytoplankton (dimensionless)
!
  implicit none

  double precision,intent(inout) :: input

  call simulation%set_splas_phyto_mort(input)

end subroutine set_phytoplankton_mort

!------------------------------------------------------------------------------!

function get_phytoplankton_mort() result(get)
! Wrapper for get_splas_phyto_mort
!
! @return get : mortality fraction for phytoplankton
!
  double precision :: get
  get = simulation%get_splas_phyto_mort()

end function get_phytoplankton_mort

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

subroutine set_phyto_sinking(input)
! Wrapper for set_splas_phyto_sink
!
! @param input: sinking mortality for phytoplankton
!
  implicit none

  double precision,intent(inout) :: input

  call simulation%set_splas_phyto_sink(input)

end subroutine set_phyto_sinking

!------------------------------------------------------------------------------!

function get_phyto_sinking() result(get)
! Wrapper for get_splas_phyto_sink
!
! @return get : sinking mortality for phytoplankton
!
  double precision :: get
  get = simulation%get_splas_phyto_sink()

end function get_phyto_sinking

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

subroutine set_zooplankton_dim(input)
! Wrapper for set_splas_zoo_dim
!
! @param input: number of size classes for zooplankton (dimension of array)
!
  implicit none

  integer,intent(inout) :: input

  call simulation%set_splas_zoo_dim(input)

end subroutine set_zooplankton_dim

!------------------------------------------------------------------------------!

function get_zooplankton_dim() result(get)
! Wrapper for get_splas_zoo_dim
!
! @return get : number of size classes of zooplankton
!
  integer :: get
  get = simulation%get_splas_zoo_dim()

end function get_zooplankton_dim

!------------------------------------------------------------------------------!

subroutine set_zooplankton_xmin(input)
! Wrapper for set_splas_zoo_xmin
!
! @param input: minimum size for zooplankton
!
  implicit none

  double precision,intent(inout) :: input

  call simulation%set_splas_zoo_xmin(input)

end subroutine set_zooplankton_xmin

!------------------------------------------------------------------------------!

function get_zooplankton_xmin() result(get)
! Wrapper for get_splas_zoo_xmin
!
! @return get : minimum size for zooplankton
!
  double precision :: get
  get = simulation%get_splas_zoo_xmin()

end function get_zooplankton_xmin

!------------------------------------------------------------------------------!

subroutine set_zooplankton_xmax(input)
! Wrapper for set_splas_zoo_xmax
!
! @param input: maximum size for zooplankton
!
  implicit none

  double precision,intent(inout) :: input

  call simulation%set_splas_zoo_xmax(input)

end subroutine set_zooplankton_xmax

!------------------------------------------------------------------------------!

function get_zooplankton_xmax() result(get)
! Wrapper for get_splas_zoo_xmax
!
! @return get : maximum size for zooplankton
!
  double precision :: get
  get = simulation%get_splas_zoo_xmax()

end function get_zooplankton_xmax

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

subroutine set_zooplankton_growth(input)
! Wrapper for set_splas_zoo_growth
!
! @param input: growth efficiency for zooplankton (dimensionless)
!
  implicit none

  double precision,intent(inout) :: input

  call simulation%set_splas_zoo_growth(input)

end subroutine set_zooplankton_growth

!------------------------------------------------------------------------------!

function get_zooplankton_growth() result(get)
! Wrapper for get_splas_zoo_growth
!
! @return get : growth efficiency for zooplankton (dimensionless)
!
  double precision :: get
  get = simulation%get_splas_zoo_growth()

end function get_zooplankton_growth

!------------------------------------------------------------------------------!

subroutine set_zooplankton_egestion(input)
! Wrapper for set_splas_zoo_eg
!
! @param input: egestion fraction for zooplankton (dimensionless)
!
  implicit none

  double precision,intent(inout) :: input

  call simulation%set_splas_zoo_eg(input)

end subroutine set_zooplankton_egestion

!------------------------------------------------------------------------------!

function get_zooplankton_egestion() result(get)
! Wrapper for get_splas_zoo_eg
!
! @return get : egestion fraction for zooplankton (dimensionless)
!
  double precision :: get
  get = simulation%get_splas_zoo_eg()

end function get_zooplankton_egestion

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

subroutine set_zoooplankton_mort(input)
! Wrapper for set_splas_zoo_mort
!
! @param input: mortality for zooplankton ( NOT dimensionless, not a fraction)
!
  implicit none

  double precision,intent(inout) :: input

  call simulation%set_splas_zoo_mort(input)

end subroutine set_zoooplankton_mort

!------------------------------------------------------------------------------!

function get_zooplankton_mort() result(get)
! Wrapper for get_splas_zoo_mort
!
! @return get : mortality for zooplankton ( NOT dimensionless, not a fraction)
!
  double precision :: get
  get = simulation%get_splas_zoo_mort()

end function get_zooplankton_mort

!------------------------------------------------------------------------------!

subroutine set_zoooplankton_half_sat(input)
! Wrapper for set_splas_zoo_hs
!
! @param input: prey half saturation level for zooplankton
!
  implicit none

  double precision,intent(inout) :: input

  call simulation%set_splas_zoo_hs(input)

end subroutine set_zoooplankton_half_sat

!------------------------------------------------------------------------------!

function get_zooplankton_half_sat() result(get)
! Wrapper for get_splas_zoo_hs
!
! @return get : prey half saturation level for zooplankton
!
  double precision :: get
  get = simulation%get_splas_zoo_hs()

end function get_zooplankton_half_sat

!------------------------------------------------------------------------------!

subroutine set_zoooplankton_size_tol(input)
! Wrapper for set_splas_z_preysize_tolerance
!
! @param input: prey size tolerance (choosyness) for zooplankton
!
  implicit none

  double precision,intent(inout) :: input

  call simulation%set_splas_z_preysize_tolerance(input)

end subroutine set_zoooplankton_size_tol

!------------------------------------------------------------------------------!

function get_zooplankton_size_tol() result(get)
! Wrapper for get_splas_z_preysize_tolerance
!
! @return get : prey size tolerance (choosyness) for zooplankton
!
  double precision :: get
  get = simulation%get_splas_z_preysize_tolerance()

end function get_zooplankton_size_tol

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
subroutine set_irradiance_profile(input_a,dim2,dim3)

  implicit none

  integer,intent(in) :: dim2,dim3
  double precision,intent(inout),dimension(dim2,dim3) :: input_a
  
  call simulation%set_splas_irr_profile(input_a,dim2,dim3)

end subroutine set_irradiance_profile

subroutine set_i_val(input,dim)

  implicit none
  integer,intent(in) :: dim
  double precision,intent(inout),dimension(dim) ::input

  call simulation%set_ivalue(input,dim)

end subroutine set_i_val

!------------------------------------------------------------------------------!
function get_irradiance_x(idx2,idx3) result(get)

  integer,intent(in) :: idx2,idx3
  double precision :: get
  get = simulation%get_splas_irr_profile(idx2,idx3)

end function get_irradiance_x

!------------------------------------------------------------------------------!
function get_irradiance_profile(dm2,dm3) result(get)

  integer :: idx2,idx3
  integer, intent(in) :: dm2,dm3
  double precision, dimension(dm2,dm3) :: get

  do idx3 = 1,dm3
    do idx2 = 1,dm2
      get(idx2,idx3) = simulation%get_splas_irr_profile(idx2,idx3)
    end do
  end do

end function get_irradiance_profile
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

subroutine set_initial_n(input)
! Wrapper for set_splas_nut_value
!
! @param input: initial N concentration
!
  implicit none

  double precision,intent(inout) :: input

  call simulation%set_splas_nut_value(input)

end subroutine set_initial_n

!------------------------------------------------------------------------------!

subroutine set_initial_p(input)
! Wrapper for set_splas_initial_p
!
! @param input: initial phytoplankton concentration
!
  implicit none

  double precision,intent(inout) :: input

  call simulation%set_splas_initial_p(input)

end subroutine set_initial_p

!------------------------------------------------------------------------------!

subroutine set_initial_z(input)
! Wrapper for set_splas_initial_z
!
! @param input: initial zooplankton concentration
!
  implicit none

  double precision,intent(inout) :: input

  call simulation%set_splas_initial_z(input)

end subroutine set_initial_z

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

subroutine set_initial_pquota(input)
! Wrapper for set_splas_initial_pquota
!
! @param input: initial phytoplankton internal quota
!
  implicit none

  double precision,intent(inout) :: input

  call simulation%set_splas_initial_pquota(input)

end subroutine set_initial_pquota

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!SIMULATION
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

subroutine set_simulation_length(input)
! Wrapper for set_splas_sim_length
!
! @param input: length of the simulation (usually in days)
!
  implicit none

  double precision,intent(inout) :: input

  call simulation%set_splas_sim_length(input)

end subroutine set_simulation_length

!------------------------------------------------------------------------------!

function get_simulation_length() result(get)
! Wrapper for get_splas_sim_length
!
! @return get : length of the simulation (days)
!
  double precision :: get
  get = simulation%get_splas_sim_length()

end function get_simulation_length

!------------------------------------------------------------------------------!

subroutine set_time_step(input)
! Wrapper for set_splas_delta_t
!
! @param input: time step of the simulation (usually in days)
!
  implicit none

  double precision,intent(inout) :: input

  call simulation%set_splas_delta_t(input)

end subroutine set_time_step

!------------------------------------------------------------------------------!

function get_time_step() result(get)
! Wrapper for get_splas_delta_t
!
! @return get : time step of the simulation (days)
!
  double precision :: get
  get = simulation%get_splas_delta_t()

end function get_time_step

!------------------------------------------------------------------------------!

subroutine set_nt()
! Wrapper for set_splas_nt
!
! @param input: number of time steps of the simulation
!
  implicit none

  call simulation%set_splas_nt()

end subroutine set_nt

!------------------------------------------------------------------------------!

function get_nt() result(get)
! Wrapper for get_splas_nt
!
! @return get : number of time steps of the simulation
!
  integer :: get
  get = simulation%get_splas_nt()

end function get_nt

!------------------------------------------------------------------------------!

subroutine set_print_step(input)
! Wrapper for set_splas_print_step
!
! @param input: writing step of the simulation
!
  implicit none

  integer,intent(inout) :: input

  call simulation%set_splas_print_step(input)

end subroutine set_print_step

!------------------------------------------------------------------------------!

function get_print_step() result(get)
! Wrapper for get_splas_print_step
!
! @return get : writing step of the simulation
!
  integer :: get
  get = simulation%get_splas_print_step()

end function get_print_step

!------------------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VERTICAL SPACE-GRID!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------------!

subroutine set_delta_z(input)
! Wrapper for set_splas_delta_z
!
! @param input: vertical grid step of the simulation (usually in meters)
!
  implicit none

  double precision,intent(inout) :: input

  call simulation%set_splas_delta_z(input)

end subroutine set_delta_z

!------------------------------------------------------------------------------!

function get_delta_z() result(get)
! Wrapper for get_splas_delta_z
!
! @return get : vertical space step of the simulation
!
  double precision :: get
  get = simulation%get_splas_delta_z()

end function get_delta_z

!------------------------------------------------------------------------------!

subroutine set_number_layers(input)
! Wrapper for set_splas_number_layers
!
! @param input: number of vertical layers
!
  implicit none

  integer,intent(inout) :: input

  call simulation%set_splas_number_layers(input)

end subroutine set_number_layers

!------------------------------------------------------------------------------!

function get_number_layers() result(get)
! Wrapper for get_splas_number_layers
!
! @return get : number of vertical layers (that imply overall depth)
!
  double precision :: get
  get = simulation%get_splas_number_layers()

end function get_number_layers

!------------------------------------------------------------------------------!

subroutine set_diffusion(input_a,dim1)
! Wrapper for set_splas_diffusion_coefficient
!
! @param input_a : array of phytoplankton diffusion coefficient
! @param dim1 : number of layers
! in general depending on depth (layers)

  implicit none

  integer,intent(in) :: dim1
  double precision,intent(inout),dimension(dim1) :: input_a

  call simulation%set_splas_diffusion_coefficient(input_a,dim1)

end subroutine set_diffusion

!------------------------------------------------------------------------------!

function get_diffusion(idx) result(get)
! Wrapper for get_splas_diffusion_coefficient
!
! @return get : diffusion coefficient at single layer
! in general dependent on depth (number of layers)

  integer,intent(in) :: idx
  double precision :: get

  get = simulation%get_splas_diffusion_coefficient(idx)

end function get_diffusion

!------------------------------------------------------------------------------!

function get_diffusion_array(dm) result(get)
! Wrapper for get_splas_diffusion_coefficient for the whole array (all layers)
!
! @return get : diffusion coefficient with depth
! in general dependent on depth (number of layers)

  integer :: idx
  integer,intent(in) :: dm
  double precision, dimension(dm) :: get

  do idx = 1,dm
  get(idx) = simulation%get_splas_diffusion_coefficient(idx)
  end do

end function get_diffusion_array

!------------------------------------------------------------------------------!

subroutine set_initial_diff(input)
! Wrapper for set_splas_initial_diff
!
! @param input: initial constant value of diffusion
!
  implicit none

  double precision,intent(inout) :: input

  call simulation%set_splas_initial_diff(input)

end subroutine set_initial_diff

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! WRAPPERS FOR OTHER OUTPUT QUANTITIES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------------!

subroutine set_nutrient_metrics(input_a,dim1,dim2,dim3)
! Wrapper for set_splas_nut_metrics
!
! @param input_a : array of nutrient limitation metrics
! @param dim1,dim2,dim3 : size of nutrient metrics array
! now depending on depth (layers)

  implicit none

  integer,intent(in) :: dim1,dim2,dim3
  double precision,intent(inout),dimension(dim1,dim2,dim3) :: input_a

  call simulation%set_splas_nut_metrics(input_a,dim1,dim2,dim3)

end subroutine set_nutrient_metrics

!------------------------------------------------------------------------------!

function get_nut_metrics_x(idx1,idx2,idx3) result(get)
! Wrapper for get_splas_nut_metrics
!
! @param idx1,idx2,idx3  : indeces of array
! @return get : value of nutrient limitation metrics at idx1 at time idx2 at layer idx3
! now depending on depth (layers)

  integer,intent(in) :: idx1,idx2,idx3
  double precision :: get
  get = simulation%get_splas_nut_metrics(idx1,idx2,idx3)

end function get_nut_metrics_x

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

function get_whole_nut_metrics(dm1,dm2,dm3) result(get)
! Wrapper for get_splas_nut_metrics that returns the whole array
!
! @param dm1,dm2,dm3  : dimension of array
! @return get : array with the values of nutrient limitation at all times at all layers
! now depending on depth (layers)

  integer :: idx1,idx2,idx3
  integer, intent(in) :: dm1,dm2,dm3
  double precision, dimension(dm1,dm2,dm3) :: get

  do idx3 = 1,dm3
    do idx2 = 1,dm2
      do idx1 = 1,dm1
        get(idx1,idx2,idx3) = simulation%get_splas_nut_metrics(idx1,idx2,idx3)
      end do
    end do
  end do

end function get_whole_nut_metrics

!------------------------------------------------------------------------------!

subroutine set_grazing_metrics(input_a,dim1,dim2,dim3)
! Wrapper for set_splas_graz_metrics
!
! @param input_a : array of grazing limitation metrics
! @param dim1,dim2,dim3 : size of grazing metrics array
! now depending on depth (layers)

  implicit none

  integer,intent(in) :: dim1,dim2,dim3
  double precision,intent(inout),dimension(dim1,dim2,dim3) :: input_a

  call simulation%set_splas_graz_metrics(input_a,dim1,dim2,dim3)

end subroutine set_grazing_metrics

!------------------------------------------------------------------------------!

function get_graz_metrics_x(idx1,idx2,idx3) result(get)
! Wrapper for get_splas_graz_metrics
!
! @param idx1,idx2,idx3  : indeces of array
! @return get : value of grazing limitation metrics at idx1 at time idx2 at layer idx3
! now depending on depth (layers)

  integer,intent(in) :: idx1,idx2,idx3
  double precision :: get
  get = simulation%get_splas_graz_metrics(idx1,idx2,idx3)

end function get_graz_metrics_x

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

function get_whole_graz_metrics(dm1,dm2,dm3) result(get)
! Wrapper for get_splas_graz_metrics that returns the whole array
!
! @param dm1,dm2 ,dm3 : dimension of array
! @return get : array with the values of grazing limitation at all times
! now depending on depth (layers)

  integer :: idx1,idx2,idx3
  integer, intent(in) :: dm1,dm2,dm3
  double precision, dimension(dm1,dm2,dm3) :: get

  do idx3 = 1,dm3
    do idx2 = 1,dm2
      do idx1 = 1,dm1
        get(idx1,idx2,idx3) = simulation%get_splas_graz_metrics(idx1,idx2,idx3)
      end do
    end do
  end do

end function get_whole_graz_metrics

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

subroutine set_phyto_time_mean(input_a,dim1,dim2)
! Wrapper for set_splas_time_mean_mass_p
!
! @param input_a : array of phytoplankton time average biomass
! @param dim1,dim2 : size of phytoplankton mean biomass array
! now depending on depth (layers)

  implicit none

  integer,intent(in) :: dim1,dim2
  double precision,intent(inout),dimension(dim1,dim2) :: input_a

  call simulation%set_splas_time_mean_mass_p(input_a,dim1,dim2)

end subroutine set_phyto_time_mean

!------------------------------------------------------------------------------!

function get_phyto_time_mean(idx1,idx2) result(get)
! Wrapper for get_splas_time_mean_mass_p
!
! @param idx1,idx2  : index of array
! @return get : value of phytoplankton mean biomass for size idx1 at layer idx2
! now depending on depth (layers)

  integer,intent(in) :: idx1,idx2
  double precision :: get

  get = simulation%get_splas_time_mean_mass_p(idx1,idx2)

end function get_phyto_time_mean

!------------------------------------------------------------------------------!

function get_phyto_time_mean_array(dm1,dm2) result(get)
! Wrapper for get_splas_time_mean_mass_p that returns the whole array
!
! @param dm1,dm2  : dimension of array
! @return get : array with the values of phytoplankton mean biomass for each size
! now depending on depth (layers)

  integer :: idx1,idx2
  integer, intent(in) :: dm1,dm2
  double precision, dimension(dm1,dm2) :: get

  do idx2 = 1,dm2
    do idx1 = 1,dm1
      get(idx1,idx2) = simulation%get_splas_time_mean_mass_p(idx1,idx2)
    end do
  end do

end function get_phyto_time_mean_array

!------------------------------------------------------------------------------!

subroutine set_zoo_time_mean(input_a,dim1,dim2)
! Wrapper for set_splas_time_mean_mass_z
!
! @param input_a : array of zooplankton time average biomass
! @param dim1,dim2 : size of zooplankton mean biomass array
! now depending on depth (layers)

  implicit none

  integer,intent(in) :: dim1,dim2
  double precision,intent(inout),dimension(dim1,dim2) :: input_a

  call simulation%set_splas_time_mean_mass_z(input_a,dim1,dim2)

end subroutine set_zoo_time_mean

!------------------------------------------------------------------------------!

function get_zoo_time_mean(idx1,idx2) result(get)
! Wrapper for get_splas_time_mean_mass_z
!
! @param idx1,idx2  : index of array
! @return get : value of zooplankton mean biomass for size idx1 at layer idx2
! now depending on depth (layers)

  integer,intent(in) :: idx1,idx2
  double precision :: get

  get = simulation%get_splas_time_mean_mass_z(idx1,idx2)

end function get_zoo_time_mean

!------------------------------------------------------------------------------!

function get_zoo_time_mean_array(dm1,dm2) result(get)
! Wrapper for get_splas_time_mean_mass_z that returns the whole array
!
! @param dm1,dm2  : dimension of array
! @return get : array with the values of zooplankton mean biomass for each size
! now depending on depth (layers)

  integer :: idx1,idx2
  integer, intent(in) :: dm1,dm2
  double precision, dimension(dm1,dm2) :: get

  do idx2 = 1,dm2
    do idx1 = 1,dm1
      get(idx1,idx2) = simulation%get_splas_time_mean_mass_z(idx1,idx2)
    end do
  end do

end function get_zoo_time_mean_array

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

subroutine set_phyto_shannon_index(input_a,dim1,dim2)
! Wrapper for set_splas_shannon_index_p
!
! @param input_a : array of phytoplankton evenness
! @param dim1,dim2 : size of phytoplankton evenness array
! now depending on depth (layers)

  implicit none

  integer,intent(in) :: dim1,dim2
  double precision,intent(inout),dimension(dim1,dim2) :: input_a

  call simulation%set_splas_shannon_index_p(input_a,dim1,dim2)

end subroutine set_phyto_shannon_index

!------------------------------------------------------------------------------!

function get_phyto_shannon_index(idx1,idx2) result(get)
! Wrapper for get_splas_shannon_index_p
!
! @param idx1,idx2  : index of array
! @return get : value of phytoplankton evenness at time idx1 in layer idx2
! now depending on depth (layers)

  integer,intent(in) :: idx1,idx2
  double precision :: get

  get = simulation%get_splas_shannon_index_p(idx1,idx2)

end function get_phyto_shannon_index

!------------------------------------------------------------------------------!

function get_phyto_sh_idx_array(dm1,dm2) result(get)
! Wrapper for get_splas_shannon_index_p that returns the whole array
!
! @param dm1,dm2  : dimension of array
! @return get : array with the values of phytoplankton evenness for each time
! now depending on depth (layers)

  integer :: idx1,idx2
  integer, intent(in) :: dm1,dm2
  double precision, dimension(dm1,dm2) :: get

  do idx2 = 1,dm2
    do idx1 = 1,dm1
      get(idx1,idx2) = simulation%get_splas_shannon_index_p(idx1,idx2)
    end do
  end do

end function get_phyto_sh_idx_array

!------------------------------------------------------------------------------!

subroutine set_zoo_shannon_index(input_a,dim1,dim2)
! Wrapper for set_splas_shannon_index_z
!
! @param input_a : array of zooplankton evenness
! @param dim1,dim2 : size of zooplankton evenness array
! now depending on depth (layers)

  implicit none

  integer,intent(in) :: dim1,dim2
  double precision,intent(inout),dimension(dim1,dim2) :: input_a

  call simulation%set_splas_shannon_index_z(input_a,dim1,dim2)

end subroutine set_zoo_shannon_index

!------------------------------------------------------------------------------!

function get_zoo_shannon_index(idx1,idx2) result(get)
! Wrapper for get_splas_shannon_index_z
!
! @param idx1,idx2  : index of array
! @return get : value of zooplankton evenness at time idx1 at layer idx2
! now depending on depth (layers)

  integer,intent(in) :: idx1,idx2
  double precision :: get

  get = simulation%get_splas_shannon_index_z(idx1,idx2)

end function get_zoo_shannon_index

!------------------------------------------------------------------------------!

function get_zoo_sh_idx_array(dm1,dm2) result(get)
! Wrapper for get_splas_shannon_index_z that returns the whole array
!
! @param dm1,dm2  : dimension of array
! @return get : array with the values of zooplankton evenness for each time
! now depending on depth (layers)

  integer :: idx1,idx2
  integer, intent(in) :: dm1,dm2
  double precision, dimension(dm1,dm2) :: get

  do idx2 = 1,dm2
    do idx1 = 1,dm1
      get(idx1,idx2) = simulation%get_splas_shannon_index_z(idx1,idx2)
    end do
  end do

end function get_zoo_sh_idx_array

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

subroutine set_phyto_biomass(input_a,dim1,dim2)
! Wrapper for set_splas_total_p_biomass
!
! @param input_a : array of phytoplankton total biomass in time
! @param dim1,dim2 : size of phytoplankton total biomass array
! now depending on depth (layers)

  implicit none

  integer,intent(in) :: dim1,dim2
  double precision,intent(inout),dimension(dim1,dim2) :: input_a

  call simulation%set_splas_total_p_biomass(input_a,dim1,dim2)

end subroutine set_phyto_biomass

!------------------------------------------------------------------------------!

function get_phyto_biomass(idx1,idx2) result(get)
! Wrapper for get_splas_total_p_biomass
!
! @param idx1,idx2  : index of array
! @return get : value of phytoplankton total biomass at time idx1 at layer idx2
! now depending on depth (layers)

  integer,intent(in) :: idx1,idx2
  double precision :: get

  get = simulation%get_splas_total_p_biomass(idx1,idx2)

end function get_phyto_biomass

!------------------------------------------------------------------------------!

function get_phyto_biomass_array(dm1,dm2) result(get)
! Wrapper for get_splas_total_p_biomass that returns the whole array
!
! @param dm1,dm2  : dimension of array
! @return get : array with the values of phytoplankton total biomass at each time
! now depending on depth (layers)

  integer :: idx1,idx2
  integer, intent(in) :: dm1,dm2
  double precision, dimension(dm1,dm2) :: get

  do idx2 = 1,dm2
    do idx1 = 1,dm1
      get(idx1,idx2) = simulation%get_splas_total_p_biomass(idx1,idx2)
    end do
  end do

end function get_phyto_biomass_array

!------------------------------------------------------------------------------!

subroutine set_zoo_biomass(input_a,dim1,dim2)
! Wrapper for set_splas_total_z_biomass
!
! @param input_a : array of zooplankton total biomass in time
! @param dim1,dim2 : size of zooplankton total biomass array
! now depending on depth (layers)

  implicit none

  integer,intent(in) :: dim1,dim2
  double precision,intent(inout),dimension(dim1,dim2) :: input_a

  call simulation%set_splas_total_z_biomass(input_a,dim1,dim2)

end subroutine set_zoo_biomass

!------------------------------------------------------------------------------!

function get_zoo_biomass(idx1,idx2) result(get)
! Wrapper for get_splas_total_z_biomass
!
! @param idx1,idx2  : index of array
! @return get : value of zooplankton total biomass at time idx1 at layer idx2
! now depending on depth (layers)

  integer,intent(in) :: idx1,idx2
  double precision :: get

  get = simulation%get_splas_total_z_biomass(idx1,idx2)

end function get_zoo_biomass

!------------------------------------------------------------------------------!

function get_zoo_biomass_array(dm1,dm2) result(get)
! Wrapper for get_splas_total_z_biomass that returns the whole array
!
! @param dm1,dm2  : dimension of array
! @return get : array with the values of zooplankton total biomass at each time
! now depending on depth (layers)

  integer :: idx1,idx2
  integer, intent(in) :: dm1,dm2
  double precision, dimension(dm1,dm2) :: get

  do idx2 = 1,dm2
    do idx1 = 1,dm1
      get(idx1,idx2) = simulation%get_splas_total_z_biomass(idx1,idx2)
    end do
  end do

end function get_zoo_biomass_array

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

function get_phyto_small_biomass(idx) result(get)
! Wrapper for get_splas_phyto_small_biomass
!
! @return get : biomass of small size phytoplankton
! now dependent on depth (number of layers)

  integer,intent(in) :: idx
  double precision :: get

  get = simulation%get_splas_phyto_small_biomass(idx)

end function get_phyto_small_biomass

!------------------------------------------------------------------------------!

function get_phyto_small_biomass_array(dm) result(get)
! Wrapper for get_splas_phyto_small_biomass for the whole array (all layers)
!
! @return get : biomass of small size phytoplankton
! now dependent on depth (number of layers)

  integer :: idx
  integer,intent(in) :: dm
  double precision, dimension(dm) :: get

  do idx = 1,dm
  get(idx) = simulation%get_splas_phyto_small_biomass(idx)
  end do

end function get_phyto_small_biomass_array

!------------------------------------------------------------------------------!

function get_phyto_total_mean_biomass(idx) result(get)
! Wrapper for get_splas_phyto_total_mean_biomass
!
! @return get : mean total biomass of phytoplankton
! now dependent on depth (number of layers)

  integer, intent(in) :: idx
  double precision :: get

  get = simulation%get_splas_phyto_total_mean_biomass(idx)

end function get_phyto_total_mean_biomass

!------------------------------------------------------------------------------!

function get_phyto_total_mean_biomass_array(dm) result(get)
! Wrapper for get_splas_phyto_total_mean_biomass for the whole array (all layers)
!
! @return get : mean total biomass of phytoplankton
! now dependent on depth (number of layers)

  integer :: idx
  integer,intent(in) :: dm
  double precision, dimension(dm) :: get

  do idx = 1,dm
  get(idx) = simulation%get_splas_phyto_total_mean_biomass(idx)
  end do

end function get_phyto_total_mean_biomass_array

!------------------------------------------------------------------------------!

function get_zoo_total_mean_biomass(idx) result(get)
! Wrapper for get_splas_zoo_total_mean_biomass
!
! @return get : mean total biomass of zooplankton
! now dependent on depth (number of layers)

  integer, intent(in) :: idx
  double precision :: get

  get = simulation%get_splas_zoo_total_mean_biomass(idx)

end function get_zoo_total_mean_biomass

!------------------------------------------------------------------------------!

function get_zoo_total_mean_biomass_array(dm) result(get)
! Wrapper for get_splas_zoo_total_mean_biomass for the whole array (all layers)
!
! @return get : mean total biomass of zooplankton
! now dependent on depth (number of layers)

  integer :: idx
  integer,intent(in) :: dm
  double precision, dimension(dm) :: get

  do idx = 1,dm
  get(idx) = simulation%get_splas_zoo_total_mean_biomass(idx)
  end do

end function get_zoo_total_mean_biomass_array

!------------------------------------------------------------------------------!

function get_phyto_mean_shannon(idx) result(get)
! Wrapper for get_splas_mean_shannon_p
!
! @return get : time mean phytoplankton shannon evenness
! now dependent on depth (number of layers)

  integer, intent(in) :: idx
  double precision :: get

  get = simulation%get_splas_mean_shannon_p(idx)

end function get_phyto_mean_shannon

!------------------------------------------------------------------------------!

function get_zoo_mean_shannon(idx) result(get)
! Wrapper for get_splas_mean_shannon_z
!
! @return get : time mean zooplankton shannon evenness
! now dependent on depth (number of layers)

  integer, intent(in) :: idx
  double precision :: get

  get = simulation%get_splas_mean_shannon_z(idx)

end function get_zoo_mean_shannon

!------------------------------------------------------------------------------!

function get_phyto_mean_shannon_array(dm) result(get)
! Wrapper for get_splas_mean_shannon_p for the whole array (all layers)
!
! @return get : time mean phytoplankton shannon evenness
! now dependent on depth (number of layers)

  integer :: idx
  integer,intent(in) :: dm
  double precision, dimension(dm) :: get

  do idx = 1,dm
  get(idx) = simulation%get_splas_mean_shannon_p(idx)
  end do

end function get_phyto_mean_shannon_array

!------------------------------------------------------------------------------!

function get_zoo_mean_shannon_array(dm) result(get)
! Wrapper for get_splas_mean_shannon_z for the whole array (all layers)
!
! @return get : time mean zooplankton shannon evenness
! now dependent on depth (number of layers)

  integer :: idx
  integer,intent(in) :: dm
  double precision, dimension(dm) :: get

  do idx = 1,dm
  get(idx) = simulation%get_splas_mean_shannon_z(idx)
  end do

end function get_zoo_mean_shannon_array

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  CALCULATION OF PRELIMINARY QUANTITIES WRAPPER
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------------!

subroutine calculate_q()
! Wrapper for subroutine calculation, that computes several quantities
!
!

  call simulation%calculation()

end subroutine calculate_q

!------------------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! WRAPPERS FOR QUANTITIES COMPUTED BY calculate_q()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------------!

function get_phytoplankton_size_xp(idx) result(get)
! Wrapper for get_splas_phyto_xp
!
! @param idx  : index of array
! @return get : value of phytoplankton size (xp) at idx
!
  integer,intent(in) :: idx
  double precision :: get

  get = simulation%get_splas_phyto_xp(idx)

end function get_phytoplankton_size_xp

!------------------------------------------------------------------------------!

function get_phytoplankton_size_xp_array(dm) result(get)
! Wrapper for get_splas_phyto_xp that returns the whole array
!
! @param dm  : dimension of array
! @return get : array with the values of phytoplankton sizes (xp)
!
  integer :: idx
  integer, intent(in) :: dm
  double precision, dimension(dm) :: get

  do idx = 1,dm
  get(idx) = simulation%get_splas_phyto_xp(idx)
  end do

end function get_phytoplankton_size_xp_array

!------------------------------------------------------------------------------!

function get_zooplankton_size_xz(idx) result(get)
! Wrapper for get_splas_zoo_xz
!
! @param idx  : index of array
! @return get : value of zooplankton size (xz) at idx
!
  integer,intent(in) :: idx
  double precision :: get

  get = simulation%get_splas_zoo_xz(idx)

end function get_zooplankton_size_xz

!------------------------------------------------------------------------------!

function get_zooplankton_size_xz_array(dm) result(get)
! Wrapper for get_splas_zoo_xz that returns the whole array
!
! @param dm  : dimension of array
! @return get : array with the values of zooplankton sizes (xz)
!
  integer :: idx
  integer, intent(in) :: dm
  double precision, dimension(dm) :: get

  do idx = 1,dm
  get(idx) = simulation%get_splas_zoo_xz(idx)
  end do

end function get_zooplankton_size_xz_array

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

function get_zp_optimum_prey_size(idx) result(get)
! Wrapper for get_splas_zp_xoptz
!
! @param idx  : index of array
! @return get : value of optimum phytoplankton size for zooplankton diet at idx
!
  integer,intent(in) :: idx
  double precision :: get

  get = simulation%get_splas_zp_xoptz(idx)

end function get_zp_optimum_prey_size

!------------------------------------------------------------------------------!

function get_zp_optimum_prey_size_array(dm) result(get)
! Wrapper for get_splas_zp_xoptz that returns the whole array
!
! @param dm  : dimension of array
! @return get : array with the values of phytoplankton optimum size to be eaten
!               by zooplankton
!
  integer :: idx
  integer, intent(in) :: dm
  double precision, dimension(dm) :: get

  do idx = 1,dm
  get(idx) = simulation%get_splas_zp_xoptz(idx)
  end do

end function get_zp_optimum_prey_size_array

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

function get_phytoplankton_growth_rate(idx) result(get)
! Wrapper for get_splas_phyto_mu_0
!
! @param idx  : index of array
! @return get : value of phytoplankton growth rate at idx
!
  integer,intent(in) :: idx
  double precision :: get

  get = simulation%get_splas_phyto_mu_0(idx)

end function get_phytoplankton_growth_rate

!------------------------------------------------------------------------------!

function get_phytoplankton_growth_rate_array(dm) result(get)
! Wrapper for get_splas_phyto_mu_0 that returns the whole array
!
! @param dm  : dimension of array
! @return get : array with the values of phytoplankton growth rate
!
  integer :: idx
  integer, intent(in) :: dm
  double precision, dimension(dm) :: get

  do idx = 1,dm
  get(idx) = simulation%get_splas_phyto_mu_0(idx)
  end do

end function get_phytoplankton_growth_rate_array

!------------------------------------------------------------------------------!

function get_zooplankton_ingestion_rate(idx) result(get)
! Wrapper for get_splas_zoo_irz
!
! @param idx  : index of array
! @return get : value of zooplankton ingestion rate at idx
!
  integer,intent(in) :: idx
  double precision :: get

  get = simulation%get_splas_zoo_irz(idx)

end function get_zooplankton_ingestion_rate

!------------------------------------------------------------------------------!

function get_zooplankton_ingestion_rate_array(dm) result(get)
! Wrapper for get_splas_zoo_irz that returns the whole array
!
! @param dm  : dimension of array
! @return get : array with the values of zooplankton ingestion rate
!
  integer :: idx
  integer, intent(in) :: dm
  double precision, dimension(dm) :: get

  do idx = 1,dm
  get(idx) = simulation%get_splas_zoo_irz(idx)
  end do

end function get_zooplankton_ingestion_rate_array

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

function get_phytoplankton_half_sat(idx) result(get)
! Wrapper for get_splas_phyto_ks
!
! @param idx  : index of array
! @return get : value of phytoplankton half saturation for nut at idx
!
  integer,intent(in) :: idx
  double precision :: get

  get = simulation%get_splas_phyto_ks(idx)

end function get_phytoplankton_half_sat

!------------------------------------------------------------------------------!

function get_phytoplankton_half_sat_array(dm) result(get)
! Wrapper for get_splas_phyto_ks that returns the whole array
!
! @param dm  : dimension of array
! @return get : array with the values of phytoplankton half saturation for nut
!
  integer :: idx
  integer, intent(in) :: dm
  double precision, dimension(dm) :: get

  do idx = 1,dm
  get(idx) = simulation%get_splas_phyto_ks(idx)
  end do

end function get_phytoplankton_half_sat_array

!------------------------------------------------------------------------------!

function get_zp_pref_matrix_elem(idx1,idx2) result(get)
! Wrapper for get_splas_zp_preference_phi, to obtain one element i,j
!
! @param idx1,idx2  : index of array
! @return get : value of preference matrix at row idx1 and column idx2
!
  integer,intent(in) :: idx1,idx2
  double precision :: get

  get = simulation%get_splas_zp_preference_phi(idx1,idx2)

end function get_zp_pref_matrix_elem

!------------------------------------------------------------------------------!

function get_zp_preference_matrix(dm1,dm2) result(get)
! Wrapper for get_splas_zp_preference_phi that returns the whole array (matrix)
!
! @param dm1,dm2 : dimension of array (should be square matrix) dm1=dm2 ...
! @return get : matrix with the values of relative preference z-p
!
  integer :: idx1,idx2
  integer, intent(in) :: dm1,dm2
  double precision, dimension(dm1,dm2) :: get

  do idx2 = 1,dm2
    do idx1 = 1,dm1
      get(idx1,idx2) = simulation%get_splas_zp_preference_phi(idx1,idx2)
    end do
  end do

end function get_zp_preference_matrix

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!!------------------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  SIMULATION SUBROUTINE WRAPPER
!
!  subroutine that integrates the differential equations calling
!  integration() from splas.f90
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------------!

subroutine integrate_eqs()
! Wrapper for integration() that integrates all the differential equations
!
!
!

  call simulation%integration()

end subroutine integrate_eqs

!------------------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
! GET FUNCTIONS WRAPPER FOR TIME-DEPENDENT QUANTITIES NEEDED IN THE EQUATIONS
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------------!

function get_nutrient_uptake_p(idx1,idx2) result(get)
! Wrapper for get_splas_nutrient_uptake
!
! @param idx1,idx2  : index of array
! @return get : value of nutrient uptake by phytoplankton at idx1,idx2
! now depending on depth (layers)

  integer,intent(in) :: idx1,idx2
  double precision :: get

  get = simulation%get_splas_nutrient_uptake(idx1,idx2)

end function get_nutrient_uptake_p

!------------------------------------------------------------------------------!

function get_nutrient_uptake_p_array(dm1,dm2) result(get)
! Wrapper for get_splas_nutrient_uptake that returns the whole array
!
! @param dm1,dm2  : dimension of array
! @return get : array with the values of nutrient uptake by phytoplankton
! now depending on depth (layers)

  integer :: idx1,idx2
  integer, intent(in) :: dm1,dm2
  double precision, dimension(dm1,dm2) :: get

  do idx2 = 1,dm2
    do idx1 = 1,dm1
      get(idx1,idx2) = simulation%get_splas_nutrient_uptake(idx1,idx2)
    end do
  end do

end function get_nutrient_uptake_p_array

!------------------------------------------------------------------------------!

function get_p_growth_rate(idx1,idx2) result(get)
! Wrapper for get_splas_p_growth_rate
!
! @param idx1,idx2  : index of array
! @return get : value of phytoplankton type idx at layer idx2 growth rate (quota model)
! now depending on depth (layers)

  integer,intent(in) :: idx1,idx2
  double precision :: get

  get = simulation%get_splas_p_growth_rate(idx1,idx2)

end function get_p_growth_rate

!------------------------------------------------------------------------------!

function get_p_growth_rate_array(dm1,dm2) result(get)
! Wrapper for get_splas_p_growth_rate that returns the whole array
!
! @param dm1,dm2 : dimension of array
! @return get : array with the values of phytoplankton growth rate (quota model)
! now depending on depth (layers)

  integer :: idx1,idx2
  integer, intent(in) :: dm1,dm2
  double precision, dimension(dm1,dm2) :: get

  do idx2 = 1,dm2
    do idx1 = 1,dm1
      get(idx1,idx2) = simulation%get_splas_p_growth_rate(idx1,idx2)
    end do
  end do

end function get_p_growth_rate_array

!------------------------------------------------------------------------------!

function get_grazing_matrix_elem(idx1,idx2,idx3) result(get)
! Wrapper for get_splas_grazing_z_p, to obtain one element i,j,k
!
! @param idx1,idx2,idx3  : index of array
! @return get : value of grazing matrix at row idx1 and column idx2 and layer idx3
! now depending on depth (layers)

  integer,intent(in) :: idx1,idx2,idx3
  double precision :: get

  get = simulation%get_splas_grazing_z_p(idx1,idx2,idx3)

end function get_grazing_matrix_elem

!------------------------------------------------------------------------------!

function get_grazing_matrix(dm1,dm2,dm3) result(get)
! Wrapper for get_splas_grazing_z_p that returns the whole array (matrix)
!
! @param dm1,dm2,dm3 : dimension of array (should be square matrix) dm1=dm2, dm3 anything
! @return get : matrix with the values of z-p grazing at each layer
!
  integer :: idx1,idx2,idx3
  integer, intent(in) :: dm1,dm2,dm3
  double precision, dimension(dm1,dm2,dm3) :: get

  do idx3 = 1,dm3
    do idx2 = 1,dm2
      do idx1 = 1,dm1
        get(idx1,idx2,idx3) = simulation%get_splas_grazing_z_p(idx1,idx2,idx3)
      end do
    end do
  end do

end function get_grazing_matrix

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

function get_total_p_predation(idx1,idx2) result(get)
! Wrapper for get_splas_total_p_predation_z
!
! @param idx1,idx2  : index of array
! @return get : value of total predation on all phytoplankton by nano-micro zoo idx1,idx2
! now depending on depth (layers)

  integer,intent(in) :: idx1,idx2
  double precision :: get

  get = simulation%get_splas_total_p_predation_z(idx1,idx2)

end function get_total_p_predation

!------------------------------------------------------------------------------!

function get_total_p_predation_array(dm1,dm2) result(get)
! Wrapper for get_splas_total_p_predation_z that returns the whole array
!
! @param dm1,dm2 : dimension of array
! @return get : array with the values of predation on all phyto classes by nano-micro zoo
! now depending on depth (layers)

  integer :: idx1,idx2
  integer, intent(in) :: dm1,dm2
  double precision, dimension(dm1,dm2) :: get

  do idx2 = 1,dm2
    do idx1 = 1,dm1
      get(idx1,idx2) = simulation%get_splas_total_p_predation_z(idx1,idx2)
    end do
  end do

end function get_total_p_predation_array

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------!

!PHYTOPLANKTON

!------------------------------------------------------------------------------!

function get_pcell_volume(idx) result(get)
! Wrapper for get_splas_p_cell_volume
!
! @param idx  : index of array
! @return get : value of phytoplankton idx type cell volume
!
  integer,intent(in) :: idx
  double precision :: get

  get = simulation%get_splas_p_cell_volume(idx)

end function get_pcell_volume

!------------------------------------------------------------------------------!

function get_pcell_volume_array(dm) result(get)
! Wrapper for get_splas_p_cell_volume that returns the whole array
!
! @param dm  : dimension of array
! @return get : array with the values of phytoplankton cell volume
!
  integer :: idx
  integer, intent(in) :: dm
  double precision, dimension(dm) :: get

  do idx = 1,dm
  get(idx) = simulation%get_splas_p_cell_volume(idx)
  end do

end function get_pcell_volume_array

!------------------------------------------------------------------------------!

function get_pv_uptake(idx) result(get)
! Wrapper for get_splas_p_v_uptake
!
! @param idx  : index of array
! @return get : value of phytoplankton idx type maximum uptake rate
!
  integer,intent(in) :: idx
  double precision :: get

  get = simulation%get_splas_p_v_uptake(idx)

end function get_pv_uptake

!------------------------------------------------------------------------------!

function get_pv_uptake_array(dm) result(get)
! Wrapper for get_splas_p_v_uptake that returns the whole array
!
! @param dm  : dimension of array
! @return get : array with the values of phytoplankton maximum uptake rate
!
  integer :: idx
  integer, intent(in) :: dm
  double precision, dimension(dm) :: get

  do idx = 1,dm
  get(idx) = simulation%get_splas_p_v_uptake(idx)
  end do

end function get_pv_uptake_array

!------------------------------------------------------------------------------!

function get_pmin_cell_quota(idx) result(get)
! Wrapper for get_splas_p_min_cell_quota
!
! @param idx  : index of array
! @return get : value of phytoplankton idx type minimum cell quota
!
  integer,intent(in) :: idx
  double precision :: get

  get = simulation%get_splas_p_min_cell_quota(idx)

end function get_pmin_cell_quota

!------------------------------------------------------------------------------!

function get_pmin_cell_quota_array(dm) result(get)
! Wrapper for get_splas_p_min_cell_quota that returns the whole array
!
! @param dm  : dimension of array
! @return get : array with the values of phytoplankton minimum cell quota
!
  integer :: idx
  integer, intent(in) :: dm
  double precision, dimension(dm) :: get

  do idx = 1,dm
  get(idx) = simulation%get_splas_p_min_cell_quota(idx)
  end do

end function get_pmin_cell_quota_array

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!! CARBON SECTION
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------------!

function get_p_carbon(idx) result(get)
! Wrapper for get_splas_phyto_carbon
!
! @param idx  : index of array
! @return get : value of phytoplankton idx type carbon content
!
  integer,intent(in) :: idx
  double precision :: get

  get = simulation%get_splas_phyto_carbon(idx)

end function get_p_carbon

!------------------------------------------------------------------------------!

function get_p_carbon_array(dm) result(get)
! Wrapper for get_splas_phyto_carbon that returns the whole array
!
! @param dm  : dimension of array
! @return get : array with the values of phytoplankton carbon content
!
  integer :: idx
  integer, intent(in) :: dm
  double precision, dimension(dm) :: get

  do idx = 1,dm
  get(idx) = simulation%get_splas_phyto_carbon(idx)
  end do

end function get_p_carbon_array

!------------------------------------------------------------------------------!

subroutine set_p_cexport(input_a,dim1,dim2)
! Wrapper for set_splas_phyto_c_export
!
! @param input_a : array of phytoplankton c export
! @param dim1,dim2 : size of phytoplankton c export array
!
  implicit none

  integer,intent(in) :: dim1,dim2
  double precision,intent(inout),dimension(dim1,dim2) :: input_a

  call simulation%set_splas_phyto_c_export(input_a,dim1,dim2)

end subroutine set_p_cexport

!------------------------------------------------------------------------------!

function get_p_cexport(idx1,idx2) result(get)
! Wrapper for get_splas_phyto_c_export
!
! @param idx1,idx2  : index of array
! @return get : value of phytoplankton c export at idx1 at time idx2
!
  integer,intent(in) :: idx1,idx2
  double precision :: get
  get = simulation%get_splas_phyto_c_export(idx1,idx2)

end function get_p_cexport

!------------------------------------------------------------------------------!

function get_p_cexport_array(dm,idx2) result(get)
! Wrapper for get_splas_phyto_c_export that returns the whole array at a certain time
!
! @param dm  : dimension of array
! @param idx2 : time index of array
! @return get : array with the values of phytoplankton c export
!
  integer :: idx
  integer, intent(in) :: dm,idx2
  double precision, dimension(dm) :: get

  do idx = 1,dm
  get(idx) = simulation%get_splas_phyto_c_export(idx,idx2)
  end do

end function get_p_cexport_array

!------------------------------------------------------------------------------!

function get_p_cexport_time(idx1,dm2) result(get)
! Wrapper for get_splas_phyto_c_export that returns one phyto c exp class at all times
!
! @param dm2 : dimension of array
! @param idx1 : index of phytoplankton size class
! @return get : array with the values of phyto c export class idx1 at all times
!
  integer :: idx
  integer, intent(in) :: dm2,idx1
  double precision, dimension(dm2) :: get

  do idx = 1,dm2
  get(idx) = simulation%get_splas_phyto_c_export(idx1,idx)
  end do

end function get_p_cexport_time

!------------------------------------------------------------------------------!

function get_whole_p_cexport(dm1,dm2) result(get)
! Wrapper for get_splas_phyto_c_export that returns the whole array
!
! @param dm1,dm2  : dimension of array
! @return get : array with the values of phytoplankton c export at all times
!
  integer :: idx1,idx2
  integer, intent(in) :: dm1,dm2
  double precision, dimension(dm1,dm2) :: get

  do idx2 = 1,dm2
    do idx1 = 1,dm1
  get(idx1,idx2) = simulation%get_splas_phyto_c_export(idx1,idx2)
    end do
  end do

end function get_whole_p_cexport

!------------------------------------------------------------------------------!

subroutine set_p_ctotalexport(input_a,dim1)
! Wrapper for set_splas_phyto_c_total_export
!
! @param input_a : array of phytoplankton c total export
! @param dim1 : size of phytoplankton c total export array
!
  implicit none

  integer,intent(in) :: dim1
  double precision,intent(inout),dimension(dim1) :: input_a

  call simulation%set_splas_phyto_c_total_export(input_a,dim1)

end subroutine set_p_ctotalexport

!------------------------------------------------------------------------------!

function get_p_ctotalexport(idx1) result(get)
! Wrapper for get_splas_phyto_c_total_export
!
! @param idx1 : index of array
! @return get : value of phytoplankton c total export at time idx1
!
  integer,intent(in) :: idx1
  double precision :: get
  get = simulation%get_splas_phyto_c_total_export(idx1)

end function get_p_ctotalexport

!------------------------------------------------------------------------------!

function get_p_ctotalexport_array(dm) result(get)
! Wrapper for get_splas_phyto_c_total_export that returns the whole time array
!
! @param dm  : dimension of array
! @return get : array with the values of phytoplankton c total export
!
  integer :: idx
  integer, intent(in) :: dm
  double precision, dimension(dm) :: get

  do idx = 1,dm
  get(idx) = simulation%get_splas_phyto_c_total_export(idx)
  end do

end function get_p_ctotalexport_array

!------------------------------------------------------------------------------!

subroutine set_z_cexport(input_a,dim1)
! Wrapper for set_splas_zoo_c_export
!
! @param input_a : array of zooplankton c total export
! @param dim1 : size of zooplankton c total export array
!
  implicit none

  integer,intent(in) :: dim1
  double precision,intent(inout),dimension(dim1) :: input_a

  call simulation%set_splas_zoo_c_export(input_a,dim1)

end subroutine set_z_cexport

!------------------------------------------------------------------------------!

function get_z_cexport(idx1) result(get)
! Wrapper for get_splas_zoo_c_export
!
! @param idx1 : index of array
! @return get : value of zooplankton c total export at time idx1
!
  integer,intent(in) :: idx1
  double precision :: get
  get = simulation%get_splas_zoo_c_export(idx1)

end function get_z_cexport

!------------------------------------------------------------------------------!

function get_z_cexport_array(dm) result(get)
! Wrapper for get_splas_zoo_c_export that returns the whole time array
!
! @param dm  : dimension of array
! @return get : array with the values of zooplankton c total export
!
  integer :: idx
  integer, intent(in) :: dm
  double precision, dimension(dm) :: get

  do idx = 1,dm
  get(idx) = simulation%get_splas_zoo_c_export(idx)
  end do

end function get_z_cexport_array

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!


end module wrapper_py
