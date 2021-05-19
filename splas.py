from splaspy import wrapper_py as splaspyw
import splasplot as spl
import numpy as np
import splasIO as iojson
import sys
import os
from sys import argv
import pickle

dir_pickle = "./SPLAS_PICKLE/"

class splas:

    def __init__(self,file):
        splaspyw.init_wrapper()
        self.io = iojson.io_json(file)
        self.io.__init__(file)

    def get_input(self):
        self.io.read()

    def sim(self):
        print("Start simulation")
        splaspyw.integrate_eqs()
        print("End simulation")

    def write(self):
        splaspyw.set_simulation_length(self.io.simulation_length)
        splaspyw.set_time_step(self.io.time_step)
        splaspyw.set_nt()
        splaspyw.set_print_step(self.io.print_step)

        ################################################################################################################################################
        self.nt = splaspyw.get_nt()
        self.ps = splaspyw.get_print_step()
        self.tdim = int(self.nt/self.ps)+1
        self.days = np.arange(0,self.io.simulation_length,self.ps*self.io.time_step)
        #IRRADIANCE
        
        self.itemp = np.asfortranarray(np.linspace(self.io.i,self.io.i + (self.io.delta),self.nt+1)) # increasing array version 
        self.izero = np.asfortranarray(np.random.normal(self.itemp,self.io.sd,(self.nt+1)))  # random version
        splaspyw.set_i_val(self.izero,(self.nt+1))
        #self.itemp = np.asfortranarray(np.genfromtxt('Irr_1h.csv',delimiter=','))
        #self.izero = splaspyw.set_i_val(self.itemp,(self.nt+1))
        ################################################################################################################################################

        #vertical space-grid
        ################################################################################################################################################
        splaspyw.set_delta_z(self.io.dz)
        splaspyw.set_number_layers(self.io.m)
        splaspyw.set_initial_diff(self.io.init_diff)
        splaspyw.set_diffusion(np.zeros(self.io.m))

        self.h = self.io.dz*float(self.io.m - 1)
        self.z = np.zeros(self.io.m)
        for k in range(self.io.m):
            self.z[k] = (float(k+1)-0.5)*self.io.dz
        ################################################################################################################################################

        splaspyw.set_nutrient_number(self.io.nut_number)
        splaspyw.set_phytoplankton_dim(self.io.phyto_num_class)
        splaspyw.set_zooplankton_dim(self.io.zoo_num_class)
        #1-D
        splaspyw.set_nutrients(np.zeros((self.io.nut_number,self.tdim,self.io.m), order='F'))
        splaspyw.set_phytoplankton(np.zeros((self.io.phyto_num_class,self.tdim,self.io.m), order='F'))
        splaspyw.set_zooplankton(np.zeros((self.io.zoo_num_class,self.tdim,self.io.m), order='F'))
        splaspyw.set_pquota(np.zeros((self.io.phyto_num_class,self.tdim,self.io.m), order='F'))
	#IRRADIANCE            fil
        splaspyw.set_irradiance_profile(np.zeros((self.tdim,self.io.m), order='F'))
        #DETRITUS (for now time x depth dimensions)
        splaspyw.set_detritus(np.zeros((self.tdim,self.io.m), order='F'))

        splaspyw.set_nutrient_metrics(np.zeros((self.io.phyto_num_class,self.tdim,self.io.m), order='F'))
        splaspyw.set_grazing_metrics(np.zeros((self.io.phyto_num_class,self.tdim,self.io.m), order='F'))

        splaspyw.set_phyto_time_mean(np.zeros((self.io.phyto_num_class,self.io.m), order='F'))
        splaspyw.set_zoo_time_mean(np.zeros((self.io.zoo_num_class,self.io.m), order='F'))
        splaspyw.set_phyto_shannon_index(np.zeros((self.tdim,self.io.m), order='F'))
        splaspyw.set_zoo_shannon_index(np.zeros((self.tdim,self.io.m), order='F'))
        splaspyw.set_phyto_biomass(np.zeros((self.tdim,self.io.m), order='F'))
        splaspyw.set_zoo_biomass(np.zeros((self.tdim,self.io.m), order='F'))
        #
        splaspyw.set_nutrient_supply(self.io.nut_supply)
        splaspyw.set_time_variation(bool(self.io.time_flag))
        #PHYTOPLANKTON
        splaspyw.set_phytoplankton_xmin(self.io.phyto_min_size)
        splaspyw.set_phytoplankton_xmax(self.io.phyto_max_size)
        splaspyw.set_phytoplankton_mort(self.io.phyto_mortality)
        splaspyw.set_phyto_sinking(self.io.phyto_sink)
        #ZOOPLANKTON
        splaspyw.set_zooplankton_xmin(self.io.zoo_min_size)
        splaspyw.set_zooplankton_xmax(self.io.zoo_max_size)
        splaspyw.set_zooplankton_growth(self.io.zoo_growth)
        splaspyw.set_zooplankton_egestion(self.io.zoo_egestion)
        splaspyw.set_zoooplankton_mort(self.io.zoo_mortality)
        splaspyw.set_zoooplankton_half_sat(self.io.zoo_hs)
        splaspyw.set_zoooplankton_size_tol(self.io.zoo_dx)

        #initial values
        splaspyw.set_initial_n(self.io.nuts_value)
        splaspyw.set_initial_p(self.io.init_p)
        splaspyw.set_initial_z(self.io.init_z)
        splaspyw.set_initial_pquota(self.io.init_pquota)

        #carbon
        splaspyw.set_p_cexport(np.zeros((self.io.phyto_num_class,self.tdim), order='F'))
        splaspyw.set_p_ctotalexport(np.zeros(self.tdim))
        splaspyw.set_z_cexport(np.zeros(self.tdim))

        #calculation
        splaspyw.calculate_q()

################################################################################################################################################

    def store(self):

        self.p_ks = splaspyw.get_phytoplankton_half_sat_array(self.io.phyto_num_class)
        self.p_mu = splaspyw.get_phytoplankton_growth_rate_array(self.io.phyto_num_class)
        self.p_mort = splaspyw.get_phytoplankton_mort()
        #new!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        self.p_vmax = splaspyw.get_pv_uptake_array(self.io.phyto_num_class)
        self.p_qmin = splaspyw.get_pmin_cell_quota_array(self.io.phyto_num_class)
        self.p_cell_vol = splaspyw.get_pcell_volume_array(self.io.phyto_num_class)
        self.psink = splaspyw.get_phyto_sinking()
        #!!!!!!!!!!!!!!!!!!!!!!!!!!
        self.p_naff = self.p_vmax/self.p_ks
        self.p_scaled_naff = self.p_naff/self.p_qmin

        self.z_ir = splaspyw.get_zooplankton_ingestion_rate_array(self.io.zoo_num_class)
        self.z_xopt = splaspyw.get_zp_optimum_prey_size_array(self.io.zoo_num_class)
        self.z_dx = splaspyw.get_zooplankton_size_tol()

        self.xp = splaspyw.get_phytoplankton_size_xp_array(self.io.phyto_num_class)
        self.xz = splaspyw.get_zooplankton_size_xz_array(self.io.zoo_num_class)

        #size-dependent sinking of phytoplankton and nutrient requirement R*
        self.phyto_sink = self.psink*(self.xp**1.17)
        #self.p_r_star = self.p_ks*((self.p_mort+self.phyto_sink)/(self.p_mu-self.p_mort-self.phyto_sink))
        #self.p_r_star = (self.p_mu*self.p_qmin*self.p_ks*(self.p_mort+self.phyto_sink))/(self.p_vmax*self.p_mu-(self.p_vmax+self.p_mu*self.p_qmin)*(self.p_mort+self.phyto_sink))

        #estimate of r_star with sinking and diffusion (qualitative)
        self.p_r_star = (self.p_mu*self.p_qmin*self.p_ks*(self.p_mort+(self.phyto_sink/self.io.dz)-(self.io.init_diff)/(self.io.dz**2)))/(self.p_vmax*self.p_mu-(self.p_vmax+self.p_mu*self.p_qmin)*(self.p_mort+(self.phyto_sink/self.io.dz)-(self.io.init_diff)/(self.io.dz**2)))

        self.phi = splaspyw.get_zp_preference_matrix(self.io.phyto_num_class,self.io.zoo_num_class)

        #######DEPTH DEPENDENT QUANTITIES###################################################################################################################

        #whole arrays (3) size-time-space
        self.phyto = splaspyw.get_whole_phyto_array(self.io.phyto_num_class,self.tdim,self.io.m)
        self.zoo = splaspyw.get_whole_zoo_array(self.io.zoo_num_class,self.tdim,self.io.m)
        self.pquota = splaspyw.get_whole_pquota(self.io.phyto_num_class,self.tdim,self.io.m)
        self.nutrients = splaspyw.get_whole_nutrients_array(self.io.nut_number,self.tdim,self.io.m)
        self.pncells = (self.phyto/self.pquota)
        self.irrprofile = splaspyw.get_irradiance_profile(self.tdim,self.io.m)
        #final time arrays (2) size-space
        self.fphyto = self.phyto[:,-1,:]
        self.fzoo = self.zoo[:,-1,:]
        self.fpquota = self.pquota[:,-1,:]
        self.fpncells = (self.fphyto/self.fpquota)
        #each nutrient at all times and layers (2) time-space
        self.Nitrogen = splaspyw.get_1_nut_all_times(1,self.tdim,self.io.m)
        self.Phosphorus = splaspyw.get_1_nut_all_times(2,self.tdim,self.io.m)
        self.Carbon = splaspyw.get_1_nut_all_times(3,self.tdim,self.io.m)
        #mean and final concentration of Nitrogen (1) at all layers -space
        self.mean_N = (1.0/self.io.simulation_length)*self.Nitrogen.sum(axis=0)*self.io.time_step*self.io.print_step
        self.final_N = self.Nitrogen[-1,:]

        #DETRITUS
        #whole array time x depth
        self.detritus = splaspyw.get_whole_detritus(self.tdim,self.io.m)
        #at final time (only depth dependence)
        self.final_detritus = self.detritus[-1,:]
        #time average (only depth dependence)
        self.mean_detritus = (1.0/self.io.simulation_length)*self.detritus.sum(axis=0)*self.io.time_step*self.io.print_step

        #time mean plankton size-spectra (time mean at each size class and each layer), (2) size-space
        self.phyto_mean = splaspyw.get_phyto_time_mean_array(self.io.phyto_num_class,self.io.m)
        self.zoo_mean = splaspyw.get_zoo_time_mean_array(self.io.zoo_num_class,self.io.m)
        #plankton total biomass (summed over size-classes), (2) time-space
        self.phyto_biomass = splaspyw.get_phyto_biomass_array(self.tdim,self.io.m)
        self.zoo_biomass = splaspyw.get_zoo_biomass_array(self.tdim,self.io.m)
        self.z_p_ratio = self.zoo_biomass/self.phyto_biomass
        #time mean z:p (1) -space
        self.ave_z_p_ratio = (1.0/self.io.simulation_length)*self.z_p_ratio.sum(axis=0)*self.io.time_step*self.io.print_step
        #picophytoplankton and total plankton time mean biomass, (1) -space
        self.small_phyto = splaspyw.get_phyto_small_biomass_array(self.io.m)
        self.phyto_total_mean_biomass = splaspyw.get_phyto_total_mean_biomass_array(self.io.m)
        self.zoo_total_mean_biomass = splaspyw.get_zoo_total_mean_biomass_array(self.io.m)

        #############################################
        #limitations whole arrays (3) size-time-space
        self.nutlim = splaspyw.get_whole_nut_metrics(self.io.phyto_num_class,self.tdim,self.io.m)
        self.grazlim = splaspyw.get_whole_graz_metrics(self.io.phyto_num_class,self.tdim,self.io.m)
        #mean and final limitations (2) size-space
        self.aven = (1.0/self.io.simulation_length)*self.nutlim.sum(axis=1)*self.io.time_step*self.io.print_step
        self.aveg = (1.0/self.io.simulation_length)*self.grazlim.sum(axis=1)*self.io.time_step*self.io.print_step
        self.nut_lim_final = self.nutlim[:,-1,:]
        self.graz_lim_final = self.grazlim[:,-1,:]
        ############################################

        #now grazing, uptakeP(N), muP(Q), total Zj predation on P DEPEND ON DEPTH, size-space
        self.grazing_zp = splaspyw.get_grazing_matrix(self.io.phyto_num_class,self.io.zoo_num_class,self.io.m)
        self.nit_uptake = splaspyw.get_nutrient_uptake_p_array(self.io.phyto_num_class,self.io.m)
        self.p_growth_rate = splaspyw.get_p_growth_rate_array(self.io.phyto_num_class,self.io.m)
        self.z_predation = splaspyw.get_total_p_predation_array(self.io.zoo_num_class,self.io.m)

        #SHANNON EVENNESS###############################
        #also shannon evenness is defined per each space layer, time-space
        self.phyto_biodiv = splaspyw.get_phyto_sh_idx_array(self.tdim,self.io.m)
        self.zoo_biodiv =  splaspyw.get_zoo_sh_idx_array(self.tdim,self.io.m)
        #time mean shannon evenness, (1) -space
        self.mean_shannon_p = splaspyw.get_phyto_mean_shannon_array(self.io.m)
        self.mean_shannon_z = splaspyw.get_zoo_mean_shannon_array(self.io.m)
        ################################################

        #initial quantities
        self.inutrients = self.Nitrogen[1,:]
        self.iphyto = self.phyto[:,1,:]
        self.izoo = self.zoo[:,1,:]
        self.ipquota = self.pquota[:,1,:]
        self.ipncells = (self.iphyto/self.ipquota)

        #diffusion depth profile
        self.diffusion_profile = splaspyw.get_diffusion_array(self.io.m)

        #######END OF DEPTH DEPENDENT QUANTITIES###################################################################################################################

        #VERTICAL SPACE AVERAGE
        self.avez_small_phyto = (1/self.h)*self.io.dz*self.small_phyto.sum()
        self.avez_p_biomass = (1/self.h)*self.io.dz*self.phyto_total_mean_biomass.sum()
        self.avez_z_biomass = (1/self.h)*self.io.dz*self.zoo_total_mean_biomass.sum()
        self.avez_shannon_p = (1/self.h)*self.io.dz*self.mean_shannon_p.sum()
        self.avez_shannon_z = (1/self.h)*self.io.dz*self.mean_shannon_z.sum()
        self.avez_N = (1/self.h)*self.io.dz*self.mean_N.sum()
        self.avez_z_p = (1/self.h)*self.io.dz*self.ave_z_p_ratio.sum()

        #####CARBON
        #self.p_carbon = 0.0109*((self.p_cell_vol)**0.991)
        #C export in time (mgC m-2 day-1)
        self.p_carbon = splaspyw.get_p_carbon_array(self.io.phyto_num_class)
        self.p_carbon_export_matrix = splaspyw.get_whole_p_cexport(self.io.phyto_num_class,self.tdim)
        self.p_total_carbon_export = splaspyw.get_p_ctotalexport_array(self.tdim)
        ########calculated from detritus dynamics
        self.z_carbon_export = splaspyw.get_z_cexport_array(self.tdim)
        ########
        self.total_carbon_export = self.p_total_carbon_export + self.z_carbon_export

        ################### PAY ATTENTION TO CALCULATION HERE ####################################################################################
        #REMOVE INITIAL CONDITIONS THAT OFFSET EXPORT
        #average export
        #average C export on whole simulation length (mgC m-2 day-1)
        self.phyto_ave_c_exp = (1.0/self.io.simulation_length)*(self.p_carbon_export_matrix.sum(axis=1) - self.p_carbon_export_matrix[:,0])*self.io.time_step*self.io.print_step
        self.phyto_ave_total_c_exp = (1.0/self.io.simulation_length)*(self.p_total_carbon_export.sum(axis=0) - self.p_total_carbon_export[0])*self.io.time_step*self.io.print_step
        #from detritus
        self.zoo_ave_total_c_exp = (1.0/self.io.simulation_length)*self.z_carbon_export.sum(axis=0)*self.io.time_step*self.io.print_step
        ###
        self.total_ave_carbon_export = (1.0/self.io.simulation_length)*(self.total_carbon_export.sum(axis=0) - self.total_carbon_export[0])*self.io.time_step*self.io.print_step
        #integrated export
        #total integrated export during simulation length (mgC m-2)
        self.phyto_integrated_c_exp = (self.p_carbon_export_matrix.sum(axis=1) - self.p_carbon_export_matrix[:,0])*self.io.time_step*self.io.print_step
        self.phyto_integrated_total_c_exp = (self.p_total_carbon_export.sum(axis=0) - self.p_total_carbon_export[0])*self.io.time_step*self.io.print_step
        #from detritus
        self.zoo_integrated_total_c_exp = self.z_carbon_export.sum(axis=0)*self.io.time_step*self.io.print_step
        ###
        self.total_integrated_carbon_export = (self.total_carbon_export.sum(axis=0) - self.total_carbon_export[0])*self.io.time_step*self.io.print_step
        ##########################################################################################################################################

        #save pickle files in directory dir_pickle
        self.write_on_file(dir_pickle)

################################################################################################################################################

    def get_simulation_output(self):

        print("   ")
        print(" NUMBER OF TIME STEPS AND VALUE OF PRINTING STEP IN THE SIMULATION ")
        print(self.nt,"  ",self.ps,"  total time dimension = ",self.tdim)

        print("   ")
        print(" PARAMETERS OF THE VERTICAL GRID (SPLAS 1-D) ")
        print(" step-size = ",self.io.dz," #layers = ",self.io.m," total depth = ",self.h)

        print("   ")
        print(" POINTS IN THE GRID ")
        print(self.z)

        print(" DIMENSION OF N,P,Z ")
        print(splaspyw.get_nutrient_number())
        print(splaspyw.get_phytoplankton_dim())
        print(splaspyw.get_zooplankton_dim())
        print(" 1ST COMPONENT OF N,P,Z ")
        print(splaspyw.get_nutrients_x(1,self.tdim,1))
        print(splaspyw.get_phytoplankton_x(1,self.tdim,1))
        print(splaspyw.get_zooplankton_x(1,self.tdim,1))
        print("   ")

        #if sys.flags.interactive:
        #os._exit(1)

       #######################################################################

        spl.get_parameter_plot(self)

        #layer-dependent
        spl.get_final_plot(self)
        spl.get_time_plot(self)
        spl.phyto_details(self)

        spl.get_depth_profiles(self)
        #

        spl.print_output(self)
        spl.get_carbon_export(self)
        #spl.ns_newplot(self)

################################################################################

    def __del__(self):

        print("splas.__del__ called")

################################################################################

###################################################################################################################################################################
#save output on file, to perform plot later...(using pickle module)
###################################################################################################################################################################

    def write_on_file(self,direct_outp):
            #
            pickle.dump(self.io, open(direct_outp + self.io.name + '_' + "json_file.p","wb"))
            pickle.dump(self.io.time_flag, open(direct_outp + self.io.name + '_' + "json_time_flag.p","wb"))
            pickle.dump(self.io.outdir, open(direct_outp + self.io.name + '_' + "json_output_directory.p","wb"))
            pickle.dump(self.io.name, open(direct_outp + self.io.name + '_' + "json_simulation_name.p","wb"))
            pickle.dump(self.io.phyto_max_size, open(direct_outp + self.io.name + '_' + "json_phyto_max_size.p","wb"))
            pickle.dump(self.io.zoo_max_size, open(direct_outp + self.io.name + '_' + "json_zoo_max_size.p","wb"))
            pickle.dump(self.io.phyto_min_size, open(direct_outp + self.io.name + '_' + "json_phyto_min_size.p","wb"))
            pickle.dump(self.io.zoo_min_size, open(direct_outp + self.io.name + '_' + "json_zoo_min_size.p","wb"))
            pickle.dump(self.io.phyto_num_class, open(direct_outp + self.io.name + '_' + "json_phyto_dim.p","wb"))
            pickle.dump(self.io.zoo_num_class, open(direct_outp + self.io.name + '_' + "json_zoo_dim.p","wb"))
            pickle.dump(self.io.simulation_length, open(direct_outp + self.io.name + '_' + "json_simulation_length.p","wb"))

            pickle.dump(self.nt, open(direct_outp + self.io.name + '_' + "nt.p","wb"))
            pickle.dump(self.ps, open(direct_outp + self.io.name + '_' + "printing_step.p","wb"))
            pickle.dump(self.tdim, open(direct_outp + self.io.name + '_' + "time_dimension.p","wb"))
            pickle.dump(self.days, open(direct_outp + self.io.name + '_' + "days.p","wb"))

            pickle.dump(self.io.dz, open(direct_outp + self.io.name + '_' + "delta_z.p","wb"))
            pickle.dump(self.io.m, open(direct_outp + self.io.name + '_' + "number_layers.p","wb"))
            pickle.dump(self.h, open(direct_outp + self.io.name + '_' + "total_depth.p","wb"))
            pickle.dump(self.z, open(direct_outp + self.io.name + '_' + "z_points.p","wb"))
            pickle.dump(self.diffusion_profile, open(direct_outp + self.io.name + '_' + "diff_vs_depth.p","wb"))

            pickle.dump(self.p_mort, open(direct_outp + self.io.name + '_' + "phyto_mortality.p","wb"))
            pickle.dump(self.psink, open(direct_outp + self.io.name + '_' + "phyto_sinking_rate.p","wb"))
            pickle.dump(self.z_dx, open(direct_outp + self.io.name + '_' + "zoo_selectivity.p","wb"))

            pickle.dump(self.p_ks, open(direct_outp + self.io.name + '_' + "phyto_half_saturation.p","wb"))
            pickle.dump(self.p_mu, open(direct_outp + self.io.name + '_' + "phyto_growth_rate.p","wb"))
            pickle.dump(self.p_r_star, open(direct_outp + self.io.name + '_' + "phyto_n_requirement.p","wb"))
            pickle.dump(self.p_vmax, open(direct_outp + self.io.name + '_' + "phyto_max_uptake_rate.p","wb"))
            pickle.dump(self.p_qmin, open(direct_outp + self.io.name + '_' + "phyto_min_quota.p","wb"))
            pickle.dump(self.p_cell_vol, open(direct_outp + self.io.name + '_' + "phyto_cell_vol.p","wb"))
            pickle.dump(self.p_naff, open(direct_outp + self.io.name + '_' + "phyto_n_affinity.p","wb"))
            pickle.dump(self.p_scaled_naff, open(direct_outp + self.io.name + '_' + "phyto_n_scal_affinity.p","wb"))

            pickle.dump(self.z_ir, open(direct_outp + self.io.name + '_' + "zoo_ingestion_rate.p","wb"))
            pickle.dump(self.z_xopt, open(direct_outp + self.io.name + '_' + "zoo_optimum_prey.p","wb"))

            pickle.dump(self.xp, open(direct_outp + self.io.name + '_' + "xp_phyto_ESD.p","wb"))
            pickle.dump(self.xz, open(direct_outp + self.io.name + '_' + "xz_zoo_ESD.p","wb"))

            #######
            pickle.dump(self.fphyto, open(direct_outp + self.io.name + '_' + "phyto_final_spectrum.p","wb"))
            pickle.dump(self.fzoo, open(direct_outp + self.io.name + '_' + "zoo_final_spectrum.p","wb"))
            pickle.dump(self.fpquota, open(direct_outp + self.io.name + '_' + "phyto_quota_final_spectrum.p","wb"))
            pickle.dump(self.fpncells, open(direct_outp + self.io.name + '_' + "phyto_ncells_final_spectrum.p","wb"))

            pickle.dump(self.nut_lim_final, open(direct_outp + self.io.name + '_' + "phyto_nutlim_final.p","wb"))
            pickle.dump(self.graz_lim_final, open(direct_outp + self.io.name + '_' + "phyto_grazlim_final.p","wb"))

            pickle.dump(self.grazing_zp, open(direct_outp + self.io.name + '_' + "grazing_zp.p","wb"))

            pickle.dump(self.phi, open(direct_outp + self.io.name + '_' + "phi_zp.p","wb"))

            pickle.dump(self.phyto, open(direct_outp + self.io.name + '_' + "phytoplankton.p","wb"))
            pickle.dump(self.zoo, open(direct_outp + self.io.name + '_' + "zooplankton.p","wb"))
            pickle.dump(self.pquota, open(direct_outp + self.io.name + '_' + "phytoplankton_quota.p","wb"))
            pickle.dump(self.pncells, open(direct_outp + self.io.name + '_' + "phytoplankton_ncells.p","wb"))

            pickle.dump(self.nutlim, open(direct_outp + self.io.name + '_' + "n_limitation.p","wb"))
            pickle.dump(self.grazlim, open(direct_outp + self.io.name + '_' + "graz_limitation.p","wb"))

            pickle.dump(self.aven, open(direct_outp + self.io.name + '_' + "n_limitation_taverage.p","wb"))
            pickle.dump(self.aveg, open(direct_outp + self.io.name + '_' + "graz_limitation_taverage.p","wb"))

            pickle.dump(self.Nitrogen, open(direct_outp + self.io.name + '_' + "nitrogen.p","wb"))
            pickle.dump(self.Phosphorus, open(direct_outp + self.io.name + '_' + "phosphorus.p","wb"))
            pickle.dump(self.Carbon, open(direct_outp + self.io.name + '_' + "carbon.p","wb"))

            pickle.dump(self.phyto_biodiv, open(direct_outp + self.io.name + '_' + "phyto_biodiv.p","wb"))
            pickle.dump(self.zoo_biodiv, open(direct_outp + self.io.name + '_' + "zoo_biodiv.p","wb"))

            pickle.dump(self.phyto_mean, open(direct_outp + self.io.name + '_' + "phyto_mean_spectrum.p","wb"))
            pickle.dump(self.zoo_mean, open(direct_outp + self.io.name + '_' + "zoo_mean_spectrum.p","wb"))

            pickle.dump(self.phyto_biomass, open(direct_outp + self.io.name + '_' + "phyto_biomass.p","wb"))
            pickle.dump(self.zoo_biomass, open(direct_outp + self.io.name + '_' + "zoo_biomass.p","wb"))

            pickle.dump(self.nit_uptake, open(direct_outp + self.io.name + '_' + "nitrogen_uptake.p","wb"))
            pickle.dump(self.p_growth_rate, open(direct_outp + self.io.name + '_' + "growth_rate.p","wb"))

            pickle.dump(self.z_predation, open(direct_outp + self.io.name + '_' + "z_total_predation.p","wb"))

            pickle.dump(self.small_phyto, open(direct_outp + self.io.name + '_' + "picophytoplankton.p","wb"))
            pickle.dump(self.phyto_total_mean_biomass, open(direct_outp + self.io.name + '_' + "phytoplankton_total_mean_biomass.p","wb"))
            pickle.dump(self.zoo_total_mean_biomass, open(direct_outp + self.io.name + '_' + "zooplankton_total_mean_biomass.p","wb"))
            pickle.dump(self.mean_shannon_p, open(direct_outp + self.io.name + '_' + "mean_shannon_p.p","wb"))
            pickle.dump(self.mean_shannon_z, open(direct_outp + self.io.name + '_' + "mean_shannon_z.p","wb"))
            pickle.dump(self.mean_N, open(direct_outp + self.io.name + '_' + "mean_N.p","wb"))
            pickle.dump(self.final_N, open(direct_outp + self.io.name + '_' + "final_N.p","wb"))
            pickle.dump(self.irrprofile, open(direct_outp + self.io.name + '_' + "irradiance_profile.p","wb"))
            #DETRITUS
            pickle.dump(self.detritus, open(direct_outp + self.io.name + '_' + "detritus.p","wb"))
            pickle.dump(self.final_detritus, open(direct_outp + self.io.name + '_' + "final_detritus.p","wb"))
            pickle.dump(self.mean_detritus, open(direct_outp + self.io.name + '_' + "mean_detritus.p","wb"))
            #####

            bio_dict = {"picophytoplankton_biomass":self.avez_small_phyto,
                        "phytoplankton_total_mean_biomass":self.avez_p_biomass,
                        "zooplankton_total_mean_biomass":self.avez_z_biomass,
                        "phytoplankton_mean_shannon_evenness":self.avez_shannon_p,
                        "zooplankton_mean_shannon_evenness":self.avez_shannon_z,
                        "Nitrogen_mean_concentration":self.avez_N,
                        "mean_z_p_ratio":self.avez_z_p}
            pickle.dump(bio_dict, open(direct_outp + self.io.name + '_' + "D_biomass_biodiv.p","wb"))

            pickle.dump(self.z_p_ratio, open(direct_outp + self.io.name + '_' + "z_p_ratio.p","wb"))
            pickle.dump(self.ave_z_p_ratio, open(direct_outp + self.io.name + '_' + "ave_z_p_ratio.p","wb"))

            #####CARBON
            carbon_export_dict =    {"PHYTOPLANKTON AVERAGE CARBON EXPORT (mgC m-2 day-1)":self.phyto_ave_total_c_exp,
                                     "DETRITUS AVERAGE CARBON EXPORT (mgC m-2 day-1)":self.zoo_ave_total_c_exp,
                                     "TOTAL AVERAGE CARBON EXPORT (mgC m-2 day-1)":self.total_ave_carbon_export,
                                     "PHYTOPLANKTON INTEGRATED CARBON EXPORT (mgC m-2)":self.phyto_integrated_total_c_exp,
                                     "DETRITUS INTEGRATED CARBON EXPORT (mgC m-2)":self.zoo_integrated_total_c_exp,
                                     "TOTAL INTEGRATED CARBON EXPORT (mgC m-2)":self.total_integrated_carbon_export}
            pickle.dump(carbon_export_dict, open(direct_outp + self.io.name + '_' + "D_carbon_export.p","wb"))

            pickle.dump(self.p_carbon, open(direct_outp + self.io.name + '_' + "phytoplankton_carbon.p","wb"))
            pickle.dump(self.p_carbon_export_matrix, open(direct_outp + self.io.name + '_' + "phyto_carbon_export_matrix.p","wb"))
            pickle.dump(self.p_total_carbon_export, open(direct_outp + self.io.name + '_' + "phyto_carbon_export.p","wb"))
            pickle.dump(self.z_carbon_export, open(direct_outp + self.io.name + '_' + "detritus_carbon_export.p","wb"))

            pickle.dump(self.total_carbon_export, open(direct_outp + self.io.name + '_' + "total_carbon_export.p","wb"))
            pickle.dump(self.phyto_ave_c_exp, open(direct_outp + self.io.name + '_' + "phyto_carbon_export_spectrum_taverage.p","wb"))
            pickle.dump(self.phyto_integrated_c_exp, open(direct_outp + self.io.name + '_' + "phyto_carbon_export_spectrum_tintegrated.p","wb"))

            #write also biomass
            f1 = open("./OUTPUT/" + self.io.name + '_' + 'BIOMASS.txt','w')
            f1.write("{}".format(bio_dict))
            f1.close()
            #write carbon export on file when storing
            f2 = open("./OUTPUT/" + self.io.name + '_' + 'C_export.txt','w')
            f2.write("{}".format(carbon_export_dict))
            f2.close()
            #
################################################################################

###################################################################################################################################################################
#load output from file, to perform plot decoupled from simulation...(using pickle module)
###################################################################################################################################################################

    def load_file(self,direct,name):

        self.io = pickle.load(open(direct + name + '_' + "json_file.p","rb"))
        self.io.time_flag = pickle.load(open(direct + name + '_' + "json_time_flag.p","rb"))
        self.io.outdir = pickle.load(open(direct + name + '_' + "json_output_directory.p","rb"))
        self.io.name = pickle.load(open(direct + name + '_' + "json_simulation_name.p","rb"))
        self.io.phyto_max_size = pickle.load(open(direct + name + '_' + "json_phyto_max_size.p","rb"))
        self.io.zoo_max_size = pickle.load(open(direct + name + '_' + "json_zoo_max_size.p","rb"))
        self.io.phyto_min_size = pickle.load(open(direct + name + '_' + "json_phyto_min_size.p","rb"))
        self.io.zoo_min_size = pickle.load(open(direct + name + '_' + "json_zoo_min_size.p","rb"))
        self.io.phyto_num_class = pickle.load(open(direct + name + '_' + "json_phyto_dim.p","rb"))
        self.io.zoo_num_class = pickle.load(open(direct + name + '_' + "json_zoo_dim.p","rb"))
        self.io.simulation_length = pickle.load(open(direct + name + '_' + "json_simulation_length.p","rb"))

        self.nt = pickle.load(open(direct + name + '_' + "nt.p","rb"))
        self.ps = pickle.load(open(direct + name + '_' + "printing_step.p","rb"))
        self.tdim = pickle.load(open(direct + name + '_' + "time_dimension.p","rb"))
        self.days = pickle.load(open(direct + name + '_' + "days.p","rb"))

        self.io.dz = pickle.load(open(direct + name + '_' + "delta_z.p","rb"))
        self.io.m = pickle.load(open(direct + name + '_' + "number_layers.p","rb"))
        self.h = pickle.load(open(direct + name + '_' + "total_depth.p","rb"))
        self.z = pickle.load(open(direct + name + '_' + "z_points.p","rb"))
        self.diffusion_profile = pickle.load(open(direct + name + '_' + "diff_vs_depth.p","rb"))

        self.p_mort = pickle.load(open(direct + name + '_' + "phyto_mortality.p","rb"))
        self.psink = pickle.load(open(direct + name + '_' + "phyto_sinking_rate.p","rb"))
        self.z_dx = pickle.load(open(direct + name + '_' + "zoo_selectivity.p","rb"))

        self.p_ks = pickle.load(open(direct + name + '_' + "phyto_half_saturation.p","rb"))
        self.p_mu = pickle.load(open(direct + name + '_' + "phyto_growth_rate.p","rb"))
        self.p_r_star = pickle.load(open(direct + name + '_' + "phyto_n_requirement.p","rb"))
        self.p_vmax = pickle.load(open(direct + name + '_' + "phyto_max_uptake_rate.p","rb"))
        self.p_qmin = pickle.load(open(direct + name + '_' + "phyto_min_quota.p","rb"))
        self.p_cell_vol = pickle.load(open(direct + name + '_' + "phyto_cell_vol.p","rb"))
        self.p_naff = pickle.load(open(direct + name + '_' + "phyto_n_affinity.p","rb"))
        self.p_scaled_naff = pickle.load(open(direct + name + '_' + "phyto_n_scal_affinity.p","rb"))

        self.z_ir = pickle.load(open(direct + name + '_' + "zoo_ingestion_rate.p","rb"))
        self.z_xopt = pickle.load(open(direct + name + '_' + "zoo_optimum_prey.p","rb"))

        self.xp = pickle.load(open(direct + name + '_' + "xp_phyto_ESD.p","rb"))
        self.xz = pickle.load(open(direct + name + '_' + "xz_zoo_ESD.p","rb"))

        self.fphyto = pickle.load(open(direct + name + '_' + "phyto_final_spectrum.p","rb"))
        self.fzoo = pickle.load(open(direct + name + '_' + "zoo_final_spectrum.p","rb"))
        self.fpquota = pickle.load(open(direct + name + '_' + "phyto_quota_final_spectrum.p","rb"))
        self.fpncells = pickle.load(open(direct + name + '_' + "phyto_ncells_final_spectrum.p","rb"))

        self.nut_lim_final = pickle.load(open(direct + name + '_' + "phyto_nutlim_final.p","rb"))
        self.graz_lim_final = pickle.load(open(direct + name + '_' + "phyto_grazlim_final.p","rb"))

        self.grazing_zp = pickle.load(open(direct + name + '_' + "grazing_zp.p","rb"))

        self.phi = pickle.load(open(direct + name + '_' + "phi_zp.p","rb"))

        self.phyto = pickle.load(open(direct + name + '_' + "phytoplankton.p","rb"))
        self.zoo = pickle.load(open(direct + name + '_' + "zooplankton.p","rb"))
        self.pquota = pickle.load(open(direct + name + '_' + "phytoplankton_quota.p","rb"))

        self.pncells = pickle.load(open(direct + name + '_' + "phytoplankton_ncells.p","rb"))

        self.nutlim = pickle.load(open(direct + name + '_' + "n_limitation.p","rb"))
        self.grazlim = pickle.load(open(direct + name + '_' + "graz_limitation.p","rb"))

        self.aven = pickle.load(open(direct + name + '_' + "n_limitation_taverage.p","rb"))
        self.aveg = pickle.load(open(direct + name + '_' + "graz_limitation_taverage.p","rb"))

        self.Nitrogen = pickle.load(open(direct + name + '_' + "nitrogen.p","rb"))
        self.Phosphorus = pickle.load(open(direct + name + '_' + "phosphorus.p","rb"))
        self.Carbon = pickle.load(open(direct + name + '_' + "carbon.p","rb"))

        self.phyto_biodiv = pickle.load(open(direct + name + '_' + "phyto_biodiv.p","rb"))
        self.zoo_biodiv = pickle.load(open(direct + name + '_' + "zoo_biodiv.p","rb"))

        self.phyto_mean = pickle.load(open(direct + name + '_' + "phyto_mean_spectrum.p","rb"))
        self.zoo_mean = pickle.load(open(direct + name + '_' + "zoo_mean_spectrum.p","rb"))

        self.phyto_biomass = pickle.load(open(direct + name + '_' + "phyto_biomass.p","rb"))
        self.zoo_biomass = pickle.load(open(direct + name + '_' + "zoo_biomass.p","rb"))

        self.z_p_ratio = pickle.load(open(direct + name + '_' + "z_p_ratio.p","rb"))
        self.ave_z_p_ratio = pickle.load(open(direct + name + '_' + "ave_z_p_ratio.p","rb"))

        self.nit_uptake = pickle.load(open(direct + name + '_' + "nitrogen_uptake.p","rb"))
        self.p_growth_rate = pickle.load(open(direct + name + '_' + "growth_rate.p","rb"))

        self.z_predation = pickle.load(open(direct + name + '_' + "z_total_predation.p","rb"))

        self.small_phyto = pickle.load(open(direct + name + '_' + "picophytoplankton.p","rb"))
        self.phyto_total_mean_biomass = pickle.load(open(direct + name + '_' + "phytoplankton_total_mean_biomass.p","rb"))
        self.zoo_total_mean_biomass = pickle.load(open(direct + name + '_' + "zooplankton_total_mean_biomass.p","rb"))
        self.mean_shannon_p = pickle.load(open(direct + name + '_' + "mean_shannon_p.p","rb"))
        self.mean_shannon_z = pickle.load(open(direct + name + '_' + "mean_shannon_z.p","rb"))
        self.mean_N = pickle.load(open(direct + name + '_' + "mean_N.p","rb"))
        self.final_N = pickle.load(open(direct + name + '_' + "final_N.p","rb"))

        self.p_carbon = pickle.load(open(direct + name + '_' + "phytoplankton_carbon.p","rb"))
        self.p_carbon_export_matrix = pickle.load(open(direct + name + '_' + "phyto_carbon_export_matrix.p","rb"))
        self.p_total_carbon_export = pickle.load(open(direct + name + '_' + "phyto_carbon_export.p","rb"))
        self.z_carbon_export = pickle.load(open(direct + name + '_' + "detritus_carbon_export.p","rb"))

        self.total_carbon_export = pickle.load(open(direct + name + '_' + "total_carbon_export.p","rb"))
        self.phyto_ave_c_exp = pickle.load(open(direct + name + '_' + "phyto_carbon_export_spectrum_taverage.p","rb"))
        self.phyto_integrated_c_exp = pickle.load(open(direct + name + '_' + "phyto_carbon_export_spectrum_tintegrated.p","rb"))
        self.irrprofile = pickle.load(open(direct + name + '_' + "irradiance_profile.p","rb"))
        #DETRITUS
        self.detritus = pickle.load(open(direct + name + '_' + "detritus.p","rb"))
        self.final_detritus = pickle.load(open(direct + name + '_' + "final_detritus.p","rb"))
        self.mean_detritus = pickle.load(open(direct + name + '_' + "mean_detritus.p","rb"))
        #####

        self.bio_dict = pickle.load(open(direct + name + '_' + "D_biomass_biodiv.p","rb"))
        self.carbon_export_dict = pickle.load(open(direct + name + '_' + "D_carbon_export.p","rb"))

        f1 = open(name + '_' + 'biomass.txt','w')
        f2 = open(name + '_' + 'C_export.txt','w')

        f1.write("{}".format(self.bio_dict))
        f2.write("{}".format(self.carbon_export_dict))

        f1.close()
        f2.close()

################################################################################

if __name__ == "__main__":

    script, namesim = argv
    splas.load_file(splas,dir_pickle,namesim)

    spl.get_parameter_plot(splas)
    #
    spl.get_final_plot(splas)
    spl.get_time_plot(splas)
    spl.phyto_details(splas)
    spl.get_depth_profiles(splas)
    #
    spl.get_carbon_export(splas)

################################################################################
