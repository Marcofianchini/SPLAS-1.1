import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
#import pylab as p
import mpl_toolkits.mplot3d as p3
from matplotlib.font_manager import FontProperties

####################################################################################################################################

def get_parameter_plot(inclass):

    print("Plot of model parameters (ks,mu_0,r*,ir,x_opt)")

    X = inclass.xp
    Y = inclass.xz

    #p half saturation
    plt.plot(X,inclass.p_ks,'lime')
    plt.grid()
    plt.axis()
    plt.xlim([0,inclass.io.phyto_max_size])
    plt.xlabel('Size (micro-meters)')
    #plt.title('Parameters')
    plt.ylim([0,(np.amax(inclass.p_ks)+1)]) #
    plt.ylabel('Half-Saturation Constant (uMN)')
    #fontpar = FontProperties()
    #fontpar.set_size('small')
    #plt.legend(['half-saturation (uMN)','R* (uMN)'], loc = 2, prop=fontpar, bbox_to_anchor = (-0.01, 1.12), ncol=2)
    plt.savefig(inclass.io.outdir + '/SPLAS_PARAMETERS/' + inclass.io.name + '_' + 'Phyto_HS.png', dpi = 250)
    #plt.show()
    plt.clf()

    #N equilibrium concentration (dP/dt=0-->Q* and then dQ/dt=0-->N*(Q*))
    plt.plot(X,inclass.p_r_star,'m')
    plt.grid()
    plt.axis()
    plt.xlim([0,inclass.io.phyto_max_size])
    plt.xlabel('Size (micro-meters)')
    plt.title('N equilibrium concentration (for dQ/dt = 0)')
    plt.ylim([0,(np.amax(inclass.p_r_star)+1)]) #
    plt.ylabel('N* (uMN)')
    #fontpar2 = FontProperties()
    #fontpar2.set_size('small')
    #plt.legend(['half-saturation (uMN)','R* (uMN)'], loc = 2, prop=fontpar2, bbox_to_anchor = (-0.01, 1.12), ncol=2)
    plt.savefig(inclass.io.outdir + '/SPLAS_PARAMETERS/' + inclass.io.name + '_' + 'Phyto_N_eq.png', dpi = 250)
    #plt.show()
    plt.clf()

    #p max growth rate
    plt.plot(X,inclass.p_mu,'orange')
    plt.grid()
    plt.axis()
    plt.axis([0,inclass.io.phyto_max_size+1,0,np.amax(inclass.p_mu)+0.2])
    plt.xlabel('Size (micro-meters)')
    #plt.title('Parameters')
    plt.ylabel('Maximum Growth Rate (day-1)')
    #fontpar = FontProperties()
    #fontpar.set_size('small')
    #plt.legend([], loc = 2, prop=fontpar, bbox_to_anchor = (-0.01, 1.12), ncol=3)
    plt.savefig(inclass.io.outdir + '/SPLAS_PARAMETERS/' + inclass.io.name + '_' + 'Phyto_Growth_Rate.png', dpi = 250)
    #plt.show()
    plt.clf()

    #p max growth rate (detail for picophytoplankton)
    plt.plot(X,inclass.p_mu,'orange')
    plt.grid()
    plt.axis()
    plt.axis([0,2,0,np.amax(inclass.p_mu)+0.2])
    plt.xlabel('Size (micro-meters)')
    plt.title('picophytoplankton')
    plt.ylabel('Maximum Growth Rate (day-1)')
    plt.savefig(inclass.io.outdir + '/SPLAS_PARAMETERS/' + inclass.io.name + '_' + 'picophytoplankton_Growth_Rate.png', dpi = 250)
    plt.clf()

    #z optimum prey size
    plt.plot(Y,inclass.z_xopt,'b')
    plt.errorbar(Y,inclass.z_xopt,yerr=10**(inclass.z_dx))
    #plt.errorbar(Y,inclass.z_xopt,yerr=inclass.z_dx)
    plt.grid()
    plt.axis()
    plt.axis([0,inclass.io.zoo_max_size+2,0,inclass.io.phyto_max_size+3])
    plt.xlabel('Zooplankton Size (micro-meters)')
    #plt.title('Parameters')
    plt.ylabel('Prey Size (micro-meters)')
    #plt.legend(['ingestion rate (day-1)','optimum prey size (micro-meters)'], loc = 2, prop=fontpar, bbox_to_anchor = (-0.01, 1.12), ncol=2)
    plt.savefig(inclass.io.outdir + '/SPLAS_PARAMETERS/' + inclass.io.name + '_' + 'Zoo_OPS.png', dpi = 250)
    #plt.show()
    plt.clf()

    #z ingestion rate
    plt.plot(Y,inclass.z_ir,'g')
    plt.grid()
    plt.axis()
    plt.axis([0,inclass.io.zoo_max_size+2,0,np.amax(inclass.z_ir)+0.2])
    plt.xlabel('Zooplankton Size (micro-meters)')
    #plt.title('Parameters')
    plt.ylabel('Maximum Ingestion Rate (day-1)')
    #plt.legend(['ingestion rate (day-1)','optimum prey size (micro-meters)'], loc = 2, prop=fontpar, bbox_to_anchor = (-0.01, 1.12), ncol=2)
    plt.savefig(inclass.io.outdir + '/SPLAS_PARAMETERS/' + inclass.io.name + '_' + 'Zoo_Ingestion_Rate.png', dpi = 250)
    #plt.show()
    plt.clf()

    #p sinking and natural mortality terms
    one = np.ones(inclass.io.phyto_num_class)
    sink = inclass.psink*(X**1.17)
    plt.plot(X,sink,'m',X,inclass.p_mort*one,'r--')
    plt.grid()
    plt.axis()
    plt.axis([inclass.io.phyto_min_size,inclass.io.phyto_max_size,0,np.amax(sink)+0.07])
    plt.xlabel('Phytoplankton Size (micro-meters)')
    #plt.title('Parameters')
    plt.ylabel('Phytoplankton Mortality')
    fontpar3 = FontProperties()
    fontpar3.set_size('small')
    plt.legend(['Sinking (day-1)','Size independent mortality (day-1)'], loc = 2, prop=fontpar3, bbox_to_anchor = (-0.01, 1.12), ncol=2)
    plt.savefig(inclass.io.outdir + '/SPLAS_PARAMETERS/' + inclass.io.name + '_' + 'Phyto_Mortality.png', dpi = 250)
    #plt.show()
    plt.clf()

    #p uptake rate and min internal cellular quota
    plt.plot(X,inclass.p_vmax,'lime',X,inclass.p_qmin,'b')
    plt.grid()
    plt.axis()
    plt.axis([inclass.io.phyto_min_size,inclass.io.phyto_max_size,0,np.amax(inclass.p_vmax)])
    plt.xlabel('Phytoplankton Size (micro-meters)')
    #plt.title('Parameters')
    plt.ylabel('Phytoplankton Parameters')
    fontpar4 = FontProperties()
    fontpar4.set_size('small')
    plt.legend(['uptake rate (mmolN day-1 cell-1)','min internal quota (mmolN cell-1)'], loc = 2, prop=fontpar4, bbox_to_anchor = (-0.01, 1.12), ncol=2)
    plt.savefig(inclass.io.outdir + '/SPLAS_PARAMETERS/' + inclass.io.name + '_' + 'Phyto_NEW_parameters.png', dpi = 250)
    #plt.show()
    plt.clf()

####################################################################################################

    #uptake affinity and scaled uptake affinity
    plt.plot(np.log10(X),np.log10(inclass.p_naff),'c',np.log10(X),np.log10(inclass.p_scaled_naff),'b--')
    plt.grid()
    plt.axis()
    plt.axis([np.log10(inclass.io.phyto_min_size),np.log10(inclass.io.phyto_max_size),-12.0,12.0])
    plt.xlabel('Log10[Phytoplankton Size (micro-meters)]')
    #plt.title('Parameters')
    plt.ylabel('Log10[Phytoplankton Uptake Affinity]')
    fontpar5 = FontProperties()
    fontpar5.set_size('small')
    plt.legend(['Log10 [N affinity] [m3 day-1 cell-1]','Log10 [N scaled affinity] [m3 mmolN-1 cell-1]'], loc = 2, prop=fontpar5, bbox_to_anchor = (-0.01, 1.12), ncol=2)
    plt.savefig(inclass.io.outdir + '/SPLAS_PARAMETERS/' + inclass.io.name + '_' + 'Phyto_AFFINITY.png', dpi = 250)
    #plt.show()
    plt.clf()

    #cellular carbon content
    plt.plot(X,inclass.p_carbon,'r')
    plt.grid()
    plt.axis()
    plt.axis([inclass.io.phyto_min_size,inclass.io.phyto_max_size,0,np.amax(inclass.p_carbon)+0.05])
    plt.xlabel('Phytoplankton Size (micro-meters)')
    #plt.title('Parameters')
    plt.ylabel('Phytoplankton C content (pgC cell-1)')
    #fontpar6 = FontProperties()
    #fontpar6.set_size('small')
    #plt.legend(['uptake rate (umolN day-1 cell-1)','min internal quota (umolN cell-1)'], loc = 2, prop=fontpar6, bbox_to_anchor = (-0.01, 1.12), ncol=2)
    plt.savefig(inclass.io.outdir + '/SPLAS_PARAMETERS/' + inclass.io.name + '_' + 'Phyto_CARBON.png', dpi = 250)
    #plt.show()
    plt.clf()

##################be careful! All together
    plt.plot(X,inclass.p_ks,'lime',X,inclass.p_mu,'orange',X,inclass.p_r_star,'m--')
    plt.grid()
    plt.axis()
    plt.axis([inclass.io.phyto_min_size,inclass.io.phyto_max_size,0,np.amax(inclass.p_mu)+0.3])
    plt.xlabel('Size (micro-meters)')
    #plt.title('Parameters')
    plt.ylabel('Phytoplankton Parameters and Limitations')
    fontpar2 = FontProperties()
    fontpar2.set_size('small')
    plt.legend(['half-saturation (uMN)','max growth rate (day-1)','R* (uMN)'], loc = 2, prop=fontpar2, bbox_to_anchor = (-0.01, 1.12), ncol=3)
    plt.savefig(inclass.io.outdir + '/SPLAS_PARAMETERS/' + inclass.io.name + '_' + 'Phyto_TRADE_OFF.png', dpi = 250)
    #plt.show()
    plt.clf()

####################################################################################################################################

    plt.close('all')

def get_final_plot(inclass):
#now will depend on layer! (DEPTH)

    for k in range(inclass.io.m):
        I = np.arange(1,inclass.io.phyto_num_class+1,1)
        pp = inclass.fphyto[:,k]*1000/np.amax(inclass.fphyto[:,k])*1000
        zz = inclass.fzoo[:,k]*1000/np.amax(inclass.fzoo[:,k])*1000
        plt.plot(I,pp,'gs',I,zz,'r^')
        plt.plot(I,pp,'g--',I,zz,'r--')
        plt.grid()
        plt.axis()
        plt.axis([1,inclass.io.phyto_num_class,0.0,1.01])
        plt.xlabel('Size-class Index')
        plt.ylabel('Normalized Biomass')
        plt.title('Layer ' + str(k) + 'Rel Size Spectrum')
        fontparpz = FontProperties()
        fontparpz.set_size('small')
        plt.legend(['P','Z'], loc = 2, prop=fontparpz, bbox_to_anchor = (-0.01, 1.12), ncol=2)
        plt.savefig(inclass.io.outdir + '/P_VS_Z_REL/' + inclass.io.name + '_' + 'P_vs_Z_REL_layer'+ str(k) +'_.png', dpi = 250)
        #plt.show()
        plt.clf()

        pp2 = inclass.fphyto[:,k]
        zz2 = inclass.fzoo[:,k]
        plt.plot(I,pp2,'gs',I,zz2,'r^')
        plt.plot(I,pp2,'g--',I,zz2,'r--')
        plt.grid()
        plt.axis()
        plt.xlim([1,inclass.io.phyto_num_class])
        #plt.axis([1,inclass.io.phyto_num_class,0.0,np.amax(pp2)+0.01])
        plt.xlabel('Size-class Index')
        plt.ylabel('Biomass (uMN)')
        plt.title('Layer ' + str(k) + 'Size Spectrum')
        fontparpz2 = FontProperties()
        fontparpz2.set_size('small')
        plt.legend(['P','Z'], loc = 2, prop=fontparpz2, bbox_to_anchor = (-0.01, 1.12), ncol=2)
        plt.savefig(inclass.io.outdir + '/P_VS_Z/' + inclass.io.name + '_' + 'P_vs_Z_layer'+ str(k) +'_.png', dpi = 250)
        #plt.show()
        plt.clf()

        print("Plot of final time quantities, after ",inclass.io.simulation_length," days")

        X = inclass.xp
        Y = inclass.xz

        #phytoplankton final time plot
        plt.plot(X,inclass.fphyto[:,k],'g')
        #plt.plot(X,inclass.fphyto[:,k],'k--')
        plt.grid()
        plt.axis()
        plt.xlim([0,inclass.io.phyto_max_size+1])
        plt.xlabel('Size (micro-meters)')
        plt.ylabel('Biomass (uMN)')
        plt.title('Layer ' + str(k) + 'P Size Spectrum')
        #plt.legend(['P'], loc = 2, bbox_to_anchor = (-0.02, 1.02))
        plt.savefig(inclass.io.outdir + '/SPLAS_SIZE_SPECTRA/' + inclass.io.name + '_' + 'P_layer'+ str(k) +'_.png', dpi = 250)
        #plt.show()
        plt.clf()

        #zooplankton final time plot
        plt.plot(Y,inclass.fzoo[:,k],'r')
        #plt.plot(Y,inclass.fzoo[:,k],'k--')
        plt.grid()
        plt.axis()
        plt.xlim([0,inclass.io.zoo_max_size+2])
        plt.xlabel('Size (micro-meters)')
        plt.ylabel('Biomass (uMN)')
        plt.title('Layer ' + str(k) + 'Z Size Spectrum')
        #plt.legend(['Z'], loc = 2, bbox_to_anchor = (-0.02, 1.02))
        plt.savefig(inclass.io.outdir + '/SPLAS_SIZE_SPECTRA/' + inclass.io.name + '_' + 'Z_layer'+ str(k) +'_.png', dpi = 250)
        #plt.show()
        plt.clf()

    ###########################################################

        #phytoplankton grazing and nutrient limitations at final time
        plt.plot(X,inclass.nut_lim_final[:,k],'c',X,inclass.graz_lim_final[:,k],'y--')
        plt.grid()
        plt.axis()
        plt.xlim([inclass.io.phyto_min_size,inclass.io.phyto_max_size])
        plt.xlabel('Size (micro-meters)')
        plt.ylabel('Limitation')
        plt.title('Layer ' + str(k) + 'P Limitation (Final Time)')
        fontngl = FontProperties()
        fontngl.set_size('small')
        plt.legend(['Nutrient','Grazing'], loc = 2, prop = fontngl, bbox_to_anchor = (-0.01, 1.12), ncol = 2)
        plt.savefig(inclass.io.outdir + '/SPLAS_P_LIMIT/' + inclass.io.name + '_' + 'P_finalNGL_layer'+ str(k) +'_.png', dpi = 250)
        #plt.show()
        plt.clf()

    #########################################################################################################################################

        X, Y = np.meshgrid(X, Y)

    #########################################################################################################################################

    ######################################################################### ZOOPLANKTON

        #2D colorplot of Z Grazing on P
        Zboh = inclass.grazing_zp[:,:,k]/(np.amax(inclass.grazing_zp)*0.01)
        Z = Zboh.transpose()
        cmap= plt.get_cmap('PiYG')
        #norm = plt.BoundaryNorm(MaxNlocator(nbins=15).tick_values(np.amin(Z),np.amax(Z)), ncolors=cmap.N, clip=True)
        figgraz = plt.figure()
        plt.pcolormesh(X,Y,Z,cmap=cmap)
        plt.xlabel('Size of Phytoplankton (micro-meters)')
        plt.ylabel('Size of Zooplankton (micro-meters)')
        plt.axis([inclass.io.phyto_min_size,inclass.io.phyto_max_size,inclass.io.zoo_min_size,inclass.io.zoo_max_size])
        plt.title('Layer ' + str(k) + 'Z Grazing on P (Final Time)')
        gbar = plt.colorbar()
        gbar.set_label('Relative Grazing Flux')
        figgraz.savefig(inclass.io.outdir + '/SPLAS_GRAZING/' + inclass.io.name + '_' + 'ZgrazP_2D_layer'+ str(k) +'_.png',dpi=250)
        #plt.show()
        plt.clf()

################################################################################################################
    #2D colorplot of Z Preference for P
    ZF = inclass.phi
    figpref = plt.figure()
    plt.pcolormesh(Y,X,ZF.transpose(),vmin=0.0,vmax=1.0)
    plt.ylabel('Size of Phytoplankton (micro-meters)')
    plt.xlabel('Size of Zooplankton (micro-meters)')
    plt.axis([inclass.io.zoo_min_size,inclass.io.zoo_max_size,inclass.io.phyto_min_size,inclass.io.phyto_max_size])
    plt.title('Zooplankton Preference for Phytoplankton')
    prbar = plt.colorbar()
    prbar.set_label('Preference ')
    figpref.savefig(inclass.io.outdir + '/SPLAS_PARAMETERS/' + inclass.io.name + '_' + 'ZpreferenceforP_2D.png',dpi=250)
    #plt.show()
    plt.clf()

################################################################################

    plt.close('all')

####################################################################################################################################

def get_time_plot(inclass):
#now will depend on layer! (DEPTH)

    for k in range(inclass.io.m):
            print("Color-plot of Phytoplankton Biomass Size Spectrum in Time")

            #2D colorplot of Phytoplankton
            XC = inclass.days
            YC = inclass.xp
            XC, YC = np.meshgrid(XC, YC)
            ZC = inclass.phyto[:,:,k]
            figp = plt.figure()
            plt.pcolormesh(XC,YC,ZC, vmin=0,vmax=np.amax(ZC) )
            plt.xlabel('Time (days)')
            plt.ylabel('Size (micro-meters)')
            plt.axis([0,inclass.io.simulation_length,inclass.io.phyto_min_size,inclass.io.phyto_max_size])
            plt.title('Layer ' + str(k) + 'P Size Spectrum')
            pbar = plt.colorbar()
            pbar.set_label('Biomass (uMN)')
            figp.savefig(inclass.io.outdir + '/SPLAS_SIZE_SPECTRA/' + inclass.io.name + '_' + 'P_2D_layer'+ str(k) +'_.png',dpi=250)
            #plt.show()
            plt.clf()

            print("Color-plot of Zooplankton Biomass Size Spectrum in Time")

            #2D colorplot of Zooplankton
            XD = inclass.days
            YD = inclass.xz
            XD, YD = np.meshgrid(XD, YD)
            ZD = inclass.zoo[:,:,k]
            figz = plt.figure()
            plt.pcolormesh(XD,YD,ZD, vmin=0, vmax= np.amax(ZD))
            plt.xlabel('Time (days)')
            plt.ylabel('Size (micro-meters)')
            plt.axis([0,inclass.io.simulation_length,inclass.io.zoo_min_size,inclass.io.zoo_max_size])
            plt.title('Layer ' + str(k) + 'Z Size Spectrum')
            zbar = plt.colorbar()
            zbar.set_label('Biomass (uMN)')
            figz.savefig(inclass.io.outdir + '/SPLAS_SIZE_SPECTRA/' + inclass.io.name + '_' + 'Z_2D_layer'+ str(k) +'_.png',dpi=250)
            #plt.show()
            plt.clf()

            print("Biomass Size Spectrum in Time, NORMALIZED Color-plot with respect to maximum at the time")

            #2D RELATIVE TO MAXIMUM AT THAT TIME colorplot of Phytoplankton
            Xrp = inclass.days
            Yrp = inclass.xp
            Xrp, Yrp = np.meshgrid(Xrp, Yrp)
            Zrp = inclass.phyto[:,:,k]*1000/np.amax(inclass.phyto[:,:,k],0)*1000
            figrp = plt.figure()
            plt.pcolormesh(Xrp,Yrp,Zrp,vmin=0.0,vmax=1.0)
            plt.xlabel('Time (days)')
            plt.ylabel('Size (micro-meters)')
            plt.axis([0,inclass.io.simulation_length,inclass.io.phyto_min_size,inclass.io.phyto_max_size])
            plt.title('Layer ' + str(k) + 'P Relative Size Spectrum')
            prbar = plt.colorbar()
            prbar.set_label('Relative Fraction of Biomass')
            figrp.savefig(inclass.io.outdir + '/SPLAS_SIZE_SPECTRA/' + inclass.io.name + '_' + 'P_rel_2D_layer'+ str(k) +'_.png',dpi=250)
            #plt.show()
            plt.clf()

            #2D RELATIVE TO MAXIMUM AT THAT TIME colorplot of Zooplankton
            Xrz = inclass.days
            Yrz = inclass.xz
            Xrz, Yrz = np.meshgrid(Xrz, Yrz)
            Zrz = inclass.zoo[:,:,k]*1000/np.amax(inclass.zoo[:,:,k],0)*1000
            figrz = plt.figure()
            plt.pcolormesh(Xrz,Yrz,Zrz,vmin=0.0,vmax=1.0)
            plt.xlabel('Time (days)')
            plt.ylabel('Size (micro-meters)')
            plt.axis([0,inclass.io.simulation_length,inclass.io.zoo_min_size,inclass.io.zoo_max_size])
            plt.title('Layer ' + str(k) + 'Z Relative Size Spectrum')
            zrbar = plt.colorbar()
            zrbar.set_label('Relative Fraction of Biomass')
            figrz.savefig(inclass.io.outdir + '/SPLAS_SIZE_SPECTRA/' + inclass.io.name + '_' + 'Z_rel_2D_layer'+ str(k) +'_.png',dpi=250)
            #plt.show()
            plt.clf()

        ################################################################################

            print("Plot of Nutrients in time")

            #nitrogen vs time
            N = inclass.Nitrogen[:,k]
            plt.plot(inclass.days,N,color='lightcoral')
            plt.grid()
            plt.axis()
            plt.axis([0,inclass.io.simulation_length,0,np.amax(N)+0.3])
            plt.xlabel('Time (days)')
            plt.ylabel('N Concentration (uMN)')
            plt.title('Layer ' + str(k) + '[N] in time')
            #plt.legend(['N'], loc = 2, bbox_to_anchor = (-0.02, 1.02))
            plt.savefig(inclass.io.outdir + '/SPLAS_BIOMASS/' + inclass.io.name + '_' + 'N_t_layer'+ str(k) +'_.png', dpi = 250)
            #plt.show()
            plt.clf()

        #########################################################################################fil
            print("Plot of Irradiance in time")

            #irradiance vs time
            I = inclass.irrprofile[:,k]
            plt.plot(inclass.days,I,color='lightcoral')
            plt.grid()
            plt.axis()
            plt.axis([0,inclass.io.simulation_length,0,np.amax(I)+0.3])
            plt.xlabel('Time (days)')
            plt.ylabel('Irradiance (W/m^2)')
            plt.title('Layer' + str(k) + 'Irradiance in time')
            plt.savefig(inclass.io.outdir + '/SPLAS_IRRADIANCE/' + inclass.io.name + '_' + 'I_t_layer' + str(k) + '_.png', dpi = 250)

        #########################################################################################

        ################### BE CAREFUL ############################## BE CAREFUL ################

            print("Color-plot of Phytoplankton Nutrient Limitation in Time")

            #2D colorplot of Phytoplankton Nutrient Limitation
            XNL = inclass.days
            YNL = inclass.xp
            XNL, YNL = np.meshgrid(XNL, YNL)
            ZNL = inclass.nutlim[:,:,k]
            fignl = plt.figure()
            plt.pcolormesh(XNL,YNL,ZNL,vmin=0.0,vmax=1.0)
            plt.xlabel('Time (days)')
            plt.ylabel('Size (micro-meters)')
            plt.axis([0,inclass.io.simulation_length,inclass.io.phyto_min_size,inclass.io.phyto_max_size])
            plt.title('Layer ' + str(k) + 'P Nutrient Limitation')
            pbar = plt.colorbar()
            pbar.set_label('Nutrient Limitation')
            fignl.savefig(inclass.io.outdir + '/SPLAS_P_LIMIT/' + inclass.io.name + '_' + 'P_N_Limit_2D_layer'+ str(k) +'_.png',dpi=250)
            #plt.show()
            plt.clf()

            print("Color-plot of Phytoplankton Grazing Limitation in Time")

            #2D colorplot of Phytoplankton Grazing Limitation
            XGL = inclass.days
            YGL = inclass.xp
            XGL, YGL = np.meshgrid(XGL, YGL)
            ZGL = inclass.grazlim[:,:,k]
            figgl = plt.figure()
            plt.pcolormesh(XGL,YGL,ZGL,vmin=0.0,vmax=np.amax(ZGL))
            plt.xlabel('Time (days)')
            plt.ylabel('Size (micro-meters)')
            plt.axis([0,inclass.io.simulation_length,inclass.io.phyto_min_size,inclass.io.phyto_max_size])
            plt.title('Layer ' + str(k) + 'P Grazing Limitation')
            pbar = plt.colorbar()
            pbar.set_label('Grazing Limitation')
            figgl.savefig(inclass.io.outdir + '/SPLAS_P_LIMIT/' + inclass.io.name + '_' + 'P_G_Limit_2D_layer'+ str(k) +'_.png',dpi=250)
            #plt.show()
            plt.clf()

            #DIRECTLY IN PYTHON, CALCULATION OF TIME AVERAGE NUTRIENT AND GRAZING LIMITATIONS
            #phytoplankton mean grazing and nutrient limitations
            plt.plot(inclass.xp,inclass.aven[:,k],'lime',inclass.xp,inclass.aveg[:,k],'m--')
            plt.grid()
            plt.axis()
            plt.xlim([inclass.io.phyto_min_size,inclass.io.phyto_max_size])
            plt.xlabel('Size (micro-meters)')
            plt.ylabel('Time Mean Limitation')
            plt.title('Layer ' + str(k) + 'Phytoplankton Limitation')
            fontang = FontProperties()
            fontang.set_size('small')
            plt.legend(['Nutrient','Grazing'], loc = 2, prop = fontang, bbox_to_anchor = (-0.01, 1.12), ncol = 2)
            plt.savefig(inclass.io.outdir + '/SPLAS_P_LIMIT/' + inclass.io.name + '_' + 'P_avetNGL_layer'+ str(k) +'_.png', dpi = 250)
            #plt.show()
            plt.clf()

        ################### BE CAREFUL ############################## BE CAREFUL ################

            print("Plot of Phytoplankton Diversity in time")

            #phytoplankton shannon index vs time
            PSI = inclass.phyto_biodiv[:,k]
            plt.plot(inclass.days,PSI)
            plt.grid()
            plt.axis()
            plt.xlim([0,inclass.io.simulation_length])
            plt.xlabel('Time (days)')
            plt.ylabel('Shannon Evenness')
            plt.title('Layer ' + str(k) + 'Phytoplankton Diversity in Time')
            #plt.legend(['PSI'], loc = 2, bbox_to_anchor = (-0.02, 1.02))
            plt.savefig(inclass.io.outdir + '/SPLAS_DIVERSITY/' + inclass.io.name + '_' + 'PSE_time_layer'+ str(k) +'_.png', dpi = 250)
            #plt.show()
            plt.clf()

            print("Plot of Zooplankton Diversity in time")

            #zooplankton shannon index vs time
            ZSI = inclass.zoo_biodiv[:,k]
            plt.plot(inclass.days,ZSI)
            plt.grid()
            plt.axis()
            plt.xlim([0,inclass.io.simulation_length])
            plt.xlabel('Time (days)')
            plt.ylabel('Shannon Evenness')
            plt.title('Layer ' + str(k) + 'Zooplankton Diversity in Time')
            #plt.legend(['ZSI'], loc = 2, bbox_to_anchor = (-0.02, 1.02))
            plt.savefig(inclass.io.outdir + '/SPLAS_DIVERSITY/' + inclass.io.name + '_' + 'ZSE_time_layer'+ str(k) +'_.png', dpi = 250)
            #plt.show()
            plt.clf()

            XP = inclass.xp
            XZ = inclass.xz

            print("Plot of Phytoplankton Time Mean Biomass Spectrum")

            #phytoplankton time mean size spectrum plot
            plt.plot(XP,inclass.phyto_mean[:,k],'g')
            #plt.plot(XP,inclass.phyto_mean[:,k],'k--')
            plt.grid()
            plt.axis()
            plt.xlim([0,inclass.io.phyto_max_size+1])
            plt.xlabel('Size (micro-meters)')
            plt.ylabel('Time Mean Biomass (uMN)')
            plt.title('Layer ' + str(k) + 'Phytoplankton Time Mean Size Spectrum')
            #plt.legend(['< P >'], loc = 2, bbox_to_anchor = (-0.02, 1.02))
            plt.savefig(inclass.io.outdir + '/SPLAS_SIZE_SPECTRA/' + inclass.io.name + '_' + 'P_Time_Mean_layer'+ str(k) +'_.png', dpi = 250)
            #plt.show()
            plt.clf()

            print("Plot of Zooplankton Time Mean Biomass Spectrum")

            #zooplankton time mean size spectrum plot
            plt.plot(XZ,inclass.zoo_mean[:,k],'r')
            #plt.plot(XZ,inclass.zoo_mean[:,k],'k--')
            plt.grid()
            plt.axis()
            plt.xlim([0,inclass.io.zoo_max_size+2])
            plt.xlabel('Size (micro-meters)')
            plt.ylabel('Time Mean Biomass (uMN)')
            plt.title('Layer ' + str(k) + 'Zooplankton Time Mean Size Spectrum')
            #plt.legend(['< Z >'], loc = 2, bbox_to_anchor = (-0.02, 1.02))
            plt.savefig(inclass.io.outdir + '/SPLAS_SIZE_SPECTRA/' + inclass.io.name + '_' + 'Z_Time_Mean_layer'+ str(k) +'_.png', dpi = 250)
            #plt.show()
            plt.clf()

            print("Plot of Phytoplankton Total Biomass in time")

            #phytoplankton total biomass vs time
            PTB = inclass.phyto_biomass[:,k]
            plt.plot(inclass.days,PTB,color='g')
            plt.grid()
            plt.axis()
            plt.axis([0,inclass.io.simulation_length,0,np.amax(PTB)+0.3])
            plt.xlabel('Time (days)')
            plt.ylabel('Total Biomass (uMN)')
            plt.title('Layer ' + str(k) + 'Phytoplankton Biomass in Time')
            #plt.legend(['PTB'], loc = 2, bbox_to_anchor = (-0.01, 1.1))
            plt.savefig(inclass.io.outdir + '/SPLAS_BIOMASS/' + inclass.io.name + '_' + 'P_total_t_layer'+ str(k) +'_.png', dpi = 250)
            #plt.show()
            plt.clf()

            print("Plot of Zooplankton Total Biomass in time")

            #zooplankton total biomass vs time
            ZTB = inclass.zoo_biomass[:,k]
            plt.plot(inclass.days,ZTB,color='r')
            plt.grid()
            plt.axis()
            plt.axis([0,inclass.io.simulation_length,0,np.amax(ZTB)+0.3])
            plt.xlabel('Time (days)')
            plt.ylabel('Total Biomass (uMN)')
            plt.title('Layer ' + str(k) + 'Zooplankton Biomass in Time')
            #plt.legend(['ZTB'], loc = 2, bbox_to_anchor = (-0.01, 1.1))
            plt.savefig(inclass.io.outdir + '/SPLAS_BIOMASS/' + inclass.io.name + '_' + 'Z_total_t_layer'+ str(k) +'_.png', dpi = 250)
            #plt.show()
            plt.clf()

            print("Plot of Zooplankton to Phytoplankton Biomass Ratio in time")

            #zooplankton over phytoplankton total biomass vs time
            R = inclass.z_p_ratio[:,k]
            plt.plot(inclass.days,R,color='orange')
            plt.grid()
            plt.axis()
            plt.axis([0,inclass.io.simulation_length,0,np.amax(R)+0.1])
            plt.xlabel('Time (days)')
            plt.ylabel('Z:P biomass ratio')
            plt.title('Layer ' + str(k) + 'Z/P Ratio in Time')
            #plt.legend(['R'], loc = 2, bbox_to_anchor = (-0.01, 1.1))
            plt.savefig(inclass.io.outdir + '/SPLAS_BIOMASS/' + inclass.io.name + '_' + 'z_p_ratio_layer'+ str(k) +'_.png', dpi = 250)
            #plt.show()
            plt.clf()

        ################################################################################

    plt.close('all')

####################################################################################################################################

def phyto_details(inclass):
#now will depend on layer! (DEPTH)

    for k in range(inclass.io.m):
            X = inclass.xp

            #phytoplankton cell quota final time plot
            plt.plot(X,inclass.fpquota[:,k],'go')
            plt.plot(X,inclass.fpquota[:,k],'k--')
            plt.grid()
            plt.axis()
            plt.xlim([0,inclass.io.phyto_max_size+1])
            plt.xlabel('Size (micro-meters)')
            #plt.ylabel('Internal Quota (micro-moles of N cell-1)')
            plt.title('Layer ' + str(k) + 'Cell Quota (mmolN cell-1)')
            #plt.title('Phytoplankton Size Spectrum')
            plt.savefig(inclass.io.outdir + '/SPLAS_P_DYN/' + inclass.io.name + '_' + 'P_fin_QUOTA_layer' + str(k) + '_.png', dpi = 250)
            #plt.show()
            plt.clf()

            #phytoplankton number of cells final time plot
            plt.plot(np.log10(X),np.log10(inclass.fpncells[:,k]),'ro')
            plt.plot(np.log10(X),np.log10(inclass.fpncells[:,k]),'k--')
            plt.grid()
            plt.axis()
            plt.xlim([np.log10(inclass.io.phyto_min_size),np.log10(inclass.io.phyto_max_size)])
            plt.xlabel('Log10 of Size (micro-meters)')
            plt.ylabel('Log10 of Density of Cells (cells m-3)')
            plt.title('Layer ' + str(k) + 'P Size Spectrum')
            plt.savefig(inclass.io.outdir + '/SPLAS_P_DYN/' + inclass.io.name + '_' + 'P_fin_log_CELLS_layer' + str(k) + '_.png', dpi = 250)
            #plt.show()
            plt.clf()

            plt.plot(X,inclass.fpncells[:,k],'ro')
            plt.plot(X,inclass.fpncells[:,k],'k--')
            plt.grid()
            plt.axis()
            plt.xlim([inclass.io.phyto_min_size,inclass.io.phyto_max_size])
            plt.xlabel('Size (micro-meters)')
            plt.ylabel('Density of Cells (cells m-3)')
            plt.title('Layer ' + str(k) + 'P Size Spectrum')
            plt.savefig(inclass.io.outdir + '/SPLAS_P_DYN/' + inclass.io.name + '_' + 'P_fin_CELLS_layer' + str(k) + '_.png', dpi = 250)
            #plt.show()
            plt.clf()

            #phytoplankton uptake and growth

            plt.plot(X,inclass.p_growth_rate[:,k],'lime')
            plt.grid()
            plt.axis()
            plt.xlim([0,inclass.io.phyto_max_size+1])
            plt.xlabel('Size (micro-meters)')
            plt.ylabel('Growth Rate (day-1)')
            plt.title('Layer ' + str(k) + 'Growth Rate')
            plt.savefig(inclass.io.outdir + '/SPLAS_P_DYN/' + inclass.io.name + '_' + 'P_fin_GROWTH_layer' + str(k) + '_.png', dpi = 250)
            #plt.show()
            plt.clf()

            plt.plot(X,inclass.p_vmax*inclass.nit_uptake[:,k],'m--')
            plt.grid()
            plt.axis()
            plt.xlim([0,inclass.io.phyto_max_size+1])
            plt.xlabel('Size (micro-meters)')
            plt.ylabel('Uptake Rate (mmolN day-1 cell-1)')
            plt.title('Layer ' + str(k) + ' Uptake Rate')
            plt.savefig(inclass.io.outdir + '/SPLAS_P_DYN/' + inclass.io.name + '_' + 'P_fin_UPTAKE_layer' + str(k) + '_.png', dpi = 250)
            #plt.show()
            plt.clf()

        ################################################################################

    plt.close('all')

####################################################################################################################################

def print_output(inclass):

    print("   ")
    print("simulation " + inclass.io.name)
    print("output files in " + inclass.io.outdir)
    print("   ")

    print("Average picophytoplankton biomass (uMN)")
    inclass.avez_small_phyto = (1/inclass.h)*inclass.io.dz*inclass.small_phyto.sum()
    print("Average phytoplankton biomass (uMN)")
    inclass.avez_p_biomass = (1/inclass.h)*inclass.io.dz*inclass.phyto_total_mean_biomass.sum()
    print("Average zooplankton biomass (uMN)")
    inclass.avez_z_biomass = (1/inclass.h)*inclass.io.dz*inclass.zoo_total_mean_biomass.sum()
    print("Average Shannon Evenness (P,Z)")
    inclass.avez_shannon_p = (1/inclass.h)*inclass.io.dz*inclass.mean_shannon_p.sum()
    inclass.avez_shannon_z = (1/inclass.h)*inclass.io.dz*inclass.mean_shannon_z.sum()
    print("Average N concentration (uMN)")
    inclass.avez_N = (1/inclass.h)*inclass.io.dz*inclass.mean_N.sum()
    print("Mean final N concentration (uMN)")
    inclass.avez_z_p = (1/inclass.h)*inclass.io.dz*inclass.ave_z_p_ratio.sum()

    #print(" PHYTOPLANKTON WINNERS ")
    #idx = np.where(inclass.phyto_mean >= 0.1)
    #np.concatenate(idx)
    #print(idx)
    #print(inclass.xp[idx])
    #print(" ZOOPLANKTON WINNERS ")
    #idx2 = np.where(inclass.zoo_mean >= 0.01)
    #np.concatenate(idx2)
    #print(idx2)
    #print(inclass.xz[idx2])

    print("PHYTOPLANKTON FINAL CARBON EXPORT (mgC m-2 day-1)")
    print(inclass.p_total_carbon_export[-1])
    print("DETRITUS FINAL CARBON EXPORT (mgC m-2 day-1)")
    print(inclass.z_carbon_export[-1])
    print("TOTAL FINAL CARBON EXPORT (mgC m-2 day-1)")
    print(inclass.total_carbon_export[-1])

    print("PHYTOPLANKTON AVERAGE CARBON EXPORT (mgC m-2 day-1)")
    print(inclass.phyto_ave_total_c_exp)
    print("DETRITUS AVERAGE CARBON EXPORT (mgC m-2 day-1)")
    print(inclass.zoo_ave_total_c_exp)
    print("TOTAL AVERAGE CARBON EXPORT (mgC m-2 day-1)")
    print(inclass.total_ave_carbon_export)

    print("FINAL DETRITUS")
    print(inclass.final_detritus)
    print("MEAN DETRITUS")
    print(inclass.mean_detritus)

   # print("PHYTOPLANKTON INTEGRATED CARBON EXPORT (mgC m-2)")
   # print(inclass.phyto_integrated_total_c_exp)
   # print("DETRITUS INTEGRATED CARBON EXPORT (mgC m-2)")
   # print(inclass.zoo_integrated_total_c_exp)
   # print("TOTAL INTEGRATED CARBON EXPORT (mgC m-2)")
   # print(inclass.total_integrated_carbon_export)

####################################################################################################################################

def get_depth_profiles(inclass):
#depth profiles of quantities

#invert y axis to show depth profile
    #ax = plt.gca()
    #ax.invert_yaxis()
    #plt.gca().invert_yaxis()
#simplest method, use matplotlib.plot.ylim
    #------> plt.ylim(max,min)
    #plt.ylim(plt.ylim()[::-1])

    #plt.title('',y=1.08)

################################################################################

#2-D colorplots
################################################################################################################################################################
#plankton biomass vs time and depth, for each size-class

    for i in range(inclass.xp.size):
            Xrz = inclass.days
            Yrz = inclass.z
            Xrz, Yrz = np.meshgrid(Xrz, Yrz)
            Z = inclass.phyto[i,:,:]/np.amax(inclass.phyto[i,:,:],0)
            Zrz = Z.transpose()
            figrz = plt.figure()
            plt.pcolormesh(Xrz,Yrz,Zrz,vmin=0.0,vmax=1.0)
            plt.xlabel('Time (days)')
            plt.ylabel('Depth (m)')
            plt.axis([0,inclass.io.simulation_length,inclass.z[-1],inclass.z[0]])
            plt.title('Size-class ' + str(i) + 'P Relative Biomass')
            zrbar = plt.colorbar()
            zrbar.set_label('Relative Biomass')
            figrz.savefig(inclass.io.outdir + '/DEPTH_P/' + inclass.io.name + '_' + 'P_depth_size'+ str(i) +'_.png',dpi=250)
            #plt.show()
            plt.clf()
    plt.close()

    for i in range(inclass.xz.size):
            Xrz = inclass.days
            Yrz = inclass.z
            Xrz, Yrz = np.meshgrid(Xrz, Yrz)
            Z = inclass.zoo[i,:,:]/np.amax(inclass.zoo[i,:,:],0)
            Zrz = Z.transpose()
            figrz = plt.figure()
            plt.pcolormesh(Xrz,Yrz,Zrz,vmin=0.0,vmax=1.0)
            plt.xlabel('Time (days)')
            plt.ylabel('Depth (m)')
            plt.axis([0,inclass.io.simulation_length,inclass.z[-1],inclass.z[0]])
            plt.title('Size-class ' + str(i) + 'Z Relative Biomass')
            zrbar = plt.colorbar()
            zrbar.set_label('Relative Biomass')
            figrz.savefig(inclass.io.outdir + '/DEPTH_Z/' + inclass.io.name + '_' + 'Z_depth_size'+ str(i) +'_.png',dpi=250)
            #plt.show()
            plt.clf()
    plt.close()

################################################################################
#time-mean and final size-spectrum vs depth

#irradiance fil
    print("Irradiance Profile")
    Xrz = inclass.days
    Yrz = inclass.z
    Xrz, Yrz = np.meshgrid (Xrz, Yrz)
    I = inclass.irrprofile[:,:]
    Zrz = I.transpose()
    figrz = plt.figure()
    plt.pcolormesh(Xrz,Yrz,Zrz,vmin=0.0,vmax=np.amax(I))
    plt.ylabel('Depth (m)')
    plt.xlabel('Time (days)')
    plt.axis([0,inclass.io.simulation_length,inclass.z[-1],inclass.z[0]])
    plt.title('Irradiance_profile')
    zrbar = plt.colorbar()
    zrbar.set_label ('Irradiance W/m^2')
    figrz.savefig(inclass.io.outdir + '/SPLAS_IRRADIANCE/' + inclass.io.name + '_' + 'Irradiance Profile' + '_.png', dpi=250)
    plt.clf()
    plt.close()
#irradiance fil
    print("Irradiance Profile")
    Xrz = inclass.days
    Yrz = inclass.z
    Xrz, Yrz = np.meshgrid (Xrz, Yrz)
    I = inclass.irrprofile[:,:]
    Zrz = I.transpose()
    figrz = plt.figure()
    plt.pcolormesh(Xrz,Yrz,Zrz)
    plt.ylabel('Depth (m)')
    plt.xlabel('Time (days)')
    plt.axis([0,inclass.io.simulation_length,inclass.z[-1],inclass.z[0]])
    plt.title('Irradiance_profile')
    zrbar = plt.colorbar()
    zrbar.set_label ('Irradiance W/m^2')
    figrz.savefig(inclass.io.outdir + '/SPLAS_IRRADIANCE/' + inclass.io.name + '_' + 'Irradiance Profile2' + '_.png', dpi=250)
    plt.clf()
    plt.close()
#irradiance fil
    print("Ln Irradiance Profile")
    Xrz = inclass.days
    Yrz = inclass.z
    Xrz, Yrz = np.meshgrid (Xrz, Yrz)
    I = inclass.irrprofile[:,:]
    Zrz = I.transpose()
    figrz = plt.figure()
    plt.pcolormesh(Xrz,Yrz,Zrz,vmin=0.0,vmax=np.amax(np.log(I+1)))
    plt.ylabel('Depth (m)')
    plt.xlabel('Time (days)')
    plt.axis([0,inclass.io.simulation_length,inclass.z[-1],inclass.z[0]])
    plt.title('Irradiance_profile')
    zrbar = plt.colorbar()
    zrbar.set_label ('Irradiance ln(W/m^2)')
    figrz.savefig(inclass.io.outdir + '/SPLAS_IRRADIANCE/' + inclass.io.name + '_' + 'Irradiance Profile4' + '_.png', dpi=250)
    plt.clf()
    plt.close()
#irradiance fil
    print("sqrt Irradiance Profile")
    Xrz = inclass.days
    Yrz = inclass.z
    Xrz, Yrz = np.meshgrid (Xrz, Yrz)
    I = inclass.irrprofile[:,:]
    Zrz = I.transpose()
    figrz = plt.figure()
    plt.pcolormesh(Xrz,Yrz,Zrz,vmin=0.0,vmax=np.amax(np.sqrt(I)),shading='auto')
    plt.ylabel('Depth (m)')
    plt.xlabel('Time (days)')
    plt.axis([0,inclass.io.simulation_length,inclass.z[-1],inclass.z[0]])
    plt.title('Irradiance_profile')
    zrbar = plt.colorbar()
    zrbar.set_label ('Irradiance strength sqrt(W/m^2)')
    figrz.savefig(inclass.io.outdir + '/SPLAS_IRRADIANCE/' + inclass.io.name + '_' + 'Irradiance Profile3' + '_.png', dpi=250)
    plt.clf()
    plt.close()
#####################################################################################################
    Xrz = inclass.xp
    Yrz = inclass.z
    Xrz, Yrz = np.meshgrid(Xrz, Yrz)
    Z = inclass.phyto_mean/np.amax(inclass.phyto_mean,0)
    Zrz = Z.transpose()
    figrz = plt.figure()
    plt.pcolormesh(Xrz,Yrz,Zrz,vmin=0.0,vmax=1.0)
    plt.ylabel('Depth (m)')
    plt.xlabel('Size (micro-meters)')
    plt.axis([inclass.io.phyto_min_size,inclass.io.phyto_max_size,inclass.z[-1],inclass.z[0]])
    plt.title('P Time Mean Size Spectrum')
    zrbar = plt.colorbar()
    zrbar.set_label('Relative Biomass')
    figrz.savefig(inclass.io.outdir + '/DEPTH_SIZE_SPECTRA/' + inclass.io.name + '_' + 'P_time_mean_SS.png',dpi=250)
    #plt.show()
    plt.clf()
    plt.close()

    Xrz = inclass.xz
    Yrz = inclass.z
    Xrz, Yrz = np.meshgrid(Xrz, Yrz)
    Z = inclass.zoo_mean/np.amax(inclass.zoo_mean,0)
    Zrz = Z.transpose()
    figrz = plt.figure()
    plt.pcolormesh(Xrz,Yrz,Zrz,vmin=0.0,vmax=1.0)
    plt.ylabel('Depth (m)')
    plt.xlabel('Size (micro-meters)')
    plt.axis([inclass.io.zoo_min_size,inclass.io.zoo_max_size,inclass.z[-1],inclass.z[0]])
    plt.title('Z Time Mean Size Spectrum')
    zrbar = plt.colorbar()
    zrbar.set_label('Relative Biomass')
    figrz.savefig(inclass.io.outdir + '/DEPTH_SIZE_SPECTRA/' + inclass.io.name + '_' + 'Z_time_mean_SS.png',dpi=250)
    #plt.show()
    plt.clf()
    plt.close()

    Xrz = inclass.xp
    Yrz = inclass.z
    Xrz, Yrz = np.meshgrid(Xrz, Yrz)
    Z = inclass.fphyto/np.amax(inclass.fphyto,0)
    Zrz = Z.transpose()
    figrz = plt.figure()
    plt.pcolormesh(Xrz,Yrz,Zrz,vmin=0.0,vmax=1.0)
    plt.ylabel('Depth (m)')
    plt.xlabel('Size (micro-meters)')
    plt.axis([inclass.io.phyto_min_size,inclass.io.phyto_max_size,inclass.z[-1],inclass.z[0]])
    plt.title('P Final Size Spectrum')
    zrbar = plt.colorbar()
    zrbar.set_label('Relative Biomass')
    figrz.savefig(inclass.io.outdir + '/DEPTH_SIZE_SPECTRA/' + inclass.io.name + '_' + 'P_final_SS.png',dpi=250)
    #plt.show()
    plt.clf()
    plt.close()

    Xrz = inclass.xz
    Yrz = inclass.z
    Xrz, Yrz = np.meshgrid(Xrz, Yrz)
    Z = inclass.fzoo/np.amax(inclass.fzoo,0)
    Zrz = Z.transpose()
    figrz = plt.figure()
    plt.pcolormesh(Xrz,Yrz,Zrz,vmin=0.0,vmax=1.0)
    plt.ylabel('Depth (m)')
    plt.xlabel('Size (micro-meters)')
    plt.axis([inclass.io.zoo_min_size,inclass.io.zoo_max_size,inclass.z[-1],inclass.z[0]])
    plt.title('Z Final Size Spectrum')
    zrbar = plt.colorbar()
    zrbar.set_label('Relative Biomass')
    figrz.savefig(inclass.io.outdir + '/DEPTH_SIZE_SPECTRA/' + inclass.io.name + '_' + 'Z_final_SS.png',dpi=250)
    #plt.show()
    plt.clf()
    plt.close()

################################################################################
#Nitrogen concentration in time vs depth

    Xrz = inclass.days
    Yrz = inclass.z
    Xrz, Yrz = np.meshgrid(Xrz, Yrz)
    Z = inclass.Nitrogen/np.amax(inclass.Nitrogen,0)
    Zrz = Z.transpose()
    figrz = plt.figure()
    plt.pcolormesh(Xrz,Yrz,Zrz,vmin=0.0,vmax=1.0)
    plt.xlabel('Time (days)')
    plt.ylabel('Depth (m)')
    plt.axis([0,inclass.io.simulation_length,inclass.z[-1],inclass.z[0]])
    plt.title('N(t) depth profile')
    zrbar = plt.colorbar()
    zrbar.set_label('Relative Concentration')
    figrz.savefig(inclass.io.outdir + '/DEPTH_N/' + inclass.io.name + '_' + 'N_time_depth_profile.png',dpi=250)
    #plt.show()
    plt.clf()
    plt.close()

################################################################################

    #inclass.fpquota = inclass.pquota[:,-1,:]
    #inclass.fpncells = (inclass.fphyto/inclass.fpquota)

################################################################################################################################################################

#1-D profiles
################################################################################

    #Nitrogen time mean and final concentration ################

    N = inclass.mean_N
    plt.plot(N,inclass.z,color='lightcoral')
    plt.grid()
    plt.axis()
    plt.ylim([inclass.z[-1],inclass.z[0]])
    plt.xlabel('N mean concentration (uMN)')
    plt.ylabel('Depth (m)')
    plt.title('[N] mean profile')
    #plt.legend(['N'], loc = 2, bbox_to_anchor = (-0.02, 1.02))
    plt.savefig(inclass.io.outdir + '/DEPTH_PROFILE/' + inclass.io.name + '_' + 'mean_N.png', dpi = 250)
    #plt.show()
    plt.clf()
    plt.close()

    N = inclass.final_N
    plt.plot(N,inclass.z,'m')
    plt.grid()
    plt.axis()
    plt.ylim([inclass.z[-1],inclass.z[0]])
    plt.xlabel('N final concentration (uMN)')
    plt.ylabel('Depth (m)')
    plt.title('[N] final profile')
    #plt.legend(['N'], loc = 2, bbox_to_anchor = (-0.02, 1.02))
    plt.savefig(inclass.io.outdir + '/DEPTH_PROFILE/' + inclass.io.name + '_' + 'final_N.png', dpi = 250)
    #plt.show()
    plt.clf()
    plt.close()

    #plankton biomass #########################################

    ZP = inclass.ave_z_p_ratio
    plt.plot(ZP,inclass.z,'orange')
    plt.grid()
    plt.axis()
    plt.ylim([inclass.z[-1],inclass.z[0]])
    plt.xlabel('Mean Z:P Ratio')
    plt.ylabel('Depth (m)')
    plt.title('Mean Z:P Ratio profile')
    #plt.legend(['N'], loc = 2, bbox_to_anchor = (-0.02, 1.02))
    plt.savefig(inclass.io.outdir + '/DEPTH_PROFILE/' + inclass.io.name + '_' + 'mean_Z_P_ratio.png', dpi = 250)
    #plt.show()
    plt.clf()
    plt.close()

    small = inclass.small_phyto/inclass.phyto_total_mean_biomass
    plt.plot(small,inclass.z,'lime')
    plt.grid()
    plt.axis()
    plt.ylim([inclass.z[-1],inclass.z[0]])
    plt.xlabel('pico:total P ratio')
    plt.ylabel('Depth (m)')
    plt.title('pico:total P ratio profile')
    #plt.legend(['N'], loc = 2, bbox_to_anchor = (-0.02, 1.02))
    plt.savefig(inclass.io.outdir + '/DEPTH_PROFILE/' + inclass.io.name + '_' + 'small_phyto.png', dpi = 250)
    #plt.show()
    plt.clf()
    plt.close()

    p = inclass.phyto_total_mean_biomass
    plt.plot(p,inclass.z,'g')
    plt.grid()
    plt.axis()
    plt.ylim([inclass.z[-1],inclass.z[0]])
    plt.xlabel('P Biomass (uMN)')
    plt.ylabel('Depth (m)')
    plt.title('Total P Mean Profile')
    #plt.legend(['N'], loc = 2, bbox_to_anchor = (-0.02, 1.02))
    plt.savefig(inclass.io.outdir + '/DEPTH_PROFILE/' + inclass.io.name + '_' + 'tot_phyto.png', dpi = 250)
    #plt.show()
    plt.clf()
    plt.close()

    zz = inclass.zoo_total_mean_biomass
    plt.plot(zz,inclass.z,'r')
    plt.grid()
    plt.axis()
    plt.ylim([inclass.z[-1],inclass.z[0]])
    plt.xlabel('Z Biomass (uMN)')
    plt.ylabel('Depth (m)')
    plt.title('Total Z Mean Profile')
    #plt.legend(['N'], loc = 2, bbox_to_anchor = (-0.02, 1.02))
    plt.savefig(inclass.io.outdir + '/DEPTH_PROFILE/' + inclass.io.name + '_' + 'tot_zoo.png', dpi = 250)
    #plt.show()
    plt.clf()
    plt.close()

    ###################################DETRITUS#################################

    #detritus time mean and final concentration ################

    mdet = inclass.mean_detritus
    plt.plot(mdet,inclass.z,color='orange')
    plt.grid()
    plt.axis()
    plt.ylim([inclass.z[-1],inclass.z[0]])
    plt.xlabel('Detritus mean concentration (uMN)')
    plt.ylabel('Depth (m)')
    plt.title('[Detritus] mean profile')
    #plt.legend(['N'], loc = 2, bbox_to_anchor = (-0.02, 1.02))
    plt.savefig(inclass.io.outdir + '/DEPTH_PROFILE/' + inclass.io.name + '_' + 'mean_DETRITUS.png', dpi = 250)
    #plt.show()
    plt.clf()
    plt.close()

    fdet = inclass.final_detritus
    plt.plot(fdet,inclass.z,'r')
    plt.grid()
    plt.axis()
    plt.ylim([inclass.z[-1],inclass.z[0]])
    plt.xlabel('Detritus final concentration (uMN)')
    plt.ylabel('Depth (m)')
    plt.title('[Detritus] final profile')
    #plt.legend(['N'], loc = 2, bbox_to_anchor = (-0.02, 1.02))
    plt.savefig(inclass.io.outdir + '/DEPTH_PROFILE/' + inclass.io.name + '_' + 'final_DETRITUS.png', dpi = 250)
    #plt.show()
    plt.clf()
    plt.close()

    ############################################################################

    #shannon evenness #########################################

    PSE = inclass.mean_shannon_p
    plt.plot(PSE,inclass.z,'b')
    plt.grid()
    plt.axis()
    plt.ylim([inclass.z[-1],inclass.z[0]])
    plt.xlabel('P Shannon Evenness')
    plt.ylabel('Depth (m)')
    plt.title('P Diversity Profile')
    #plt.legend(['N'], loc = 2, bbox_to_anchor = (-0.02, 1.02))
    plt.savefig(inclass.io.outdir + '/DEPTH_PROFILE/' + inclass.io.name + '_' + 'PSE.png', dpi = 250)
    #plt.show()
    plt.clf()
    plt.close()

    ZSE = inclass.mean_shannon_z
    plt.plot(ZSE,inclass.z,'k')
    plt.grid()
    plt.axis()
    plt.ylim([inclass.z[-1],inclass.z[0]])
    plt.xlabel('Z Shannon Evenness')
    plt.ylabel('Depth (m)')
    plt.title('Z Diversity Profile')
    #plt.legend(['N'], loc = 2, bbox_to_anchor = (-0.02, 1.02))
    plt.savefig(inclass.io.outdir + '/DEPTH_PROFILE/' + inclass.io.name + '_' + 'ZSE.png', dpi = 250)
    #plt.show()
    plt.clf()
    plt.close()

    #diffusion profile #########################################

    plt.plot(inclass.diffusion_profile,inclass.z,'r--')
    plt.grid()
    plt.axis()
    plt.ylim([inclass.z[-1],inclass.z[0]])
    plt.xlabel('Diffusion Coefficient (m2/s)')
    plt.ylabel('Depth (m)')
    plt.title('Diffusion Coefficient Profile')
    #plt.legend(['N'], loc = 2, bbox_to_anchor = (-0.02, 1.02))
    plt.savefig(inclass.io.outdir + '/DEPTH_PROFILE/' + inclass.io.name + '_' + 'DIFF.png', dpi = 250)
    #plt.show()
    plt.clf()
    plt.close()

#####################################################################################################################################################
#final and time mean size spectra all in one plot (4 plots --> p and z, final and mean)

    X = inclass.xp
    Y = inclass.xz

    #################################
    plt.rc('xtick',labelsize=5)
    plt.rc('ytick',labelsize=5)
    #################################

    fig, ax = plt.subplots(nrows = inclass.io.m-2, ncols = 1, sharex = True, sharey = True)

    for k in range(1,inclass.io.m-1):
        plt.subplot(inclass.io.m-2,1,k)
        plt.plot(X,inclass.fphyto[:,k],'g')
        plt.grid()
        plt.axis()
        plt.xlim([0,inclass.io.phyto_max_size+1])
        plt.title('Layer ' + str(k) + 'P Size Spectrum', fontsize=5)
        plt.subplots_adjust(hspace=0.55)

    fig.text(0.5, 0.04, 'Size (um)', ha='center', va='center')
    fig.text(0.06, 0.5, 'Biomass (uMN)', ha='center', va='center', rotation='vertical')
    plt.savefig(inclass.io.outdir + '/DEPTH_PROFILE/' + inclass.io.name + '_' + 'P_final_ss_depth_.png', dpi = 250)
    plt.clf()
    plt.close()

    #################################

    fig, ax = plt.subplots(nrows = inclass.io.m-2, ncols = 1, sharex = True, sharey = True)

    for k in range(1,inclass.io.m-1):
        plt.subplot(inclass.io.m-2,1,k)
        plt.plot(Y,inclass.fzoo[:,k],'r')
        plt.grid()
        plt.axis()
        plt.xlim([0,inclass.io.zoo_max_size+2])
        plt.title('Layer ' + str(k) + 'Z Size Spectrum', fontsize=5)
        plt.subplots_adjust(hspace=0.35)

    fig.text(0.5, 0.04, 'Size (um)', ha='center', va='center')
    fig.text(0.06, 0.5, 'Biomass (uMN)', ha='center', va='center', rotation='vertical')
    plt.savefig(inclass.io.outdir + '/DEPTH_PROFILE/' + inclass.io.name + '_' + 'Z_final_ss_depth_.png', dpi = 250)
    plt.clf()
    plt.close()

    #################################

    fig, ax = plt.subplots(nrows = inclass.io.m-2, ncols = 1, sharex = True, sharey = True)

    for k in range(1,inclass.io.m-1):
        plt.subplot(inclass.io.m-2,1,k)
        plt.plot(X,inclass.phyto_mean[:,k],'g')
        plt.grid()
        plt.axis()
        plt.xlim([0,inclass.io.phyto_max_size+1])
        plt.title('Layer ' + str(k) + 'P Mean Size Spectrum', fontsize=5)
        plt.subplots_adjust(hspace=0.55)

    fig.text(0.5, 0.04, 'Size (um)', ha='center', va='center')
    fig.text(0.06, 0.5, 'Biomass (uMN)', ha='center', va='center', rotation='vertical')
    plt.savefig(inclass.io.outdir + '/DEPTH_PROFILE/' + inclass.io.name + '_' + 'P_mean_ss_depth_.png', dpi = 250)
    plt.clf()
    plt.close()

    #################################

    fig, ax = plt.subplots(nrows = inclass.io.m-2, ncols = 1, sharex = True, sharey = True)

    for k in range(1,inclass.io.m-1):
        plt.subplot(inclass.io.m-2,1,k)
        plt.plot(Y,inclass.zoo_mean[:,k],'r')
        plt.grid()
        plt.axis()
        plt.xlim([0,inclass.io.zoo_max_size+2])
        plt.title('Layer ' + str(k) + 'Z Mean Size Spectrum', fontsize=5)
        plt.subplots_adjust(hspace=0.35)

    fig.text(0.5, 0.04, 'Size (um)', ha='center', va='center')
    fig.text(0.06, 0.5, 'Biomass (uMN)', ha='center', va='center', rotation='vertical')
    plt.savefig(inclass.io.outdir + '/DEPTH_PROFILE/' + inclass.io.name + '_' + 'Z_mean_ss_depth_.png', dpi = 250)
    plt.clf()
    plt.close()

    #################################

################################################################################

    fig, ax = plt.subplots(nrows = inclass.io.m-1, ncols = 1, sharex = True, sharey = True)

    for k in range(1,inclass.io.m):
        plt.subplot(inclass.io.m-1,1,k)
        plt.plot(X,(inclass.fphyto[:,k]-inclass.fphyto[:,k-1]),'b')
        plt.grid()
        plt.axis()
        plt.xlim([0,inclass.io.phyto_max_size+1])
        plt.title('Layer ' + str(k) + ' - Layer ' + str(k-1) + ' P mass', fontsize=5)
        plt.subplots_adjust(hspace=0.55)

    fig.text(0.5, 0.04, 'Size (um)', ha='center', va='center')
    fig.text(0.06, 0.5, 'Biomass difference (uMN)', ha='center', va='center', rotation='vertical')
    plt.savefig(inclass.io.outdir + '/DEPTH_PROFILE/' + inclass.io.name + '_' + 'DIFFERENCE_P_ss_depth_.png', dpi = 250)
    plt.clf()
    plt.close()

    for k in range(1,inclass.io.m):
        plt.subplot(inclass.io.m-1,1,k)
        plt.plot(Y,(inclass.fzoo[:,k]-inclass.fzoo[:,k-1]),'m')
        plt.grid()
        plt.axis()
        plt.xlim([0,inclass.io.zoo_max_size+2])
        plt.title('Layer ' + str(k) + ' - Layer ' + str(k-1) + ' Z mass', fontsize=5)
        plt.subplots_adjust(hspace=0.35)

    fig.text(0.5, 0.04, 'Size (um)', ha='center', va='center')
    fig.text(0.06, 0.5, 'Biomass difference (uMN)', ha='center', va='center', rotation='vertical')
    plt.savefig(inclass.io.outdir + '/DEPTH_PROFILE/' + inclass.io.name + '_' + 'DIFFERENCE_Z_ss_depth_.png', dpi = 250)
    plt.clf()
    plt.close()

################################################################################

    plt.close('all')

####################################################################################################################################

def get_carbon_export(inclass):

    matplotlib.rcParams.update(matplotlib.rcParamsDefault)
################################################################################

    print("Plot of average phytoplankton carbon export")

    #plot of average phytoplankton carbon export for each size class
    X = inclass.xp
    plt.plot(X,inclass.phyto_ave_c_exp,'yo')
    plt.plot(X,inclass.phyto_ave_c_exp,'k--')
    plt.grid()
    plt.axis()
    plt.xlim([0,inclass.io.phyto_max_size+1])
    plt.xlabel('Size (micro-meters)')
    plt.ylabel('Mean C Export (mgC m-2 day-1)')
    plt.title('Phytoplankton Mean Carbon Export')
    plt.savefig(inclass.io.outdir + '/SPLAS_CARBON/' + inclass.io.name + '_' + 'P_MEAN_CARBON_EXPORT.png', dpi = 250)
    #plt.show()
    plt.clf()

    print("Plot of integrated phytoplankton carbon export")

    #plot of integrated phytoplankton carbon export for each size class
    plt.plot(X,inclass.phyto_integrated_c_exp,'co')
    plt.plot(X,inclass.phyto_integrated_c_exp,'k--')
    plt.grid()
    plt.axis()
    plt.xlim([0,inclass.io.phyto_max_size+1])
    plt.xlabel('Size (micro-meters)')
    plt.ylabel('Total C Export (mgC m-2)')
    plt.title('Phytoplankton Integrated Carbon Export')
    plt.savefig(inclass.io.outdir + '/SPLAS_CARBON/' + inclass.io.name + '_' + 'P_INTEGRATED_CARBON_EXPORT.png', dpi = 250)
    #plt.show()
    plt.clf()

    print("Plot of Phytoplankton total Carbon export in time")

    #total phytoplankton carbon export vs time
    pexpc = inclass.p_total_carbon_export[1:]
    plt.plot(inclass.days[1:],pexpc,color='g')
    plt.grid()
    plt.axis()
    plt.xlim([0,inclass.io.simulation_length])
    plt.xlabel('Time (days)')
    plt.ylabel('Phytoplankton C Export (mgC m-2 day-1)')
    plt.title('Phytoplankton C Export in time')
    plt.savefig(inclass.io.outdir + '/SPLAS_CARBON/' + inclass.io.name + '_' + 'PHYTO_CARBON_EXPORT_time.png', dpi = 250)
    #plt.show()
    plt.clf()

    print("Plot of Detritus Carbon export in time")

    #total detritus carbon export vs time
    zexpc = inclass.z_carbon_export
    plt.plot(inclass.days,zexpc,color='r')
    plt.grid()
    plt.axis()
    plt.xlim([0,inclass.io.simulation_length])
    plt.xlabel('Time (days)')
    plt.ylabel('Detritus C Export (mgC m-2 day-1)')
    plt.title('Detritus C Export in time')
    plt.savefig(inclass.io.outdir + '/SPLAS_CARBON/' + inclass.io.name + '_' + 'DETRITUS_CARBON_EXPORT_time.png', dpi = 250)
    #plt.show()
    plt.clf()

    print("Plot of TOTAL (Sinking Phytoplankton and Detritus) Carbon export in time")

    #total carbon export vs time
    expc = inclass.total_carbon_export[1:]
    plt.plot(inclass.days[1:],expc,color='b')
    plt.grid()
    plt.axis()
    plt.xlim([0,inclass.io.simulation_length])
    plt.ylim([0,200.0])
    plt.xlabel('Time (days)')
    plt.ylabel('C Export (mgC m-2 day-1)')
    plt.title('C Export in time')
    plt.savefig(inclass.io.outdir + '/SPLAS_CARBON/' + inclass.io.name + '_' + 'TOT_ECO_CARBON_EXPORT_time.png', dpi = 250)
    #plt.show()
    plt.clf()

################################################################################

    plt.close('all')

####################################################################################################################################

####################################################################################################################################
####################################################################################################################################
