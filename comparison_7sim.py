#import matplotlib
#matplotlib.use('Agg')
import numpy as np
from matplotlib import colors
import matplotlib.pyplot as plt
from matplotlib import cm
#import pylab as p
import mpl_toolkits.mplot3d as p3
from matplotlib.font_manager import FontProperties
from sys import argv
import splas
import pickle
from tabulate import tabulate


dir_pickle = "./SPLAS_PICKLE/"

dbm_layer=[]
####################################
# DBM DEPTH 
def dbm_depth(sim1,sim2,sim3,sim4,sim5,sim6,sim7):
    sim=[sim1,sim2,sim3,sim4,sim5,sim6,sim7]
    pb= np.array([np.mean(sim1.phyto_biomass,axis=0),
                  np.mean(sim2.phyto_biomass,axis=0),
                  np.mean(sim3.phyto_biomass,axis=0),
                  np.mean(sim4.phyto_biomass,axis=0),
                  np.mean(sim5.phyto_biomass,axis=0),
                  np.mean(sim6.phyto_biomass,axis=0),
                  np.mean(sim7.phyto_biomass,axis=0)])
    #plt.plot()
    plt.xlabel('Biomass (uMN)')
    plt.ylabel('Depth (m)')
    plt.axis([0,np.amax(pb)+0.02,200,0])
    for j in range(7):
        plt.plot(pb[j], sim[j].z,label= str(sim[j].io.name))
    plt.legend(bbox_to_anchor=(1,1), fontsize = 'small',loc="upper left")
    plt.tight_layout()
    plt.title( 'P Biomass profile')
    plt.savefig('./mesh/'+'profile/'+'Phyto_'+'biomass_profile'+'.png',dpi=250)
    plt.close()
    
    zb= np.array([np.mean(sim1.zoo_biomass,axis=0),
                  np.mean(sim2.zoo_biomass,axis=0),
                  np.mean(sim3.zoo_biomass,axis=0),
                  np.mean(sim4.zoo_biomass,axis=0),
                  np.mean(sim5.zoo_biomass,axis=0),
                  np.mean(sim6.zoo_biomass,axis=0),
                  np.mean(sim7.zoo_biomass,axis=0)])
    plt.plot()
    plt.axis([0,np.amax(zb)+0.02,200,0])
    plt.xlabel('Biomass (uMN)')
    plt.ylabel('Depth (m)') 
    for j in range(7):
        plt.plot(zb[j], sim[j].z)
    plt.legend(bbox_to_anchor=(1,1), fontsize = 'small',loc="upper left",labels=labe)
    plt.tight_layout()
    plt.title( 'Z Biomass profile')
    plt.savefig('./mesh/'+'profile/'+'Zoo_'+'biomass_profile'+'.png',dpi=250)
    plt.close()
    
    plt.plot()
    plt.ylim([200,0])
    plt.xlabel('Biomass (uMN)')
    plt.ylabel('Depth (m)')
    for j in range(7):
        plt.plot((zb[j]+pb[j]), sim[j].z)
    plt.legend(bbox_to_anchor=(1,1), fontsize = 'small',loc="upper left",labels=labe)
    plt.tight_layout()
    plt.title( 'Total Biomass profile')
    plt.savefig('./mesh/'+'profile/'+'Total_'+'biomass_profile'+'.png',dpi=250)
    plt.close()

def diversity(sim1,sim2,sim3,sim4,sim5,sim6,sim7):
    sim=[sim1,sim2,sim3,sim4,sim5,sim6,sim7]
    pd= np.array([np.mean(sim1.phyto_biodiv,axis=0),
                  np.mean(sim2.phyto_biodiv,axis=0),
                  np.mean(sim3.phyto_biodiv,axis=0),
                  np.mean(sim4.phyto_biodiv,axis=0),
                  np.mean(sim5.phyto_biodiv,axis=0),
                  np.mean(sim6.phyto_biodiv,axis=0),
                  np.mean(sim7.phyto_biodiv,axis=0)])
    plt.plot()
    plt.xlabel('shannon index')
    plt.ylabel('Depth (m)')
    plt.axis([0,1.1,200,0])
    for j in range(7):
        plt.plot(pd[j], sim[j].z)
    plt.legend(bbox_to_anchor=(1,1), fontsize = 'small',loc="upper left",labels=labe)
    plt.tight_layout()
    plt.title( 'P diversity profile')
    plt.savefig('./mesh/'+'Diversity/'+'Phyto_'+'diversity_profile'+'.png',dpi=250)
    plt.close()
    
    zd= np.array([np.mean(sim1.zoo_biodiv,axis=0),
                  np.mean(sim2.zoo_biodiv,axis=0),
                  np.mean(sim3.zoo_biodiv,axis=0),
                  np.mean(sim4.zoo_biodiv,axis=0),
                  np.mean(sim5.zoo_biodiv,axis=0),
                  np.mean(sim6.zoo_biodiv,axis=0),
                  np.mean(sim7.zoo_biodiv,axis=0)])
    plt.plot()
    plt.axis([0,1.1,200,0])
    plt.xlabel('shannon index')
    plt.ylabel('Depth (m)') 
    for j in range(7):
        plt.plot(zd[j], sim[j].z)
    plt.legend(bbox_to_anchor=(1,1), fontsize = 'small',loc="upper left",labels=labe)
    plt.tight_layout()
    plt.title( 'Z diversity profile')
    plt.savefig('./mesh/'+'Diversity/'+'Zoo_'+'diversity_profile'+'.png',dpi=250)
    plt.close()
    
    plt.plot()
    plt.axis([0,1.1,200,0])
    plt.xlabel('shannon index')
    plt.ylabel('Depth (m)')
    for j in range(7):
        plt.plot((zd[j]+pd[j])/2, sim[j].z)
    plt.legend(bbox_to_anchor=(1,1), fontsize = 'small',loc="upper left",labels=labe)
    plt.tight_layout()
    plt.title( 'Total mean diversity')
    plt.savefig('./mesh/'+'Diversity/'+'Total_'+'diversity_profile'+'.png',dpi=250)
    plt.close()
    
   

####################################
#N PROFILE
def N_profile(sim1,sim2,sim3,sim4,sim5,sim6,sim7):
    sim=[sim1,sim2,sim3,sim4,sim5,sim6,sim7]
    N = np.array([np.mean(sim1.Nitrogen,axis=0),
                  np.mean(sim2.Nitrogen,axis=0),
                  np.mean(sim3.Nitrogen,axis=0),
                  np.mean(sim4.Nitrogen,axis=0),
                  np.mean(sim5.Nitrogen,axis=0),
                  np.mean(sim6.Nitrogen,axis=0),
                  np.mean(sim7.Nitrogen,axis=0)])
    plt.plot()
    plt.axis([0,np.amax(N)+0.5,200,0])
    for i in range(0,len(N)):
        plt.plot(N[i], sim[i].z)
    plt.legend(bbox_to_anchor=(1,1), fontsize = 'small',loc="upper left",labels=labe)
    plt.tight_layout()
    plt.title('N profile')
    plt.savefig('./mesh/'+'profile/'+ 'N_profile'+'.png',dpi=250)
    plt.close()
####################################
# Temporal screenshots 
def screenshot(sim1,sim2,sim3,sim4,sim5,sim6,sim7):
        calendar=['0','1yr','2yr','3yr','4yr','5yr','6yr','7yr','8yr','9yr','10yr','11yr']
        fig, axis = plt.subplots(2,3)
        fig.suptitle('Phyto_at_' + calendar[t[0]])
        
        
        images = []
        X= sim1.xp
        Y= sim1.z
        x2= sim2.xp
        y2=sim2.z
        #X,Y = np.meshgrid(x,y)
        #X2,Y2 = np.meshgrid(x2,y2)
        ZP=np.array([np.transpose(sim2.phyto[:,t[1],:]*1000/(np.amax(sim2.phyto[:,t[1],:],0)*1000)),
                     np.transpose(sim3.phyto[:,t[1],:]*1000/(np.amax(sim3.phyto[:,t[1],:],0)*1000)),
                     np.transpose(sim4.phyto[:,t[1],:]*1000/(np.amax(sim4.phyto[:,t[1],:],0)*1000)),
                     np.transpose(sim5.phyto[:,t[1],:]*1000/(np.amax(sim5.phyto[:,t[1],:],0)*1000)),
                     np.transpose(sim6.phyto[:,t[1],:]*1000/(np.amax(sim6.phyto[:,t[1],:],0)*1000)),
                     np.transpose(sim7.phyto[:,t[1],:]*1000/(np.amax(sim7.phyto[:,t[1],:],0)*1000))])
            
        images.append(axis[0,0].contourf(X,Y,ZP[0], cmap=plt.get_cmap('rainbow')))
        axis[0,0].label_outer()
        images.append(axis[0,1].contourf(X,Y,ZP[1], cmap=plt.get_cmap('rainbow')))
        axis[0,1].label_outer()
        images.append(axis[0,2].contourf(X,Y,ZP[2], cmap=plt.get_cmap('rainbow')))
        axis[0,2].label_outer()
        images.append(axis[1,0].contourf(X,Y,ZP[3], cmap=plt.get_cmap('rainbow')))
        axis[1,0].label_outer()
        images.append(axis[1,1].contourf(X,Y,ZP[4], cmap=plt.get_cmap('rainbow')))
        axis[1,1].label_outer()
        images.append(axis[1,2].contourf(X,Y,ZP[5], cmap=plt.get_cmap('rainbow')))
        axis[1,2].label_outer()
        
        axis[0,0].set_ylim(200,0)
        axis[0,1].set_ylim(200,0)
        axis[0,2].set_ylim(200,0)
        axis[1,0].set_ylim(200,0)
        axis[1,1].set_ylim(200,0)
        axis[1,2].set_ylim(200,0)
        axis[0,0].set_xlim(0, 100)
        axis[0,1].set_xlim(0, 100)
        axis[0,2].set_xlim(0, 100)
        axis[1,0].set_xlim(0, 100)
        axis[1,1].set_xlim(0, 100)
        axis[1,2].set_xlim(0, 100)
        axis[0,0].set_title(labe[1], fontsize='x-small')
        axis[0,1].set_title(labe[2], fontsize='x-small')
        axis[0,2].set_title(labe[3], fontsize='x-small')
        axis[1,0].set_title(labe[4], fontsize='x-small')
        axis[1,1].set_title(labe[5], fontsize='x-small')
        axis[1,2].set_title(labe[6], fontsize='x-small')
        
                  
        #axis[0,1].label_outer()
    
        #axis[1,0].label_outer()
        vmin = min(0 for image in images)
        vmax = max(1 for image in images)
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
        #for im in images:
        #    im.set_norm(norm)

        fig.colorbar(images[0], ax=axis, orientation='vertical', label='Biomass (uMN)')
        #def update(changed_image):
        #    for im in images:
        #        if (changed_image.get_cmap() != im.get_cmap()
        #                or changed_image.get_clim() != im.get_clim()):
        #            im.set_cmap(changed_image.get_cmap())
        #            im.set_clim(changed_image.get_clim())
        #
        #
        #for im in images:
        #    im.callbacksSM.connect('changed', update)
        #    #cfbar = plt.colorbar()
        #    #plt.xlabel('Size')
        #    #plt.ylabel('Depth')
        #    #cfbar.set_label('Biomass')
        ##fig.tight_layout()
        ##fig.add_subplot(111, frameon=False)
        ##plt.tick_params(labelbottom='off',labelleft='off')
        #
        ##plt.xlabel('Size (μm)')
        ##plt.ylabel('Depth (m)') 
        #
        #
        plt.setp(axis, xticks=np.arange(0,100,20),yticks=np.arange(0,200,25))
        plt.savefig('./mesh/'+'Time/'+'phyto_t_' + calendar[t[0]]+'.png',dpi=250)
        plt.close()

def time_N(sim1,sim2,sim3,sim4,sim5,sim6,sim7):
        calendar=['0','1yr','2yr','3yr','4yr','5yr','6yr','7yr','8yr','9yr','10yr','11yr']
        fig, axis = plt.subplots(2,3)
        fig.suptitle('[N] in time')
        
        
        images = []
        X= sim1.days
        Y= sim1.z
        x2= sim2.xp
        y2=sim2.z
        #X,Y = np.meshgrid(x,y)
        #X2,Y2 = np.meshgrid(x2,y2)
        ZP=np.array([np.transpose(sim2.Nitrogen[:,:]),
                     np.transpose(sim3.Nitrogen[:,:]),
                     np.transpose(sim4.Nitrogen[:,:]),
                     np.transpose(sim5.Nitrogen[:,:]),
                     np.transpose(sim6.Nitrogen[:,:]),
                     np.transpose(sim7.Nitrogen[:,:])])
        
        images.append(axis[0,0].contourf(X,Y,ZP[0], cmap=plt.get_cmap('rainbow')))
        axis[0,0].label_outer()
        images.append(axis[0,1].contourf(X,Y,ZP[1], cmap=plt.get_cmap('rainbow')))
        axis[0,1].label_outer()
        images.append(axis[0,2].contourf(X,Y,ZP[2], cmap=plt.get_cmap('rainbow')))
        axis[0,2].label_outer()
        images.append(axis[1,0].contourf(X,Y,ZP[3], cmap=plt.get_cmap('rainbow')))
        axis[1,0].label_outer()
        images.append(axis[1,1].contourf(X,Y,ZP[4], cmap=plt.get_cmap('rainbow')))
        axis[1,1].label_outer()
        images.append(axis[1,2].contourf(X,Y,ZP[5], cmap=plt.get_cmap('rainbow')))
        axis[1,2].label_outer()
        
        axis[0,0].set_ylim(200,0)
        axis[0,1].set_ylim(200,0)
        axis[0,2].set_ylim(200,0)
        axis[1,0].set_ylim(200,0)
        axis[1,1].set_ylim(200,0)
        axis[1,2].set_ylim(200,0)
        # axis[0,0].set_xlim(0, 100)
        # axis[0,1].set_xlim(0, 100)
        # axis[0,2].set_xlim(0, 100)
        # axis[1,0].set_xlim(0, 100)
        # axis[1,1].set_xlim(0, 100)
        # axis[1,2].set_xlim(0, 100)
        axis[0,0].set_title(labe[1], fontsize='x-small')
        axis[0,1].set_title(labe[2], fontsize='x-small')
        axis[0,2].set_title(labe[3], fontsize='x-small')
        axis[1,0].set_title(labe[4], fontsize='x-small')
        axis[1,1].set_title(labe[5], fontsize='x-small')
        axis[1,2].set_title(labe[6], fontsize='x-small')
        
                  
        #axis[0,1].label_outer()
    
        #axis[1,0].label_outer()
        vmin = min(image.get_array().min() for image in images)
        vmax = max(image.get_array().max() for image in images)
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
        for im in images:
            im.set_norm(norm)

        fig.colorbar(images[0], ax=axis, orientation='vertical', fraction=.1, label='[N] (uMN)')
        def update(changed_image):
            for im in images:
                if (changed_image.get_cmap() != im.get_cmap()
                        or changed_image.get_clim() != im.get_clim()):
                    im.set_cmap(changed_image.get_cmap())
                    im.set_clim(changed_image.get_clim())


        for im in images:
            im.callbacksSM.connect('changed', update)
            #cfbar = plt.colorbar()
            #plt.xlabel('Size')
            #plt.ylabel('Depth')
            #cfbar.set_label('Biomass')
        #fig.tight_layout()
        #fig.add_subplot(111, frameon=False)
        #plt.tick_params(labelbottom='off',labelleft='off')
        
        #plt.xlabel('Size (μm)')
        #plt.ylabel('Depth (m)') 
        
        
        #plt.setp(axis, xticks=np.arange(0,100,20),yticks=np.arange(0,200,25))
        plt.savefig('./mesh/'+'Time/'+'n_t_' + calendar[t[0]]+'.png',dpi=250)
        plt.close()


# Temporal screenshots 
def screenshot_diff(sim1,sim2,sim3,sim4,sim5,sim6,sim7):
        calendar=['0','1yr','2yr','3yr','4yr','5yr','6yr','7yr','8yr','9yr','10yr','11yr']
        sim=[sim1,sim2,sim3,sim4,sim5,sim6,sim7]
        fig, axis = plt.subplots(2,3)
               
        fig.suptitle('Phyto_variation_at_' + calendar[t[0]])
        
        images = []
        X= sim1.xp
        Y= sim1.z
        x2= sim2.xp
        y2=sim2.z
        #X,Y = np.meshgrid(x,y)
        #X2,Y2 = np.meshgrid(x2,y2)
        ZP=np.array([np.transpose(sim1.phyto[:,t[1],:]*1000/(np.amax(sim1.phyto[:,t[1],:],0)*1000)),
                     np.transpose(sim2.phyto[:,t[1],:]*1000/(np.amax(sim2.phyto[:,t[1],:],0)*1000)),
                     np.transpose(sim3.phyto[:,t[1],:]*1000/(np.amax(sim3.phyto[:,t[1],:],0)*1000)),
                     np.transpose(sim4.phyto[:,t[1],:]*1000/(np.amax(sim4.phyto[:,t[1],:],0)*1000)),
                     np.transpose(sim5.phyto[:,t[1],:]*1000/(np.amax(sim5.phyto[:,t[1],:],0)*1000)),
                     np.transpose(sim6.phyto[:,t[1],:]*1000/(np.amax(sim6.phyto[:,t[1],:],0)*1000)),
                     np.transpose(sim7.phyto[:,t[1],:]*1000/(np.amax(sim7.phyto[:,t[1],:],0)*1000))])
    
        images.append(axis[0,0].contourf(X,Y,(ZP[1]-ZP[0]), cmap=plt.get_cmap('coolwarm')))
        axis[0,0].label_outer()
        images.append(axis[0,1].contourf(X,Y,(ZP[2]-ZP[0]), cmap=plt.get_cmap('coolwarm')))
        axis[0,1].label_outer()
        images.append(axis[0,2].contourf(X,Y,(ZP[3]-ZP[0]), cmap=plt.get_cmap('coolwarm')))
        axis[0,2].label_outer()
        images.append(axis[1,0].contourf(X,Y,(ZP[4]-ZP[0]), cmap=plt.get_cmap('coolwarm')))
        axis[1,0].label_outer()
        images.append(axis[1,1].contourf(X,Y,(ZP[5]-ZP[0]), cmap=plt.get_cmap('coolwarm')))
        axis[1,1].label_outer()
        images.append(axis[1,2].contourf(X,Y,(ZP[6]-ZP[0]), cmap=plt.get_cmap('coolwarm')))
        axis[1,2].label_outer()
        
        axis[0,0].set_ylim(200,0)
        axis[0,1].set_ylim(200,0)
        axis[0,2].set_ylim(200,0)
        axis[1,0].set_ylim(200,0)
        axis[1,1].set_ylim(200,0)
        axis[1,2].set_ylim(200,0)
        axis[0,0].set_xlim(0, 100)
        axis[0,1].set_xlim(0, 100)
        axis[0,2].set_xlim(0, 100)
        axis[1,0].set_xlim(0, 100)
        axis[1,1].set_xlim(0, 100)
        axis[1,2].set_xlim(0, 100)
        axis[0,0].set_title(labe[1], fontsize='x-small')
        axis[0,1].set_title(labe[2], fontsize='x-small')
        axis[0,2].set_title(labe[3], fontsize='x-small')
        axis[1,0].set_title(labe[4], fontsize='x-small')
        axis[1,1].set_title(labe[5], fontsize='x-small')
        axis[1,2].set_title(labe[6], fontsize='x-small')
                   
        #axis[0,1].label_outer()
    
        #axis[1,0].label_outer()
        #vmin = min(0 for image in images)
        #vmax = max(1 for image in images)
        #norm = colors.Normalize(vmin=vmin, vmax=vmax)
        #for im in images:
        #    im.set_norm(norm)

        fig.colorbar(images[0], ax=axis, orientation='vertical', fraction=.1, label='Biomass variation(uMN)')
        #def update(changed_image):
        #    for im in images:
        #        if (changed_image.get_cmap() != im.get_cmap()
        #                or changed_image.get_clim() != im.get_clim()):
        #            im.set_cmap(changed_image.get_cmap())
        #            im.set_clim(changed_image.get_clim())
        #
        #
        #for im in images:
        #    im.callbacksSM.connect('changed', update)
        #    #cfbar = plt.colorbar()
        #    #plt.title(str(sim[im]))
        #    #plt.ylabel('Depth')
        #    #cfbar.set_label('Biomass')
        #     
        #    
        ##fig.tight_layout()
        ##fig.add_subplot(111, frameon=False)
        ##plt.tick_params(labelbottom='off',labelleft='off')
        #
        ##plt.xlabel('Size (μm)')
        ##plt.ylabel('Depth (m)')
        tic=int((sim1.tdim-300)/10)        
        
        calendar=["0",'1yr','2yr','3yr','4yr','5yr','6yr','7yr','8yr','9yr','10yr'] 
        plt.setp(axis, xticks=np.arange(0,100,20),yticks=np.arange(0,200,25))    
        plt.savefig('./mesh/'+'Tdiff/'+'phyto_t_' + calendar[t[0]] +'.png',dpi=250)
        plt.close() 

def time_diff_N(sim1,sim2,sim3,sim4,sim5,sim6,sim7):
        calendar=['0','1yr','2yr','3yr','4yr','5yr','6yr','7yr','8yr','9yr','10yr','11yr']
        sim=[sim1,sim2,sim3,sim4,sim5,sim6,sim7]
        fig, axis = plt.subplots(2,3)
               
        fig.suptitle('[N] variation' )
        
        images = []
        X= sim1.days
        Y= sim1.z
        x2= sim2.xp
        y2=sim2.z
        #X,Y = np.meshgrid(x,y)
        #X2,Y2 = np.meshgrid(x2,y2)
        ZP=np.array([np.transpose(sim1.Nitrogen[:,:]),
                     np.transpose(sim2.Nitrogen[:,:]),
                     np.transpose(sim3.Nitrogen[:,:]),
                     np.transpose(sim4.Nitrogen[:,:]),
                     np.transpose(sim5.Nitrogen[:,:]),
                     np.transpose(sim6.Nitrogen[:,:]),
                     np.transpose(sim7.Nitrogen[:,:])])
    
        images.append(axis[0,0].contourf(X,Y,(ZP[1]-ZP[0]), cmap=plt.get_cmap('coolwarm')))
        axis[0,0].label_outer()
        images.append(axis[0,1].contourf(X,Y,(ZP[2]-ZP[0]), cmap=plt.get_cmap('coolwarm')))
        axis[0,1].label_outer()
        images.append(axis[0,2].contourf(X,Y,(ZP[3]-ZP[0]), cmap=plt.get_cmap('coolwarm')))
        axis[0,2].label_outer()
        images.append(axis[1,0].contourf(X,Y,(ZP[4]-ZP[0]), cmap=plt.get_cmap('coolwarm')))
        axis[1,0].label_outer()
        images.append(axis[1,1].contourf(X,Y,(ZP[5]-ZP[0]), cmap=plt.get_cmap('coolwarm')))
        axis[1,1].label_outer()
        images.append(axis[1,2].contourf(X,Y,(ZP[6]-ZP[0]), cmap=plt.get_cmap('coolwarm')))
        axis[1,2].label_outer()
        
        axis[0,0].set_ylim(200,0)
        axis[0,1].set_ylim(200,0)
        axis[0,2].set_ylim(200,0)
        axis[1,0].set_ylim(200,0)
        axis[1,1].set_ylim(200,0)
        axis[1,2].set_ylim(200,0)
        # axis[0,0].set_xlim(0, 100)
        # axis[0,1].set_xlim(0, 100)
        # axis[0,2].set_xlim(0, 100)
        # axis[1,0].set_xlim(0, 100)
        # axis[1,1].set_xlim(0, 100)
        # axis[1,2].set_xlim(0, 100)
        axis[0,0].set_title(labe[1], fontsize='x-small')
        axis[0,1].set_title(labe[2], fontsize='x-small')
        axis[0,2].set_title(labe[3], fontsize='x-small')
        axis[1,0].set_title(labe[4], fontsize='x-small')
        axis[1,1].set_title(labe[5], fontsize='x-small')
        axis[1,2].set_title(labe[6], fontsize='x-small')
                   
        #axis[0,1].label_outer()
    
        #axis[1,0].label_outer()
        vmin = min(image.get_array().min() for image in images)
        vmax = max(image.get_array().max() for image in images)
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
        for im in images:
            im.set_norm(norm)

        fig.colorbar(images[0], ax=axis, orientation='vertical', fraction=.1, label='[N] variation(uMN)')
        def update(changed_image):
            for im in images:
                if (changed_image.get_cmap() != im.get_cmap()
                        or changed_image.get_clim() != im.get_clim()):
                    im.set_cmap(changed_image.get_cmap())
                    im.set_clim(changed_image.get_clim())


        for im in images:
            im.callbacksSM.connect('changed', update)
            #cfbar = plt.colorbar()
            #plt.title(str(sim[im]))
            #plt.ylabel('Depth')
            #cfbar.set_label('Biomass')
             
            
        #fig.tight_layout()
        #fig.add_subplot(111, frameon=False)
        #plt.tick_params(labelbottom='off',labelleft='off')
        
        #plt.xlabel('Size (μm)')
        #plt.ylabel('Depth (m)')
        tic=int((sim1.tdim-300)/10)        
        
        calendar=["0",'1yr','2yr','3yr','4yr','5yr','6yr','7yr','8yr','9yr','10yr'] 
        #plt.setp(axis, xticks=np.arange(0,100,20),yticks=np.arange(0,200,25))    
        plt.savefig('./mesh/'+'Tdiff/'+'n_t_' + '.png',dpi=250)
        plt.close()         

def diversity_diff(sim1,sim2,sim3,sim4,sim5,sim6,sim7):
        sim=[sim1,sim2,sim3,sim4,sim5,sim6,sim7]
        x=range(sim1.tdim-300)
        tic=int((sim1.tdim-300)/10)
        fig, axis = plt.subplots(2,3)
              
        fig.suptitle('Diversity_in_time_vs_default')
        
        images = []
        X= sim1.days
        Y= sim1.z
        
        #X,Y = np.meshgrid(x,y)
        #X2,Y2 = np.meshgrid(x2,y2)
        ZP=np.array([np.transpose(sim1.phyto_biodiv),
                    np.transpose(sim2.phyto_biodiv),
                    np.transpose(sim3.phyto_biodiv),
                    np.transpose(sim4.phyto_biodiv),
                    np.transpose(sim5.phyto_biodiv),
                    np.transpose(sim6.phyto_biodiv),
                    np.transpose(sim7.phyto_biodiv)])
    
        images.append(axis[0,0].contourf(X,Y,(ZP[1]-ZP[0]), cmap=plt.get_cmap('coolwarm')))
        axis[0,0].label_outer()
        images.append(axis[0,1].contourf(X,Y,(ZP[2]-ZP[0]), cmap=plt.get_cmap('coolwarm')))
        axis[0,1].label_outer()
        images.append(axis[0,2].contourf(X,Y,(ZP[3]-ZP[0]), cmap=plt.get_cmap('coolwarm')))
        axis[0,2].label_outer()
        images.append(axis[1,0].contourf(X,Y,(ZP[4]-ZP[0]), cmap=plt.get_cmap('coolwarm')))
        axis[1,0].label_outer()
        images.append(axis[1,1].contourf(X,Y,(ZP[5]-ZP[0]), cmap=plt.get_cmap('coolwarm')))
        axis[1,1].label_outer()
        images.append(axis[1,2].contourf(X,Y,(ZP[6]-ZP[0]), cmap=plt.get_cmap('coolwarm')))
        axis[1,2].label_outer()
        
        axis[0,0].set_title(labe[1], fontsize='x-small')
        axis[0,1].set_title(labe[2], fontsize='x-small')
        axis[0,2].set_title(labe[3], fontsize='x-small')
        axis[1,0].set_title(labe[4], fontsize='x-small')
        axis[1,1].set_title(labe[5], fontsize='x-small')
        axis[1,2].set_title(labe[6], fontsize='x-small')
        axis[0,0].set_ylim(200,0)
        axis[0,1].set_ylim(200,0)
        axis[0,2].set_ylim(200,0)
        axis[1,0].set_ylim(200,0)
        axis[1,1].set_ylim(200,0)
        axis[1,2].set_ylim(200,0)
        
        
                  
        #axis[0,1].label_outer()
    
        #axis[1,0].label_outer()
        vmin = min(image.get_array().min() for image in images)
        vmax = max(image.get_array().max() for image in images)
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
        for im in images:
            im.set_norm(norm)

        fig.colorbar(images[0], ax=axis, orientation='vertical', fraction=.1, label='Shannon evenness variation ')
        def update(changed_image):
            for im in images:
                if (changed_image.get_cmap() != im.get_cmap()
                        or changed_image.get_clim() != im.get_clim()):
                    im.set_cmap(changed_image.get_cmap())
                    im.set_clim(changed_image.get_clim())


        for im in images:
            im.callbacksSM.connect('changed', update)
            #plt.setp(axis, xticks=([],yticks=np.arange(0,200,25))
        #plt.setp(axis,xticks=sim1.days,xticklabels= "")   
        #plt.setp(axis, yticks=np.arange(0,200,25))
        #plt.setp(axis, xticks=np.arange(0,3650,365),label= [0,'1yr','2yr','3yr','4yr','5yr','6yr','7yr','8yr','9yr','10yr'])
        #fig.tight_layout()
        #fig.add_subplot(111, frameon=False)
        #plt.tick_params(labelbottom='off',labelleft='off')
        #fig.tight_layout()
        #plt.xlabel('Time')
        #plt.ylabel('Depth (m)')  
        
        
        #xticks(np.arange(10), )
        #fig.supxlabel('Size')
        #fig.supylabel('Depth')            
        plt.savefig('./mesh/'+'Diversity/'+'phyto_diversity_vs_default'+'.png',dpi=250)
        plt.close() 

def light_diff(sim1,sim2,sim3,sim4,sim5,sim6,sim7):
        sim=[sim1,sim2,sim3,sim4,sim5,sim6,sim7]
        
        for k in range(sim1.io.m -1):
            I=[]
            for j in range(7):
                
                #irr=sim[j].irrprofile[10:,k]
                irr=sim[j].irrprofile[10:100,k]
                I.append(irr[10:])
                #print(np.amax(I))
                
                #plt.plot(np.arange(0,sim[j].tdim-10,1), irr, label= str(sim[j].io.name))
                plt.plot(np.arange(10,100,1), irr, label= str(sim[j].io.name))
                plt.legend(loc='lower center', fontsize='x-small')
                plt.tight_layout()
                #plt.axis([0,11, 0.3*np.amin(I) , 1.1*np.amax(I)])
                #plt.xticks(np.arange(0,3650,365),labels= [0,'1yr','2yr','3yr','4yr','5yr','6yr','7yr','8yr','9yr','10yr'])
            plt.xlabel('Time')
            plt.ylabel('Irradiance (W/m^2)')
            plt.legend(bbox_to_anchor=(1,1), fontsize = 'small',loc="upper left",labels=labe)
            plt.tight_layout()
            
            plt.savefig('./mesh/'+'light/'+'Layer_' + str(k) + '.png')
            plt.clf()
        plt.close()
                
######PRINT n classes > mean

def nspec(sim1,sim2,sim3,sim4,sim5,sim6,sim7):
    sim=[sim1,sim2,sim3,sim4,sim5,sim6,sim7]       
    
    
    plt.plot()
    #plt.set_prop_cycle(colour_blind_cycle)
    
    
    pt=[]
    
    for j in range(7):
        p=[]
        for t in range(11):
            temp=sim[j].phyto_biomass[int((t)*((sim1.tdim/10)-1)),:]
            max=np.amax(temp)
            dbm=np.where(temp == max)[0]
            
            print(dbm)
            zp=np.array([sim[j].phyto[:,int((t)*((sim1.tdim/10)-1)),dbm]])
            
            mp=np.sum(zp>np.mean(zp))+0.01*j
            if t == 0:
                mp=0
            else:
                pass
            if t == 10:
                dbm_layer.append(dbm)
                #print(dbm_layer)
            else: 
                pass
            p.append(mp)
            
            
            #print(p)
        #print(dbm_layer)
        #plt.scatter(range(10),p)
        plt.plot(range(11),p, marker='o', label=[str(sim[j].io.name), 'Layer_',str(dbm)])
        print('ok,plotted!')
        print(str(dbm_layer)+'dbm_layer')
        pt.append(p)
    print(pt[1])
  
    plt.xlabel('Time')
    plt.xticks(np.arange(11),labels= [0,'1yr','2yr','3yr','4yr','5yr','6yr','7yr','8yr','9yr','10yr'], rotation='vertical')
    plt.ylabel('n. of classes> mean biomass')
    plt.xlim([-1,11])    
     
    plt.legend(bbox_to_anchor=(1,1), fontsize = 'small',loc="upper left")
    plt.tight_layout()
    plt.title( 'P nspec')
    plt.savefig('./mesh/'+'profile/'+'Phyto_'+'nspec_time''.png',dpi=250)
    plt.close()
####################################
###### zoom on dbm
def dbm_zoom(sim1,sim2,sim3,sim4,sim5,sim6,sim7):
    sim=[sim1,sim2,sim3,sim4,sim5,sim6,sim7]       
    pos = [[0,0],[0,1],[1,0],[1,1]]
    if dbm_layer.size == 1:
        for i in range(sim1.xp.size):   
            for j in range(7):
                
                
                temp=sim[j].phyto[i,100:,dbm_layer[0]]
                X = sim1.days[100:]
                
                
                plt.plot(X, temp)
            plt.ylabel('Layer'+str(dbm_layer[0]))
            plt.xlabel('Time(days)')
            plt.title('size-class' + str(i) + 'biomass')
            plt.legend(bbox_to_anchor=(1,1), labels=labe, fontsize = 'small',loc="upper left")
            plt.tight_layout()
            plt.savefig('./mesh/'+'dbm_zoom/'+'P_class_'+str(i)+'.png',dpi=250)
            plt.clf()
    else:
        for i in range(sim1.xp.size):
            fig, Axis = plt.subplots(dbm_layer.size,1)
            for l in range(dbm_layer.size):
                
                for j in range(7):
                
                
                    temp=sim[j].phyto[i,100:,dbm_layer[l]]
                    X = sim1.days[100:]
                    
                    
                    Axis[l].plot(X, temp)
                    Axis[l].set_ylabel('Layer'+str(dbm_layer[l]))
                    
            
            Axis[dbm_layer.size-1].set_xlabel('Time(days)')
             
               
            #fig.supxlabel('Time')
            #plt.xticks(np.arange(11),labels= [0,'1yr','2yr','3yr','4yr','5yr','6yr','7yr','8yr','9yr','10yr'], rotation='vertical')
            #fig.supylabel('Biomass ')
                
            fig.suptitle('size-class' + str(i) + 'biomass')
            fig.legend(bbox_to_anchor=(1,1), labels=labe, fontsize = 'small',loc="upper left")
            fig.tight_layout()
            plt.savefig('./mesh/'+'dbm_zoom/'+'P_class_'+str(i)+'.png',dpi=250)
            plt.clf()
        plt.close()

###### PRINT biomass graphs

def biomass(sim1,sim2,sim3,sim4,sim5,sim6,sim7):
    sim=[sim1,sim2,sim3,sim4,sim5,sim6,sim7]
    x=range(sim1.tdim-300)
    tic=int((sim1.tdim-300)/10)
     
    for k in range(7):
        #print(sim[k].phyto_biomass.sum(1))
        #print(sim[k].phyto_biomass.sum(0))
        y=np.array([sim[k].phyto_biomass[300:,:].sum(1),sim[k].zoo_biomass[300:,:].sum(1)])
        y2=np.sum(y,axis=0)
        
        #print(y2.shape)
        
        plt.plot(x, y2)
    plt.xlabel('Time')
    plt.ylabel('Total system biomass (uMN)')
    plt.xticks(np.arange(0,sim1.tdim-300,tic),labels= [0,'1yr','2yr','3yr','4yr','5yr','6yr','7yr','8yr','9yr','10yr'], rotation='vertical')    
    plt.legend(bbox_to_anchor=(1,1), fontsize = 'small',loc="upper left",labels=labe)
    plt.tight_layout()    
    plt.savefig('./mesh/'+'biomass/'+'System_biomass'+'.png')
    plt.close()

#####PRINT carbon graphs    

def carbon_out(sim1,sim2,sim3,sim4,sim5,sim6,sim7):
    sim=[sim1,sim2,sim3,sim4,sim5,sim6,sim7]
    print(sim1.tdim)
    #tic=int((sim1.tdim-300)/10)
    for j in range(7):
        plt.plot(sim1.days[1:],sim[j].total_carbon_export[1:], label =str(sim[j].io.name).format(j=j))
    
    plt.title('Carbon export vs time')
    
    plt.legend(bbox_to_anchor=(1,1), fontsize = 'small',loc="upper left",labels=labe)
    plt.xlabel('Time')
    plt.ylabel('Carbon export (mgC m-2 day-1)')
    plt.xlim([0,sim1.io.simulation_length])
    plt.ylim([0,np.amax(sim3.total_carbon_export[1:],0)+1])
    plt.tight_layout()
    plt.savefig('./mesh/'+'carbon/'+ 'Tot_c_export'+'.png', dpi=250)
    plt.clf()

def carbon_light(sim1,sim2,sim3,sim4,sim5,sim6,sim7):
    sim=[sim1,sim2,sim3,sim4,sim5,sim6,sim7]
    clm_tot=[]
    clf_tot=[]
    irr=[]
    
    for j in range(7):
        clm=[]
        clf=[]
        
        clf= sim[j].total_carbon_export[-1]
        clm= sim[j].total_ave_carbon_export
        irr.append(sim[j].io.i)
        clm_tot.append([clm,sim[j].io.i])
        clf_tot.append([clf,sim[j].io.i])
        
        plt.plot(irr[j],clf_tot[j][0], marker='o',linestyle='None')
    #clm_tot=np.sort(clm_tot,axis=0)
    #clf_tot=np.sort(clf_tot,axis=0)    
    plt.legend(bbox_to_anchor=(1,1), fontsize = 'small',loc="upper left",labels=labe)    
    #x_val=irr #[x[0] for x in irr]
    #y_val=np.sort([x[0] for x in clf_tot], axis=0)
    #y_val2=np.sort([x[0] for x in clm_tot],axis=0)    
    #plt.plot(x_val,y_val, 'k')
    
   
    plt.xlabel('Irradiance')
    plt.ylabel('Carbon export')
    plt.tight_layout()
    plt.savefig('./mesh/'+'carbon/'+ 'c_export_vs_I_final'+'.png', dpi=250)
    plt.clf()
    for j in range(7):
        plt.scatter(irr[j],clm_tot[j][0], marker='s', linestyle='None')
    plt.legend(bbox_to_anchor=(1,1), fontsize = 'small',loc="upper left",labels=labe)    
    #x_val=irr #[x[0] for x in irr]
    #y_val=clf_tot #[x[0] for x in clf]
    #y_val2=clm_tot #[x[0] for x in clm]   
    #plt.plot(x_val,y_val2, 'k') 
    plt.xlabel('Irradiance')
    plt.ylabel('Mean Carbon export')
    plt.tight_layout()
    plt.savefig('./mesh/'+'carbon/'+ 'mean_export_vs_I'+'.png', dpi=250)
    plt.clf()    
    plt.close()

#####PRINT sensitivity table    

def sens_table(sim1,sim2,sim3,sim4,sim5,sim6,sim7):
    sim=[sim1,sim2,sim3,sim4,sim5,sim6,sim7]
    f = open("./mesh/"+ "table/" +"sensitivity_table.txt", "a")
    xlab= [" ","I + 10%", "I + 20%", "I + 30%", "I - 10%", "I - 20%", "I - 30%"]
    head= [" ", "C_Export", "P_mass", "Z_mass", "[N]", "[OM]", "PSE", "ZSE"]
    table = []
    for j in range(1,7):
        C_Export = np.mean(100*(np.array(sim[j].total_ave_carbon_export) - np.array(sim[0].total_ave_carbon_export))/np.array(sim[0].total_ave_carbon_export))
        print(C_Export)
        P_mass = np.mean(100*(np.array(sim[j].avez_p_biomass) - np.array(sim[0].avez_p_biomass))/np.array(sim[0].avez_p_biomass))
        Z_mass = np.mean(100*(np.array(sim[j].avez_z_biomass) - np.array(sim[0].avez_z_biomass))/np.array(sim[0].avez_z_biomass))
        N = np.mean(100*(np.array(sim[j].avez_N) - np.array(sim[0].avez_N))/np.array(sim[0].avez_N))
        Detritus = np.mean(100*(np.array(sim[j].mean_detritus) - np.array(sim[0].mean_detritus))/np.array(sim[0].mean_detritus))
        PSE = np.mean(100*(np.array(sim[j].avez_shannon_p) - np.array(sim[0].avez_shannon_p))/np.array(sim[0].avez_shannon_p))
        ZSE = np.mean(100*((np.array(sim[j].avez_shannon_z) - np.array(sim[0].avez_shannon_z))/np.array(sim[0].avez_shannon_z)))
        table.append([xlab[j],C_Export,P_mass,Z_mass,N,Detritus,PSE,ZSE])
    tavola= tabulate(table,headers=head, tablefmt="plain")
    print( tavola, file= f)
    f.close()
    
def winners(sim1,sim2,sim3,sim4,sim5,sim6,sim7):
    bb=['name','biomass','p_class']
    sim=[sim1,sim2,sim3,sim4,sim5,sim6,sim7]
    f = open("./mesh/" + "winners/"+ "classifica.txt",'a')
    for i in range(sim1.io.m):
        valori=np.array([])
        classifica=[]
        classi=np.array([])
        # classifica[0].append('name')
        # classifica[1].append('biomass')
        # classifica[2].append('p_class')
        for j in range(7):
            temp=sim[j].fphyto[:,i]
            winner=np.amax(temp)
            winner_class=np.where(temp==winner)
            
            print(str(winner) + 'Winner')
            print(str(winner_class) + 'class')
            classifica.append([str(labe[j]),winner,winner_class[0]])
            #classifica[1].append(winner)
            #classifica[2].append(winner_class[0])
        print(classifica)
        df=sorted(classifica,key=lambda x:x[1], reverse=True)
        #print(df,file=f)
        #cl=sorted(classifica,key = lambda classifica:classifica[1])
        # np.insert(classifica[0],values=bb[0])
        # np.insert(classifica[1],values=bb[1])
        # np.insert(classifica[2],values=bb[2])
        #print(classifica)        
        print( 'Layer' + str(i), file=f) 
        print('NAME_SIM      BIOMASS     P_CLASS',file=f)
        print('-------------------------------------------',file=f)
        print(tabulate(df, tablefmt='plain'), sep="\n", file=f)
        print('-------------------------------------------',file=f)
    f.close()        
                

####################################
#loading files and avez calculating
def size_spectrum(sim1,sim2,sim3,sim4,sim5,sim6,sim7):
        calendar=['0','1yr','2yr','3yr','4yr','5yr','6yr','7yr','8yr','9yr','10yr','11yr']
        
        for k in range(sim1.io.m):
            fig, axis = plt.subplots(2,3)
            fig.suptitle('Phyto_spectrum x= time y=size layer_' + str(k))
            
            
            images = []
            X= sim1.days
            Y= sim1.xp
            x2= sim2.xp
            y2=sim2.z
            #X,Y = np.meshgrid(x,y)
            #X2,Y2 = np.meshgrid(x2,y2)
            ZP=np.array([(sim2.phyto[:,:,k]*1000/(np.amax(sim2.phyto[:,:,k],0)*1000)),
                         (sim3.phyto[:,:,k]*1000/(np.amax(sim3.phyto[:,:,k],0)*1000)),
                         (sim4.phyto[:,:,k]*1000/(np.amax(sim4.phyto[:,:,k],0)*1000)),
                         (sim5.phyto[:,:,k]*1000/(np.amax(sim5.phyto[:,:,k],0)*1000)),
                         (sim6.phyto[:,:,k]*1000/(np.amax(sim6.phyto[:,:,k],0)*1000)),
                         (sim7.phyto[:,:,k]*1000/(np.amax(sim7.phyto[:,:,k],0)*1000))])
                
            images.append(axis[0,0].contourf(X,Y,ZP[0],vmin=0.0,vmax=1.0, cmap=plt.get_cmap('rainbow')))
            axis[0,0].label_outer()
            images.append(axis[0,1].contourf(X,Y,ZP[1],vmin=0.0,vmax=1.0, cmap=plt.get_cmap('rainbow')))
            axis[0,1].label_outer()
            images.append(axis[0,2].contourf(X,Y,ZP[2],vmin=0.0,vmax=1.0, cmap=plt.get_cmap('rainbow')))
            axis[0,2].label_outer()
            images.append(axis[1,0].contourf(X,Y,ZP[3],vmin=0.0,vmax=1.0, cmap=plt.get_cmap('rainbow')))
            axis[1,0].label_outer()
            images.append(axis[1,1].contourf(X,Y,ZP[4],vmin=0.0,vmax=1.0, cmap=plt.get_cmap('rainbow')))
            axis[1,1].label_outer()
            images.append(axis[1,2].contourf(X,Y,ZP[5],vmin=0.0,vmax=1.0, cmap=plt.get_cmap('rainbow')))
            axis[1,2].label_outer()
            
            axis[0,0].set_ylim(0,100)
            axis[0,1].set_ylim(0,100)
            axis[0,2].set_ylim(0,100)
            axis[1,0].set_ylim(0,100)
            axis[1,1].set_ylim(0,100)
            axis[1,2].set_ylim(0,100)
            axis[0,0].set_title(labe[1], fontsize='x-small')
            axis[0,1].set_title(labe[2], fontsize='x-small')
            axis[0,2].set_title(labe[3], fontsize='x-small')
            axis[1,0].set_title(labe[4], fontsize='x-small')
            axis[1,1].set_title(labe[5], fontsize='x-small')
            axis[1,2].set_title(labe[6], fontsize='x-small')
            
                      
            #axis[0,1].label_outer()
            
            #axis[1,0].label_outer()
            
            fig.colorbar(images[0], ax=axis, orientation='vertical', fraction=.1, label='Biomass (uMN)')
           
            plt.savefig('./mesh/'+'size_spectrum/'+'phyto_t_layer' + str(k)+'.png',dpi=250)
            plt.close()

def final_size(sim1,sim2,sim3,sim4,sim5,sim6,sim7):
        calendar=['0','1yr','2yr','3yr','4yr','5yr','6yr','7yr','8yr','9yr','10yr','11yr']
        
    
        fig, axis = plt.subplots(2,3)
        fig.suptitle('Phyto_spectrum_final')
        
        
        images = []
        X= sim1.xp
        Y= sim1.z
        
        #X,Y = np.meshgrid(x,y)
        #X2,Y2 = np.meshgrid(x2,y2)
        ZP=np.array([np.transpose(sim2.fphyto*1000/(np.amax(sim2.fphyto,0)*1000)),
                     np.transpose(sim3.fphyto*1000/(np.amax(sim3.fphyto,0)*1000)),
                     np.transpose(sim4.fphyto*1000/(np.amax(sim4.fphyto,0)*1000)),
                     np.transpose(sim5.fphyto*1000/(np.amax(sim5.fphyto,0)*1000)),
                     np.transpose(sim6.fphyto*1000/(np.amax(sim6.fphyto,0)*1000)),
                     np.transpose(sim7.fphyto*1000/(np.amax(sim7.fphyto,0)*1000))])
            
        images.append(axis[0,0].contourf(X,Y,ZP[0],vmin=0.0,vmax=1.0, cmap=plt.get_cmap('rainbow')))
        axis[0,0].label_outer()
        images.append(axis[0,1].contourf(X,Y,ZP[1],vmin=0.0,vmax=1.0, cmap=plt.get_cmap('rainbow')))
        axis[0,1].label_outer()
        images.append(axis[0,2].contourf(X,Y,ZP[2],vmin=0.0,vmax=1.0, cmap=plt.get_cmap('rainbow')))
        axis[0,2].label_outer()
        images.append(axis[1,0].contourf(X,Y,ZP[3],vmin=0.0,vmax=1.0, cmap=plt.get_cmap('rainbow')))
        axis[1,0].label_outer()
        images.append(axis[1,1].contourf(X,Y,ZP[4],vmin=0.0,vmax=1.0, cmap=plt.get_cmap('rainbow')))
        axis[1,1].label_outer()
        images.append(axis[1,2].contourf(X,Y,ZP[5],vmin=0.0,vmax=1.0, cmap=plt.get_cmap('rainbow')))
        axis[1,2].label_outer()
        
        axis[0,0].set_ylim(200,0)
        axis[0,1].set_ylim(200,0)
        axis[0,2].set_ylim(200,0)
        axis[1,0].set_ylim(200,0)
        axis[1,1].set_ylim(200,0)
        axis[1,2].set_ylim(200,0)
        axis[0,0].set_title(labe[1], fontsize='x-small')
        axis[0,1].set_title(labe[2], fontsize='x-small')
        axis[0,2].set_title(labe[3], fontsize='x-small')
        axis[1,0].set_title(labe[4], fontsize='x-small')
        axis[1,1].set_title(labe[5], fontsize='x-small')
        axis[1,2].set_title(labe[6], fontsize='x-small')
        
                  
        
        
        fig.colorbar(images[0], ax=axis, orientation='vertical', fraction=.1, label='Relative fraction of Biomass')
       
        
        
       
        plt.savefig('./mesh/'+'biomass/'+'phyto_final_spectrum' + '.png',dpi=250)
        plt.close()




def load(self,direct,name):
    self.h = pickle.load(open(direct + name + '_' + "total_depth.p","rb"))
    self.phyto_biomass = pickle.load(open(direct + name + '_' + "phyto_biomass.p","rb"))
    self.zoo_biomass = pickle.load(open(direct + name + '_' + "zoo_biomass.p","rb"))
    self.xp = pickle.load(open(direct + name + '_' + "xp_phyto_ESD.p","rb"))
    self.xz = pickle.load(open(direct + name + '_' + "xz_zoo_ESD.p","rb"))
    self.z = pickle.load(open(direct + name + '_' + "z_points.p","rb"))
    self.phyto = pickle.load(open(direct + name + '_' + "phytoplankton.p","rb"))
    self.zoo = pickle.load(open(direct + name + '_' + "zooplankton.p","rb"))
    self.tdim = pickle.load(open(direct + name + '_' + "time_dimension.p","rb"))
    self.phyto_biodiv = pickle.load(open(direct + name + '_' + "phyto_biodiv.p","rb"))
    self.zoo_biodiv = pickle.load(open(direct + name + '_' + "zoo_biodiv.p","rb"))
    self.Nitrogen = pickle.load(open(direct + name + '_' + "nitrogen.p","rb"))
    self.fphyto = pickle.load(open(direct + name + '_' + "phyto_final_spectrum.p","rb"))
    self.phyto_mean = pickle.load(open(direct + name + '_' + "phyto_mean_spectrum.p","rb"))
    self.zoo_mean = pickle.load(open(direct + name + '_' + "zoo_mean_spectrum.p","rb"))
    self.total_carbon_export = pickle.load(open(direct + name + '_' + "total_carbon_export.p","rb"))
    self.days = pickle.load(open(direct + name + '_' + "days.p","rb"))
    self.irrprofile = pickle.load(open(direct + name + '_' + "irradiance_profile.p","rb"))
    self.fphyto = pickle.load(open(direct + name + '_' + "phyto_final_spectrum.p","rb"))
    self.phyto_total_mean_biomass = pickle.load(open(direct + name + '_' + "phytoplankton_total_mean_biomass.p","rb"))
    self.zoo_total_mean_biomass = pickle.load(open(direct + name + '_' + "zooplankton_total_mean_biomass.p","rb"))
    self.mean_shannon_p = pickle.load(open(direct + name + '_' + "mean_shannon_p.p","rb"))
    self.mean_shannon_z = pickle.load(open(direct + name + '_' + "mean_shannon_z.p","rb"))
    self.mean_N = pickle.load(open(direct + name + '_' + "mean_N.p","rb"))
    self.mean_detritus = pickle.load(open(direct + name + '_' + "mean_detritus.p","rb"))
    
    self.total_ave_carbon_export = (1.0/self.io.simulation_length)*(self.total_carbon_export.sum(axis=0) - self.total_carbon_export[0])*self.io.time_step*self.io.print_step
    self.avez_p_biomass = (1/self.h)*self.io.dz*self.phyto_total_mean_biomass.sum()
    self.avez_z_biomass = (1/self.h)*self.io.dz*self.zoo_total_mean_biomass.sum()
    self.avez_shannon_p = (1/self.h)*self.io.dz*self.mean_shannon_p.sum()
    self.avez_shannon_z = (1/self.h)*self.io.dz*self.mean_shannon_z.sum()
    self.avez_N = (1/self.h)*self.io.dz*self.mean_N.sum()
    #self.total_integrated_carbon_export = (self.total_carbon_export.sum(axis=0) - self.total_carbon_export[0])*(self.io.time_step)*self.io.print_step
#######open world 
labe=["I=959.825439","I + 10%", "I + 20%", "I + 30%", "I - 10%", "I - 20%", "I - 30%"]   
sim1= splas.splas("./input_ps1/input0.json")
sim1.get_input()
load(sim1,dir_pickle,sim1.io.name)

sim2= splas.splas("./input_ps1/input1.json")
sim2.get_input()
load(sim2,dir_pickle,sim2.io.name)

sim3= splas.splas("./input_ps1/input2.json")
sim3.get_input()
load(sim3,dir_pickle,sim3.io.name)

sim4= splas.splas("./input_ps1/input3.json")
sim4.get_input()
load(sim4,dir_pickle,sim4.io.name)

sim5= splas.splas("./input_ps1/input4.json")
sim5.get_input()
load(sim5,dir_pickle,sim5.io.name)

sim6= splas.splas("./input_ps1/input5.json")
sim6.get_input()
load(sim6,dir_pickle,sim6.io.name)

sim7= splas.splas("./input_ps1/input6.json")
sim7.get_input()
load(sim7,dir_pickle,sim7.io.name)

####input closed world changing irradiance in time 
# labe=["I=959.825439","I + 10%", "I + 20%", "I + 30%", "I - 10%", "I - 20%", "I - 30%"]
# sim1= splas.splas("./input_closed_world_tv/input0.json")
# sim1.get_input()
# load(sim1,dir_pickle,sim1.io.name)

# sim2= splas.splas("./input_closed_world_tv/input1.json")
# sim2.get_input()
# load(sim2,dir_pickle,sim2.io.name)


# sim3= splas.splas("./input_closed_world_tv/input2.json")
# sim3.get_input()
# load(sim3,dir_pickle,sim3.io.name)

# sim4= splas.splas("./input_closed_world_tv/input3.json")
# sim4.get_input()
# load(sim4,dir_pickle,sim4.io.name)

# sim5= splas.splas("./input_closed_world_tv/input4.json")
# sim5.get_input()
# load(sim5,dir_pickle,sim5.io.name)

# sim6= splas.splas("./input_closed_world_tv/input5.json")
# sim6.get_input()
# load(sim6,dir_pickle,sim6.io.name)

# sim7= splas.splas("./input_closed_world_tv/input6.json")
# sim7.get_input()
# load(sim7,dir_pickle,sim7.io.name)

##########
########## closed world I changing
#labe=["I=959.825439","I + 10%", "I + 20%", "I + 30%", "I - 10%", "I - 20%", "I - 30%"]
# dir_pickle = "./SPLAS_PICKLE/"
# sim1= splas.splas("./input_cw_changing/input0.json")
# sim1.get_input()
# load(sim1,dir_pickle,sim1.io.name)


# sim2= splas.splas("./input_cw_changing/increasing10.json")
# sim2.get_input()
# load(sim2,dir_pickle,sim2.io.name)


# sim3= splas.splas("./input_cw_changing/increasing20.json")
# sim3.get_input()
# load(sim3,dir_pickle,sim3.io.name)

# sim4= splas.splas("./input_cw_changing/increasing30.json")
# sim4.get_input()
# load(sim4,dir_pickle,sim4.io.name)

# sim5= splas.splas("./input_cw_changing/decreasing10.json")
# sim5.get_input()
# load(sim5,dir_pickle,sim5.io.name)

# sim6= splas.splas("./input_cw_changing/decreasing20.json")
# sim6.get_input()
# load(sim6,dir_pickle,sim6.io.name)

# sim7= splas.splas("./input_cw_changing/decreasing30.json")
# sim7.get_input()
# load(sim7,dir_pickle,sim7.io.name)


########### default closed world 
# labe=["I=959.825439","I + 10%", "I + 20%", "I + 30%", "I - 10%", "I - 20%", "I - 30%"]
# dir_pickle = "./SPLAS_PICKLE/"
# sim1= splas.splas("./input_closed_world/input0.json")
# sim1.get_input()
# load(sim1,dir_pickle,sim1.io.name)


# sim2= splas.splas("./input_closed_world/input1.json")
# sim2.get_input()
# load(sim2,dir_pickle,sim2.io.name)


# sim3= splas.splas("./input_closed_world/input2.json")
# sim3.get_input()
# load(sim3,dir_pickle,sim3.io.name)

# sim4= splas.splas("./input_closed_world/input3.json")
# sim4.get_input()
# load(sim4,dir_pickle,sim4.io.name)

# sim5= splas.splas("./input_closed_world/input4.json")
# sim5.get_input()
# load(sim5,dir_pickle,sim5.io.name)

# sim6= splas.splas("./input_closed_world/input5.json")
# sim6.get_input()
# load(sim6,dir_pickle,sim6.io.name)

# sim7= splas.splas("./input_closed_world/input6.json")
# sim7.get_input()
# load(sim7,dir_pickle,sim7.io.name)

############
#normal distribution 
# labe=["I=959.825439","I + 20%", "I + 20% sd10", "I + 20% sd20", "I - 20%", "I - 20% sd10", "I - 20% sd20"]
# sim1= splas.splas("./input_closed_world/input0.json")
# sim1.get_input()
# load(sim1,dir_pickle,sim1.io.name)

# sim2= splas.splas("./input_closed_world/input2.json")
# sim2.get_input()
# load(sim2,dir_pickle,sim2.io.name)

# sim3= splas.splas("./input_cw_changing/norm10.json")
# sim3.get_input()
# load(sim3,dir_pickle,sim3.io.name)

# sim4= splas.splas("./input_cw_changing/norm20.json")
# sim4.get_input()
# load(sim4,dir_pickle,sim4.io.name)

# sim5= splas.splas("./input_closed_world/input5.json")
# sim5.get_input()
# load(sim5,dir_pickle,sim5.io.name)

# sim6= splas.splas("./input_cw_changing/mnorm10.json")
# sim6.get_input()
# load(sim6,dir_pickle,sim6.io.name)

# sim7= splas.splas("./input_cw_changing/mnorm20.json")
# sim7.get_input()
# load(sim7,dir_pickle,sim7.io.name)


#plots
calendar=['0','1yr','2yr','3yr','4yr','5yr','6yr','7yr','8yr','9yr','10yr','11yr']
#labe=["I=959.825439","I + 10%", "I + 20%", "I + 30%", "I - 10%", "I - 20%", "I - 30%"]
##labe=[sim1.io.name,sim2.io.name,sim3.io.name,sim4.io.name,sim5.io.name,sim6.io.name,sim7.io.name]
nspec(sim1,sim2,sim3,sim4,sim5,sim6,sim7)
biomass(sim1,sim2,sim3,sim4,sim5,sim6,sim7)
dbm_depth(sim1,sim2,sim3,sim4,sim5,sim6,sim7)
diversity(sim1,sim2,sim3,sim4,sim5,sim6,sim7)
diversity_diff(sim1,sim2,sim3,sim4,sim5,sim6,sim7)
N_profile(sim1,sim2,sim3,sim4,sim5,sim6,sim7)
for t in enumerate(range(1, sim2.tdim, int(sim2.tdim/10))):
   screenshot(sim1,sim2,sim3,sim4,sim5,sim6,sim7)
   screenshot_diff(sim1,sim2,sim3,sim4,sim5,sim6,sim7)
   
   #size(sim1,sim2,sim3,sim4,sim5,sim6,sim7)
size_spectrum(sim1,sim2,sim3,sim4,sim5,sim6,sim7)
final_size(sim1,sim2,sim3,sim4,sim5,sim6,sim7)
light_diff(sim1,sim2,sim3,sim4,sim5,sim6,sim7)
carbon_out(sim1,sim2,sim3,sim4,sim5,sim6,sim7)
carbon_light(sim1,sim2,sim3,sim4,sim5,sim6,sim7)
sens_table(sim1,sim2,sim3,sim4,sim5,sim6,sim7)
dbm_layer=np.unique(dbm_layer)
dbm_layer=np.sort(dbm_layer)
print(dbm_layer) #unique + sort
dbm_zoom(sim1,sim2,sim3,sim4,sim5,sim6,sim7)
winners(sim1,sim2,sim3,sim4,sim5,sim6,sim7)
time_diff_N(sim1,sim2,sim3,sim4,sim5,sim6,sim7)
time_N(sim1,sim2,sim3,sim4,sim5,sim6,sim7)
#959,825439
plt.close('all')