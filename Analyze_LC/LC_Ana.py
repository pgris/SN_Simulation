import numpy as np
import os
import glob
import cPickle as pkl
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import UnivariateSpline
from SN_Rate import *

class LC_Ana:
    def __init__(self,data_dict={},zmin=0.,zmax=1.2,bin_z=0.01,data_dir='../../Fitted_Light_Curves_snsim'):

        self.data={}
        self.data_dir=data_dir
        self.Load_Data(data_dict)

        self.nsn_tot=0.
        self.err_tot=0.
        self.nsn_theo=0
        self.err_theo=0
        self.res_nsn=[]

        self.ms=['o','o','s','s','.','.','^','^','<','<']
        self.color=['b','r','b','r','b','r','b','r','b','r']
        data_sel={}

        for key, vals in self.data.items():
            print key
            data_sel[key]=self.Select_Data(vals)
            print len(vals),len(data_sel[key])
        
        self.Plot_Efficiency(self.data,data_sel,data_dict,'z',zmin,zmax,bin_z)
        self.Plot_N_SN(self.data,data_sel,data_dict,'z',zmin,zmax,bin_z,cumul=False)

        plt.show()

    def Load_Data(self, data_dict={}):

        for key,val in data_dict.items():
            dirmeas=self.data_dir+'/'+val.thedir+'/'+str(val.fieldid)+'/Season_'+str(val.season)
            sum_file=dirmeas+'/'+val.fieldname+'_'+str(val.fieldid)+'_X1_'+str(val.X1)+'_C_'+str(val.Color)+'_all.pkl'
            
            if os.path.exists(sum_file):
                pkl_file = open(sum_file,'rb')
                print 'loading',sum_file
                loaded_file=pkl.load(pkl_file)
                #print 'done',key,tot_resu[key],tot_resu[key].dtype,sum_file
                #idx = loaded_file['z']>=self.zmin and loaded_file['z']<self.zmax
                self.data[key]=loaded_file
            else:
                files = glob.glob(dirmeas+'/'+val.fieldname+'_'+str(val.fieldid)+'*_X1_'+str(val.X1)+'_C_'+str(val.Color)+'*.pkl')
                for fi in files:
                    pkl_file = open(fi,'rb')
                    loaded=pkl.load(pkl_file)
                    print 'loading',fi
                    #if 'mbsim' in loaded.dtype.names:
                    if not key in self.data.keys():
                        self.data[key]=loaded
                    else:
                        self.data[key]=np.vstack((self.data[key],loaded))

                pkl_out = open(sum_file,'wb')
                
                pkl.dump(self.data[key], pkl_out)
                
                pkl_out.close()


    def Select_Data(self,data):

        sel=data.copy()
        logical_test={}
    
        #rint 'before',len(sel)
        for band in 'grizy':
            logical_test[band]=np.logical_and(sel['N_bef_'+band]>=1,sel['N_aft_'+band]>=1)

        logical_and_g_r=np.logical_and(logical_test['g'],logical_test['r'])
        logical_and_r_i=np.logical_and(logical_test['r'],logical_test['i'])
        logical_and_i_z=np.logical_and(logical_test['i'],logical_test['z'])
        logical_and_z_y=np.logical_and(logical_test['z'],logical_test['y'])
                    
        #print 'logical',len(sel)
        sel=sel[np.where(np.logical_and(sel['status']=='go_fit',sel['fit_status']=='fit_ok'))]
        """
        for val in sel:
            print val['phase_first'],(val['DayMax']-val['salt2.T0'])/(1.+val['z'])
        """
        sel=sel[np.where(np.logical_and(sel['phase_first']<=-5,sel['phase_last']>=20))]
        
        #print 'all',len(sel)
        return sel

    def Histo_ratio(self,sela,selb,varname,zmin,zmax,bin_z):

        range=[zmin,zmax]
        bin_width=bin_z
        print 'alors pal',range,bin_width
        num_bins=int((range[1]-range[0])/bin_width)
        #range=[0.0,1.2]
        
        hista, bin_edgesa = np.histogram(sela[varname],bins=num_bins,range=range,weights=sela['X1_weight']*sela['Color_weight'])
        histb, bin_edgesb = np.histogram(selb[varname],bins=num_bins,range=range,weights=selb['X1_weight']*selb['Color_weight'])
        bin_center = (bin_edgesa[:-1] + bin_edgesa[1:]) / 2

        ratio=[]
        ratio_err=[]
        norm=[]
        norm_err=[]
        
        for a,b in zip(histb,hista):
            if b==0:
                ratio.append(a)
                ratio_err.append(0)
            else:
                effi=float(a)/float(b)
                ratio.append(effi)
                ratio_err.append(np.sqrt(1.-effi)*effi/np.sqrt(float(b)))
            eff=float(a)/float(np.sum(hista))
            norm.append(eff)
            norm_err.append(np.sqrt(1.-eff)*eff/np.sqrt(float(np.sum(hista))))
    
        return bin_center, ratio, ratio_err,norm,norm_err

    def Plot_Efficiency(self,data,dict_sel,data_dict,varname,min_bin,max_bin,delta_bin,obs=None):
        tot_label=[]
        figc, axc = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        for key, vals in data.items():
            self.Plot_Efficiency_Indiv(axc,vals,dict_sel[key],varname,tot_label,key,data_dict[key].colorfig,data_dict[key].season,min_bin,max_bin,delta_bin,obs)

    def Plot_Efficiency_Indiv(self,axc,sela,selb,varname,tot_label,ll,color,season,min_bin,max_bin,delta_bin,obs=None):

        filtercolors = {'u':'b', 'g':'c', 'r':'g', 'i':'y', 'z':'r', 'y':'m'}

        fontsize=12
        bin_center, ratio, ratio_err,norm,norm_err= self.Histo_ratio(sela,selb,varname,min_bin,max_bin,delta_bin)
        
        ll='Y'+str(season+1)
        axc.errorbar(bin_center,ratio, yerr=ratio_err,marker=self.ms[season], mfc=self.color[season], mec=self.color[season], ms=8, linestyle='-',color=color,label=ll)

    def Plot_N_SN(self,data,dict_sel,data_dict,varname,min_bin,max_bin,delta_bin,cumul):

        tot_label=[]
        figb, axb = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        for key, vals in data.items():
            self.Plot_N_SN_Indiv(axb,vals,dict_sel[key],varname,tot_label,key,data_dict[key].colorfig,data_dict[key].season,data_dict[key].fieldname,data_dict[key].fieldid,min_bin,max_bin,delta_bin,cumul)

        axb.set_xlabel('z')
        if cumul is True:
            axb.set_ylabel('N$_{SN}$ < z')
        else:
            axb.set_ylabel('N$_{SN}$')

        axb.set_xlim([min_bin,max_bin+0.01])

    def Plot_N_SN_Indiv(self,axb,sela,selb,varname,tot_label,ll,color,season,fieldname,fieldid,min_bin,max_bin,delta_bin,cumul=False):

        bin_center, ratio, ratio_err,norm,norm_err= self.Histo_ratio(sela,selb,varname,min_bin,max_bin,delta_bin)

        effi = interpolate.interp1d(bin_center,ratio)

        rate_name='Perrett'
        sn_rate=SN_Rate(rate=rate_name,duration=(selb['duration'][0]+0.)/365.25)

        #zz,rate,err_rate,nsn,err_nsn=sn_rate(self.zmin,self.zmax,self.bin_z)
        zz,rate,err_rate,nsn,err_nsn=sn_rate(bins=bin_center)
        
        print 'Nsn',np.sum(nsn),zz,bin_center,nsn,err_nsn,np.power(np.sum(np.power(err_nsn,2.)),0.5)
        
        nsn_season = interpolate.interp1d(zz,nsn)
        err_nsn_season = interpolate.interp1d(zz,err_nsn)

        combi=nsn_season(zz)*effi(zz)
        #yerr_combi=np.power(np.power(ratio_err*nsn_season(zz),2.)+np.power(0.1*nsn_season(zz)*effi(zz),2.),0.5)
        yerr_combi=np.power(np.power(ratio_err*nsn_season(zz),2.)+np.power(err_nsn_season(zz)*effi(zz),2.),0.5)

        print 'Number of SN Ia',season,np.sum(combi),np.power(np.sum(np.power(yerr_combi,2.)),0.5)
        N_sn=np.sum(combi)
        err_N_sn=np.power(np.sum(np.power(yerr_combi,2.)),0.5)
        self.nsn_tot+=np.sum(combi)
        self.err_tot+=np.power(err_N_sn,2.)
        self.nsn_theo+=np.sum(nsn_season(zz))
        self.err_theo+=np.sum(np.power(err_nsn_season(zz),2.))
        self.res_nsn.append((fieldname,fieldid,season,np.sum(combi),err_N_sn,np.sum(nsn_season(zz)),np.sum(np.power(err_nsn_season(zz),2.))))
        #axc[1].plot(zz,nsn,'b+')
        #axc[1].plot(zz,nsn_season(zz),marker='.', mfc='black', mec='black', ms=8, linestyle='-',color='black')

        if not cumul:
            axb.errorbar(zz,nsn_season(zz),yerr=err_nsn_season(zz),marker=self.ms[season], mfc=self.color[season], mec=self.color[season], ms=8, linestyle='--',color='black')
            axb.errorbar(zz,combi,yerr=yerr_combi,marker=self.ms[season], mfc=self.color[season], mec=self.color[season], ms=8, linestyle='-',color=color,label=ll)
        else:
            axb.errorbar(zz,np.cumsum(nsn_season(zz)),marker=self.ms[season], mfc=self.color[season], mec=self.color[season], ms=8, linestyle='--',color='black')
            axb.errorbar(zz,np.cumsum(combi),marker=self.ms[season], mfc=self.color[season], mec=self.color[season], ms=8, linestyle='-',color=color,label=ll)
