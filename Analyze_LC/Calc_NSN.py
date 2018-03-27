import numpy as np
from scipy.interpolate import griddata
from Load_X1_C import Load_X1_Color
from scipy import interpolate
from scipy.interpolate import UnivariateSpline
from SN_Rate import *
import os
import cPickle as pkl

class Calc_NSN:

    def __init__(self,fieldname,fieldid,season,effi,simu_name,selection_name):

        self.fieldname=fieldname
        self.fieldid=fieldid
        self.season=season
        self.outname='N_SN_'+simu_name+'_'+selection_name
    
        #print('there www',np.unique(effi['duration']),np.unique(effi['season']))
        self.effi=self.Interpolate(effi)

        self.nsn_tot=0.
        self.err_tot=0.
        self.nsn_theo=0
        self.err_theo=0
        self.res_nsn=[]

        self.Get_NSN()

    def Interpolate(self,sel):

        grid_x1_c=Load_X1_Color(-999.0,-999.0,1,name_low_z='Dist_X1_Color_jla_low_z.txt',name_high_z='Dist_X1_Color_jla_high_z.txt').tab
        
        duration=sel['duration'][0]
        r=[]
        for z in np.unique(sel['z']):
            idxb= sel['z']==z
            selb=sel[idxb]
            if len(selb) > 0. :
                if z < 0.1:
                    grid=grid_x1_c['low_z']
                else:
                    grid=grid_x1_c['high_z']
            grid_res=griddata(np.asarray(selb[['X1','Color']],dtype=np.float),np.asarray(selb['effi'],dtype=np.float),np.asarray\
(grid[['x1','c']], dtype=np.float), method='nearest')
            grid_res_err=griddata(np.asarray(selb[['X1','Color']],dtype=np.float),np.asarray(selb['err_effi'],dtype=np.float),np\
.asarray(grid[['x1','c']], dtype=np.float), method='nearest')
            #print grid_res                                                                                                          

            effi_comb=np.sum([grid['weight_x1'][i]*grid['weight_c'][i]*grid_res[i] for i in range(len(grid_res))])
            effi_comb_err=np.sum([grid['weight_x1'][i]*grid['weight_c'][i]*grid_res_err[i] for i in range(len(grid_res_err))])
            #print effi_comb,effi_comb_err                                                                                           
            r.append((self.fieldname,self.fieldid,self.season,z,effi_comb,effi_comb_err,duration,-999.0,-999.0))
        else:
            r.append((self.fieldname,self.fieldid,self.season,z,0.,0.,duration,-999.0,-999.0))

        return np.rec.fromrecords(r,names=['fieldname','fieldid','season','z','effi','err_effi','duration','X1','Color'])


    def Get_NSN(self):

        tab=self.effi.copy()
        tab.sort(order='z')
        #print 'ooo',tab['z'],type(tab)                                                                                                                             
        effi = interpolate.interp1d(tab['z'],tab['effi'])
        ratio_err=tab['err_effi']

        rate_name='Perrett'
        sn_rate=SN_Rate(rate=rate_name,duration=(tab['duration'][0]+0.0)/365.25)
        print('attention',self.season,tab['duration'])

        #zz,rate,err_rate,nsn,err_nsn=sn_rate(self.zmin,self.zmax,self.bin_z)                                                                                          
        zz,rate,err_rate,nsn,err_nsn=sn_rate(bins=tab['z'])

        print('Nsn',self.fieldid,self.season,np.sum(nsn),zz,tab['z'],nsn,err_nsn,np.power(np.sum(np.power(err_nsn,2.)),0.5))

        nsn_season = interpolate.interp1d(zz,nsn)
        err_nsn_season = interpolate.interp1d(zz,err_nsn)

        combi=nsn_season(zz)*effi(zz)
        #yerr_combi=np.power(np.power(ratio_err*nsn_season(zz),2.)+np.power(0.1*nsn_season(zz)*effi(zz),2.),0.5)                                                      

        yerr_combi=np.power(np.power(ratio_err*nsn_season(zz),2.)+np.power(err_nsn_season(zz)*effi(zz),2.),0.5)

        print('Number of SN Ia',self.season,np.sum(combi),np.power(np.sum(np.power(yerr_combi,2.)),0.5))

        self.zz=zz
        self.nsn_season=nsn_season
        self.err_nsn_season=err_nsn_season
        self.combi=combi
        self.yerr_combi=yerr_combi

        N_sn=np.sum(combi)
        err_N_sn=np.power(np.sum(np.power(yerr_combi,2.)),0.5)
        self.nsn_tot+=np.sum(combi)
        self.err_tot+=np.power(err_N_sn,2.)
        self.nsn_theo+=np.sum(nsn_season(zz))
        self.err_theo+=np.sum(np.power(err_nsn_season(zz),2.))
        self.res_nsn.append((self.fieldname,self.fieldid,int(self.season),self.nsn_tot,np.sqrt(self.err_tot),self.nsn_theo,np.sqrt(self.err_theo)))

        self.N_sn=N_sn
        self.err_N_sn=err_N_sn
        if len(self.res_nsn) > 0:
            rate_nsn=np.rec.fromrecords(self.res_nsn,names=['fieldname','fieldid','season','n_sn_detected','err_detected','n_sn_expected','err_expected'])
            self.Check_Dir(self.outname)
            pkl_out=open(self.outname+'/N_SN_'+self.fieldname+'_'+str(self.fieldid)+'_Season_'+str(self.season)+'.pkl','wb')
            pkl.dump(rate_nsn,pkl_out)
            pkl_out.close()
 
    def Check_Dir(self,dirout):
        if not os.path.isdir(dirout) :
            os.makedirs(dirout)

    def Plot_N_SN(self,axb,cumul=False,tot_label=[],color='r',marker='*'):

        """
        tab=self.effi.copy()
        tab.sort(order='z')
        #print 'ooo',tab['z'],type(tab)                                                                                              
        effi = interpolate.interp1d(tab['z'],tab['effi'])
        ratio_err=tab['err_effi']

        rate_name='Perrett'
        sn_rate=SN_Rate(rate=rate_name,duration=(tab['duration'][0]+0.0)/365.25)
        print 'attention',tab['duration']

        #zz,rate,err_rate,nsn,err_nsn=sn_rate(self.zmin,self.zmax,self.bin_z)                                                       
                                                                                                                           
        zz,rate,err_rate,nsn,err_nsn=sn_rate(bins=tab['z'])

        print('Nsn',self.fieldid,self.season,np.sum(nsn),zz,tab['z'],nsn,err_nsn,np.power(np.sum(np.power(err_nsn,2.)),0.5))

        nsn_season = interpolate.interp1d(zz,nsn)
        err_nsn_season = interpolate.interp1d(zz,err_nsn)

        combi=nsn_season(zz)*effi(zz)
        #yerr_combi=np.power(np.power(ratio_err*nsn_season(zz),2.)+np.power(0.1*nsn_season(zz)*effi(zz),2.),0.5)  
                                                                                                                                   
        yerr_combi=np.power(np.power(ratio_err*nsn_season(zz),2.)+np.power(err_nsn_season(zz)*effi(zz),2.),0.5)

        print('Number of SN Ia',self.season,np.sum(combi),np.power(np.sum(np.power(yerr_combi,2.)),0.5))
        N_sn=np.sum(combi)
        err_N_sn=np.power(np.sum(np.power(yerr_combi,2.)),0.5)
        self.nsn_tot+=np.sum(combi)
        self.err_tot+=np.power(err_N_sn,2.)
        self.nsn_theo+=np.sum(nsn_season(zz))
        self.err_theo+=np.sum(np.power(err_nsn_season(zz),2.))
        self.res_nsn.append((self.fieldname,self.fieldid,int(self.season),np.sum(combi),err_N_sn,np.sum(nsn_season(zz)),np.sum(np.power(err_nsn_season(zz),2.))))

        """

        nsn_str= str(int(self.nsn_theo))+'$\pm$'+str(int(np.sqrt(self.err_theo)))
        if self.season < 9:
            ll='Y'+str(self.season+1)+'   - N$_{SN Ia}$ = '+str(int(self.nsn_tot))+' $\pm$ '+str(int(np.sqrt(self.err_tot)))+' / '+nsn_str
        else:
            ll='Y'+str(self.season+1)+' - N$_{SN Ia}$ = '+str(int(self.nsn_tot))+' $\pm$ '+str(int(np.sqrt(self.err_tot)))+ ' / '+nsn_str

        if not cumul:
            axb.errorbar(self.zz,nsn_season(self.zz),yerr=self.err_nsn_season(self.zz),marker=marker, mfc=color, mec=color, ms=8, linestyle='--',color='black')
            tot_label.append(axb.errorbar(self.zz,self.combi,yerr=self.yerr_combi,marker=marker, mfc=color, mec=color, ms=8, linestyle='-',color=color,label=ll))
        else:
            axb.errorbar(self.zz,np.cumsum(self.nsn_season(self.zz)),marker=marker, mfc=color, mec=color, ms=8, linestyle='--',color='black')
            tot_label.append(axb.errorbar(self.zz,np.cumsum(self.combi),marker=marker, mfc=color, mec=color, ms=8, linestyle='-',color=color,label=ll))
        self.zmax_eff=np.max(self.zz)
