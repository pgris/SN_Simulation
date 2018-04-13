import numpy as np
from scipy.interpolate import griddata
from Load_X1_C import Load_X1_Color
from scipy import interpolate
from scipy.interpolate import UnivariateSpline
from SN_Rate import *
import os
import cPickle as pkl
from astropy.table import Table,Column,vstack
import multiprocessing
import time

class Calc_NSN:

    def __init__(self,fieldname,fieldid,season,effi,simu_name,selection_name):

        self.fieldname=fieldname
        self.fieldid=fieldid
        self.season=season
        self.outname='N_SN_'+simu_name+'_'+selection_name
    
        #print('there www',np.unique(effi['duration']),np.unique(effi['season']))
        self.effi,effi_all=self.Interpolate(effi)

        self.nsn_tot=0.
        self.err_tot=0.
        self.nsn_theo=0
        self.err_theo=0
        self.res_nsn=[]

        self.tab_theo=self.Get_NSN_old()

        print(effi_all)

        """
        tab_tot=Table()
        for x1_color in np.unique(effi_all[['X1','Color']]):
        #for x1_color in [(0.012,0.015)]:
            #print('rrr',x1_color)
            idx = (np.abs(effi_all['X1']-x1_color[0])<1.e-5)&(np.abs(effi_all['Color']-x1_color[1])<1.e-5)
            sel=effi_all[idx]
            #print('hello',sel)
            tab=self.Get_NSN(sel)
            tab_tot=vstack([tab_tot,tab])
        """

        tab_tot=self.Process_Multiproc(effi_all)
        nsn_tot=0.
        err_nsn_tot=0.

        r=[]
        for z in np.unique(tab_tot['z']):
            idx = np.abs(tab_tot['z']-z)<1.e-6
            sel=tab_tot[idx]
            r.append((z,np.sum(sel['nsn']),np.sqrt(np.sum(np.power(sel['err_nsn'],2.))),np.sum(sel['nsn_theo']),np.sqrt(np.sum(np.power(sel['err_theo'],2.)))))
            
        nsn_tot=np.sum(tab_tot['nsn'])
        #err_nsn_tot=np.sum(1./np.sum(1./np.power(tab_tot['err_nsn'],2.)))
        err_nsn_tot=np.sum(tab_tot['err_nsn']**2)
        err_theo=np.sqrt(np.sum(self.tab_theo['err_theo']**2.))
        nsn_theo=np.sum(self.tab_theo['nsn_theo'])

        self.tab_sum=np.rec.fromrecords(r,names=['z','nsn','err_nsn','nsn_theo','err_theo'])
        """
        tab_sum=Table()
        tab_sum.add_column(Column(r[:][0],name='z'))
        tab_sum.add_column(Column(r[:][1],name='nsn'))
        tab_sum.add_column(Column(r[:][2],name='err_nsn'))
        self.tab_sum=tab_sum
        """
        
        ro=[]
        ro.append((self.fieldname,self.fieldid,self.season,nsn_tot,np.sqrt(err_nsn_tot),nsn_theo,err_theo))
        print('ro',ro)
        rate_nsn=np.rec.fromrecords(ro,names=['fieldname','fieldid','season','n_sn_detected','err_detected','n_sn_expected','err_expected'])
        self.Check_Dir(self.outname)
        pkl_out=open(self.outname+'/N_SN_'+self.fieldname+'_'+str(self.fieldid)+'_Season_'+str(self.season)+'.pkl','wb')
        pkl.dump(rate_nsn,pkl_out)
        pkl_out.close()

    def Process_Multiproc(self,effi_all):

        n_per_batch=10
        x1_color=np.unique(effi_all[['X1','Color']])
        #x1_color=[(0.012,0.015),(0.012,0.015)]
        inter=range(0,len(x1_color),n_per_batch)
        inter=np.append(inter,len(x1_color))
        
        restot=None
        time_begin=time.time()
        for jo in range(len(inter)-1):

            ida=inter[jo]
            idb=inter[jo+1]
            result_queue = multiprocessing.Queue()
            #print(jo,time.time()-time_begin)
            for j in range(ida,idb):
                idx = (np.abs(effi_all['X1']-x1_color[j][0])<1.e-5)&(np.abs(effi_all['Color']-x1_color[j][1])<1.e-5)
                sel=effi_all[idx]
                #print('radec',j,len(sel),sel)
                p=multiprocessing.Process(name='Subprocess-'+str(j),target=self.Get_NSN,args=(sel,1.1,j,result_queue))
                p.start()
 
            resultdict = {}
    
            for j in range(ida,idb):
                resultdict.update(result_queue.get())

            for p in multiprocessing.active_children():
                p.join()

            for j in range(ida,idb):
                if restot is None:
                    restot=resultdict[j]
                else:
                    restot=vstack([restot,resultdict[j]])

        print('End of data loading',time.time()-time_begin)
        return restot



    def Interpolate(self,sel):

        grid_x1_c=Load_X1_Color(-999.0,-999.0,1,name_low_z='Dist_X1_Color_jla_low_z.txt',name_high_z='Dist_X1_Color_jla_high_z.txt').tab
        
        duration=sel['duration'][0]
        r=[]
        rall=[]
        #print('before',np.unique(sel['z']))
        for z in np.unique(sel['z']):
            idxb= np.abs(sel['z']-z)<1.e-5
            sel_z=sel[idxb]
            #print('hello',z,len(sel_z))
            if z < 0.1:
                grid=grid_x1_c['low_z']
            else:
                grid=grid_x1_c['high_z']

            if len(sel_z) >= 2:
               
                grid_res=griddata(np.asarray(sel_z[['X1','Color']],dtype=np.float),np.asarray(sel_z['effi'],dtype=np.float),np.asarray(grid[['x1','c']], dtype=np.float), method='nearest')
                grid_res_err=griddata(np.asarray(sel_z[['X1','Color']],dtype=np.float),np.asarray(sel_z['err_effi'],dtype=np.float),np.asarray(grid[['x1','c']], dtype=np.float), method='nearest')
                grid_duration=griddata(np.asarray(sel_z[['X1','Color']],dtype=np.float),np.asarray(sel_z['duration_z'],dtype=np.float),np.asarray(grid[['x1','c']], dtype=np.float), method='nearest')
            #print grid_res                                                                                                          

                effi_comb=np.sum([grid['weight_x1'][i]*grid['weight_c'][i]*grid_res[i] for i in range(len(grid_res))])
                effi_comb_err=np.sum([grid['weight_x1'][i]*grid['weight_c'][i]*grid_res_err[i] for i in range(len(grid_res_err))])
            #print effi_comb,effi_comb_err 
                duration_z=np.mean(sel_z['duration_z'])
                #print('duration',duration_z)
                r.append((self.fieldname,self.fieldid,self.season,z,effi_comb,effi_comb_err,duration,duration_z,-999.0,-999.0))
                #print('hello',grid[['x1','c']])
                for i in range(len(grid_res)):
                    #print('helli',grid[['x1','c']][i])
                    effi=np.asscalar(grid_res[i])
                    err_effi=np.asscalar(grid_res_err[i])
                    duration_z=np.asscalar(grid_duration[i])
                    weight_x1=np.asscalar(grid['weight_x1'][i])
                    weight_c=np.asscalar(grid['weight_c'][i])
                    
                    rall.append((z,effi,err_effi,duration_z,grid[['x1','c']][i][0][0],grid[['x1','c']][i][0][1],weight_x1,weight_c))
        
            else:
                r.append((self.fieldname,self.fieldid,self.season,z,0.,0.,duration,duration_z,-999.0,-999.0))
                for i in range(len(grid_res)):
                    duration_z=np.asscalar(grid_duration[i])
                    weight_x1=np.asscalar(grid['weight_x1'][i])
                    weight_c=np.asscalar(grid['weight_c'][i])
                    rall.append((z,0.,0.,duration_z,grid[['x1','c']][i][0][0],grid[['x1','c']][i][0][1],weight_x1,weight_c))
            
        effi_pond=np.rec.fromrecords(r,names=['fieldname','fieldid','season','z','effi','err_effi','duration','duration_z','X1','Color'])
        #print(rall)
        effi_all=np.rec.fromrecords(rall,names=['z','effi','err_effi','duration_z','X1','Color','weight_x1','weight_c'])

        tab=Table()
        tab.meta=dict(zip(['fieldname','fieldid','season','duration'],[self.fieldname,self.fieldid,self.season,duration]))
        #print('boo',effi_all['z'],len(effi_all['z']),len(rall))
        for name in effi_all.dtype.names:
            tab.add_column(Column(effi_all[name],name=name))

        return effi_pond,tab

    def Get_NSN(self,effi,zmax=1.1,j=0,output_q=None):

        tab=effi
        #print(effi)
        #tab.sort(order='z')
        #print 'ooo',tab['z'],type(tab)                                                                                                                             
        effi = interpolate.interp1d(tab['z'],tab['effi'])
        #ratio_err=tab['err_effi']
        ratio_err = interpolate.interp1d(tab['z'],tab['err_effi'])
        #print('oooo',tab['z'])
        duration_z=interpolate.interp1d(tab['z'],tab['duration_z'])
        rate_name='Perrett'
        sn_rate=SN_Rate(rate=rate_name,duration=0.)
        #print('attention',self.season,tab['duration'])

        #zz,rate,err_rate,nsn,err_nsn=sn_rate(self.zmin,self.zmax,self.bin_z)                                                                                    
        zzb=np.arange(np.min(tab['z']),np.max(tab['z']),0.01)
        #print('test duration',tab['duration_z'])
        zz,rate,err_rate,nsn,err_nsn=sn_rate(bins=zzb,duration_z=duration_z)

        #print('Nsn',self.fieldid,self.season,np.sum(nsn),err_nsn,np.power(np.sum(np.power(err_nsn,2.)),0.5))

        nsn_season = interpolate.interp1d(zz,nsn)
        err_nsn_season = interpolate.interp1d(zz,err_nsn)
        rate_zz=interpolate.interp1d(zz,rate)

        idx = zz < zmax
        zz=zz[idx]

        coeff=np.unique(tab['weight_x1'])*np.unique(tab['weight_c'])
        #coeff=0.5
        #combi=nsn_season(zz)*effi(zz)*coeff
        combi=nsn_season(zz)*effi(zz)*coeff
        #yerr_combi=np.power(np.power(ratio_err*nsn_season(zz),2.)+np.power(0.1*nsn_season(zz)*effi(zz),2.),0.5)                                                      
        
        #yerr_combi=np.power(np.power(ratio_err(zz)*rate_zz(zz)*coeff,2.)+np.power(err_nsn_season(zz)*effi(zz)*coeff,2.),0.5)
        yerr_combi=np.power(np.power(ratio_err(zz)*rate_zz(zz)*coeff,2.)+np.power(err_nsn_season(zz)*effi(zz)*coeff,2.),0.5)

        combi_th=nsn_season(zz)*coeff
        yerr_combi_th=np.power(np.power(ratio_err(zz)*rate_zz(zz)*coeff,2.)+np.power(err_nsn_season(zz)*coeff,2.),0.5)

        #print('Nsn bis',coeff,np.sum(combi),np.sqrt(np.sum(yerr_combi*yerr_combi)))
        #yerr_combi=err_nsn_season(zz)*coeff
        #print('pal',np.sum(combi),np.sqrt(np.sum(yerr_combi**2)))
        tab_out=Table()
        tab_out.add_column(Column(zz,name='z'))
        tab_out.add_column(Column(combi,name='nsn'))
        tab_out.add_column(Column(yerr_combi,name='err_nsn'))
        tab_out.add_column(Column(combi_th,name='nsn_theo'))
        tab_out.add_column(Column(yerr_combi_th,name='err_theo'))
        tab_out.add_column(Column([np.unique(tab['weight_x1'])]*len(zz),name='weight_x1'))
        tab_out.add_column(Column([np.unique(tab['weight_c'])]*len(zz),name='weight_c'))
        tab_out.add_column(Column([np.unique(tab['X1'])]*len(zz),name='X1'))
        tab_out.add_column(Column([np.unique(tab['Color'])]*len(zz),name='Color'))

        #print('finished',j)
        if output_q is not None:
            #print('returning',j)
            output_q.put({j : tab_out})
        else:
            return tab_out
        """
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
        """

    def Get_NSN_old(self,zmax=1.1):

        tab=self.effi.copy()
        tab.sort(order='z')
        #print 'ooo',tab['z'],type(tab)                                                                                                                             
        effi = interpolate.interp1d(tab['z'],tab['effi'])
        #ratio_err=tab['err_effi']
        ratio_err = interpolate.interp1d(tab['z'],tab['err_effi'])
        duration_z=interpolate.interp1d(tab['z'],tab['duration_z'])
        rate_name='Perrett'
        sn_rate=SN_Rate(rate=rate_name,duration=(tab['duration'][0]+0.0)/365.25)
        #print('attention',self.season,tab['duration'])

        #zz,rate,err_rate,nsn,err_nsn=sn_rate(self.zmin,self.zmax,self.bin_z)                                                                                    
        zzb=np.arange(np.min(tab['z']),np.max(tab['z']),0.01)
        #zzb=np.arange(0.01,0.3,0.05)
        #print('test duration',tab['duration_z'])
        zz,rate,err_rate,nsn,err_nsn=sn_rate(bins=zzb,duration_z=duration_z)

        #print('Nsn',self.fieldid,self.season,np.sum(nsn),zz,tab['z'],nsn,err_nsn,np.power(np.sum(np.power(err_nsn,2.)),0.5))

        nsn_season = interpolate.interp1d(zz,nsn)
        err_nsn_season = interpolate.interp1d(zz,err_nsn)

       
        idx = zz < zmax
        zz=zz[idx]

        combi=nsn_season(zz)*effi(zz)
        #yerr_combi=np.power(np.power(ratio_err*nsn_season(zz),2.)+np.power(0.1*nsn_season(zz)*effi(zz),2.),0.5)                                                      

        yerr_combi=np.power(np.power(ratio_err(zz)*nsn_season(zz),2.)+np.power(err_nsn_season(zz)*effi(zz),2.),0.5)

        #print('Number of SN Ia',self.season,np.sum(combi),np.power(np.sum(np.power(yerr_combi,2.)),0.5))
        """
        r=[]
        nsn_tot=np.sum(combi)
        err_tot=np.power(err_N_sn,2.)
        nsn_theo=np.sum(nsn_season(zz))
        err_theo=np.sum(np.power(err_nsn_season(zz),2.))
        r.append((self.fieldname,self.fieldid,int(self.season),nsn_tot,np.sqrt(err_tot),nsn_theo,np.sqrt(err_theo)))
        """
        tab=Table()
        tab.add_column(Column(zz,name='z'))
        tab.add_column(Column(combi,name='nsn'))
        tab.add_column(Column(yerr_combi,name='err_nsn'))
        tab.add_column(Column(nsn_season(zz),name='nsn_theo'))
        tab.add_column(Column(err_nsn_season(zz),name='err_theo'))
        return tab

        """
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
        """

    def Check_Dir(self,dirout):
        if not os.path.isdir(dirout) :
            os.makedirs(dirout)

    def zpercentile(self):
        
        #print('hello',self.zz,np.cumsum(self.combi)/np.cumsum(self.combi)[-1])
        f = interpolate.interp1d(self.zz,np.cumsum(self.combi))
        zint=np.arange(0.02,1.2,0.01)
        num_sn=f(zint)
        num_sn/=num_sn[-1]
        #print(type(num_sn/num_sn[-1]))
        idxa = (np.abs(num_sn-0.90)).argmin()
        idxb = (np.abs(num_sn-0.95)).argmin()
        #print(num_sn[idx],zint[idx])
        return zint[idxa],zint[idxb]

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

        #nsn_str= str(int(self.nsn_theo))+'$\pm$'+str(int(np.sqrt(self.err_theo)))
        """
        err_theo=np.sqrt(np.sum(self.tab_theo['err_theo']**2.))
        nsn_theo=np.sum(self.tab_theo['nsn_theo'])
        """
        err_theo=np.sqrt(np.sum(self.tab_sum['err_theo']**2.))
        nsn_theo=np.sum(self.tab_sum['nsn_theo'])
        nsn_str= str(int(nsn_theo))+'$\pm$'+str(int(err_theo))
        nsn_obs=np.sum(self.tab_sum['nsn'])
        err_obs=np.sqrt(np.sum(self.tab_sum['err_nsn']**2))
        if self.season < 9:
            #ll='Y'+str(self.season+1)+'   - N$_{SN Ia}$ = '+str(int(self.nsn_tot))+' $\pm$ '+str(int(np.sqrt(self.err_tot)))+' / '+nsn_str
            #ll='Y'+str(self.season+1)+'   - N$_{SN Ia}$ = '+str(int(self.nsn_tot))+' $\pm$ '+str(int(np.sqrt(self.err_tot)))+' / '+nsn_str
            ll='Y'+str(self.season+1)+'  - N$_{SN Ia}$ = '+str(int(nsn_obs))+ '$\pm$'+ str(int(err_obs))+ ' / '+nsn_str
        else:
            ll='Y'+str(self.season+1)+' - N$_{SN Ia}$ = '+str(int(nsn_obs))+ '$\pm$'+ str(int(err_obs))+' / '+nsn_str
            #ll='Y'+str(self.season+1)+' - N$_{SN Ia}$ = '+str(int(self.nsn_tot))+' $\pm$ '+str(int(np.sqrt(self.err_tot)))+' / '+nsn_str
        if not cumul:
            #axb.errorbar(self.zz,self.nsn_season(self.zz),yerr=self.err_nsn_season(self.zz),marker=marker, mfc=color, mec=color, ms=8, linestyle='--',color='black')
            axb.errorbar(self.tab_theo['z'],self.tab_theo['nsn_theo'],yerr=self.tab_theo['err_theo'],marker=marker, mfc=color, mec=color, ms=8, linestyle='--',color='black')
            #tot_label.append(axb.errorbar(self.zz,self.combi,yerr=self.yerr_combi,marker=marker, mfc=color, mec=color, ms=8, linestyle='-',color=color,label=ll))
            tot_label.append(axb.errorbar(self.tab_sum['z'],self.tab_sum['nsn'],yerr=self.tab_sum['err_nsn'],marker=marker, mfc=color, mec=color, ms=8, linestyle='-',color=color,label=ll))
        else:
            #axb.errorbar(self.zz,np.cumsum(self.nsn_season(self.zz)),marker=marker, mfc=color, mec=color, ms=8, linestyle='--',color='black')
            axb.errorbar(self.tab_theo['z'],np.cumsum(self.tab_theo['nsn_theo']),marker=marker, mfc=color, mec=color, ms=8, linestyle='--',color='black')         
            #tot_label.append(axb.errorbar(self.zz,np.cumsum(self.combi),marker=marker, mfc=color, mec=color, ms=8, linestyle='-',color=color,label=ll))
            tot_label.append(axb.errorbar(self.tab_sum['z'],np.cumsum(self.tab_sum['nsn']),marker=marker, mfc=color, mec=color, ms=8, linestyle='-',color=color,label=ll))
            #print('aaaaaaa',self.tab_sum['z'],np.cumsum(self.tab_sum['nsn']))
        

        self.zmax_eff=np.max(self.tab_sum['z'])
