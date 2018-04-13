import numpy as np
import os
import glob
import cPickle as pkl
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import UnivariateSpline
from SN_Rate import *
import h5py
from astropy.table import Table,vstack,Column
from scipy.spatial import distance
from scipy.interpolate import griddata
from Load_X1_C import Load_X1_Color
from matplotlib import colors as mcolors
import multiprocessing
import time
import os

class LC_Ana:
    def __init__(self,fieldname,fieldid,season,zrange,x1_color,main_dir='/sps/lsst/data/dev/pgris',data_dir='../../Fitted_Light_Curves_snsim',selection={},simu_name='',Force_reprocess=False,add_covmb=False,alpha=0.14,beta=3.1):

        self.fieldname=fieldname
        self.fieldid=fieldid
        self.season=season
        self.zrange=zrange
        self.x1_color=x1_color
        self.data_dir=main_dir+'/'+data_dir
        self.selection=selection
        #self.selection_name=selection_name
        self.simu_name=simu_name
        self.Force_reprocess=Force_reprocess
        self.add_covmb=add_covmb
        self.alpha=alpha
        self.beta=beta

        self.nsn_tot=0.
        self.err_tot=0.
        self.nsn_theo=0
        self.err_theo=0
        self.res_nsn=[]

        #self.Check_Zero_Efficiencies(tot_data)
        #self.efficiencies=self.Get_Efficiencies(tot_data)
        
        #self.Plot_Effi_new()


    def Load_Dist_Error(self):

        self.tot_data=self.Select_all(self.Load_Data_Multiproc(),self.selection)

    def Effi(self,preselection):
        self.preselection=preselection
        self.efficiencies={}
        torepro={}
        load_data=False
        for selection_name, val in self.selection.items():
            self.effi_dir='Effi'+'_'+self.simu_name+'_'+selection_name
            self.Check_Dir(self.effi_dir)
            effi_name=self.effi_dir+'/Effi_'+self.fieldname+'_'+str(self.fieldid)+'_'+str(self.season)+'_'+str(self.x1_color[0])+'_'+str(self.x1_color[1])+'.npy'
            if not os.path.isfile(effi_name) or self.Force_reprocess is True:
                torepro[selection_name]=(effi_name,True)
                load_data=True
            else:
                torepro[selection_name]=(effi_name,False)

        if load_data:
            tot_data=self.Load_Data_Multiproc()

        for selection_name, val in self.selection.items():
            if torepro[selection_name][1]:
                print('Files do not exist - reprocess'+torepro[selection_name][0])
                if preselection is not None:
                    tot_data=self.Select_all(tot_data,preselection)
                self.efficiencies[selection_name]=self.Get_Efficiencies(tot_data,self.selection[selection_name],torepro[selection_name][0]) 
            else:
                #print('File already exist')
                self.efficiencies[selection_name]=np.load(torepro[selection_name][0])
                print(self.efficiencies[selection_name]['duration_z'])
        

    def Check_Zero_Efficiencies(self,tot_data):
        zlim=dict(zip([-2.0,0.0,2.0],[0.75,1.,1.]))
        print(tot_data.dtype)
        r=[]
        sel_data=self.Select_all(tot_data,self.selection)
        for x1_color in np.unique(tot_data[['X1','Color']]):
            for z in np.unique(tot_data['z']):
                idxa=(tot_data['z']==z)&(tot_data[['X1','Color']]==x1_color)
                idxb=(sel_data['z']==z)&(sel_data[['X1','Color']]==x1_color)
                        #print(fieldid,z,len(sel_data[idxb]),len(tot_data[idxa]))
                if len(sel_data[idxb])==0 and z < zlim[x1_color[0]]:
                    print('To reprocess',self.fieldid,self.season,z,x1_color)
                    r.append((self.fieldid,self.season,z,x1_color[0],x1_color[1]))
        if len(r) > 0:
            zeros=np.rec.fromrecords(r,names=['fieldid','season','z','X1','Color'])
            print(zeros)
            np.save('Toreprocess.npy',zeros)
        else:
            print('No zero found')

    def Select_all(self,tot_data,selec):
        
        valscp=tot_data.copy()
        for ic,sel in enumerate(selec):
            valscp=self.Select_Cut(valscp,sel[0],sel[1],sel[2])
        
        return valscp


    def Save_Selected(self):

        #Load the data:
        tot_data=self.Load_Data_Multiproc()
        #loop on selection
        for selec_name, vals in self.selection.items():
            sel_data=self.Select_all(tot_data,vals)
            save_dir='Sel_'+selec_name
            print('Number of sel :',selec_name,len(sel_data),len(tot_data))
            self.Check_Dir(save_dir)
            pkl_out=open(save_dir+'/Sel_'+self.fieldname+'_'+str(self.fieldid)+'_Season_'+str(self.season)+'.pkl','wb')
            pkl.dump(sel_data,pkl_out)
            pkl_out.close()
            """
            pkl_out=open(save_dir+'/Tot_'+self.fieldname+'_'+str(self.fieldid)+'_Season_'+str(self.season)+'.pkl','wb')
            pkl.dump(tot_data,pkl_out)
            
            pkl_out.close() 
            """
        #selection_name: Sel_


    def Get_Efficiencies(self,tot_data,selec,effi_name):
        
        r=[]
        
        self.Check_Dir(self.effi_dir)

        sel_data=self.Select_all(tot_data,selec)
        for x1_color in np.unique(tot_data[['X1','Color']]):
            for z in np.unique(tot_data['z']):
                idxa=(np.abs(tot_data['z']-z)<1.e-5)&(tot_data[['X1','Color']]==x1_color)
                idxb=(np.abs(sel_data['z']-z)<1.e-5)&(sel_data[['X1','Color']]==x1_color)
            
                effi=float(len(sel_data[idxb]))/float(len(tot_data[idxa]))
                err_effi=np.sqrt(1.-effi)*effi/np.sqrt(float(len(tot_data[idxa])))
                if self.preselection is not None:
                    sel_z=tot_data[idxa]
                    idc = sel_z['phase_min_5_all']<= -5.
                    idd = sel_z['phase_max_5_all'] >= 10.
                #print('hello',z,np.min(sel_z[idc]['DayMax']),np.max(sel_z[idd]['DayMax']))
                    duration_z=np.max(sel_z[idd]['DayMax'])-np.min(sel_z[idc]['DayMax'])
                else:
                    duration_z=tot_data[idxa]['duration'][0]
                print('hello',z,tot_data[idxa]['duration'][0],duration_z)


                r.append((self.fieldname,self.fieldid,self.season,z,x1_color[0],x1_color[1],effi,err_effi,tot_data[idxa]['duration'][0],duration_z))

                #print('helli',self.season,tot_data[idxa]['duration'])
            res= np.rec.fromrecords(r,names=['fieldname','fieldid','season','z','X1','Color','effi','err_effi','duration','duration_z'])
            np.save(effi_name,res)
            return res

    def Check_Dir(self,dirout):
        if not os.path.isdir(dirout) :
            os.makedirs(dirout)

    def Load_Indiv(self,files,j,output_q=None):


        tot_data=None
        """
        dirmeas=self.data_dir+'/'+val.thedir+'/'+val.fieldname+'/'+str(val.fieldid)+'/Season_'+str(val.season)
        what=dirmeas+'/z_'+str(val.z)+'/'+val.fieldname+'_'+str(val.fieldid)+'_'+str(val.z)+'_X1_'+str(val.X1)+'_C_'+str(val.Color)+'*.hdf5'
        files = glob.glob(what)
        """
        #print 'there man',len(files),dirmeas,what
            #tot_data=None
        for fi in files:
            f = h5py.File(fi,'r')
            print('try loading',fi,len(f.keys()))
            n_lc=0
            for i,keyb in enumerate(f.keys()):
                tab=Table.read(fi, path=keyb)
                #print('tab',tab)
                tab.add_column(Column([self.fieldid]*len(tab),name='fieldid'))
                tab.add_column(Column([self.season]*len(tab),name='season'))
                sigmu=[]
                if self.add_covmb:
                    for val in tab:
                        #print(val['salt2.Covmbmb'],val['salt2.CovX1X1'],val['salt2.CovColorColor'],self.alpha,self.beta)
                        sigmu_sq=val['salt2.Covmbmb']
                        sigmu_sq+=self.alpha**2*val['salt2.CovX1X1']+self.beta**2*val['salt2.CovColorColor']
                        sigmu_sq+=2.*self.alpha*val['salt2.CovX1mb']
                        sigmu_sq+=-2.*self.alpha*self.beta*val['salt2.CovColorX1']
                        sigmu_sq+=-2.*self.beta*val['salt2.CovColormb']
                        if sigmu_sq >= 0.:
                            sigmu.append(np.sqrt(sigmu_sq))
                        else:
                            sigmu.append(0)
                    tab.add_column(Column(sigmu,name='sigma_dist_modulus'))
                    """
                    for name in tab.colnames:
                        if 'salt2' in name:
                            print name
                    """

                n_lc+=len(tab)
                if tot_data is None:
                    tot_data=tab
                else:
                    tot_data=vstack([tot_data,tab])
        #print 'number of lc',n_lc
        #print('loading',key,len(tot_data))

        if output_q is not None:
            output_q.put({j : tot_data})
        else:
            return tot_data


    def Plot_Var(self,ax,varname,color='k'):

        #ax.hist(self.tot_data[varname],histtype='step',color=color)
        ax.plot(self.tot_data['z'],self.tot_data[varname],color=color,marker='.',linestyle='None')



    def Load_Data_Multiproc(self):


        dirmeas=self.data_dir+'/'+self.fieldname+'/'+str(self.fieldid)+'/Season_'+str(self.season)

        files=[]
        if len(self.zrange) >= 1:
            for z in self.zrange:
                what=dirmeas+'/z_'+str(z)+'/'+self.fieldname+'_'+str(self.fieldid)+'_'+str(z)+'_X1_'+str(self.x1_color[0])+'_C_'+str(self.x1_color[1])+'*.hdf5'
                files+= glob.glob(what)
        else:
            what=dirmeas+'/z_*/'+self.fieldname+'_'+str(self.fieldid)+'_*_X1_*_C_*.hdf5'
            files+= glob.glob(what) 
            

        n_per_batch=10
        
        inter=range(0,len(files),n_per_batch)
        inter=np.append(inter,len(files))
        
        restot=None
        time_begin=time.time()
        for jo in range(len(inter)-1):
            ida=inter[jo]
            idb=inter[jo+1]
            result_queue = multiprocessing.Queue()
       
            for j in range(ida,idb):
                
            #print('radec',j,radec[j])
                p=multiprocessing.Process(name='Subprocess-'+str(j),target=self.Load_Indiv,args=([files[j]],j,result_queue))
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

    def Select_Cut(self,data,var,comp,val):

        tot_data=None
        
        for vv in var:
            if tot_data is None:
                tot_data=data[vv]
            else:
                tot_data+=data[vv]


        idx = comp(tot_data,val)

        #print idx
        return data[idx]

    def Histo_ratio(self,sela,selb,varname,zmin,zmax,bin_z):

        range=[zmin,zmax]
        bin_width=bin_z
        #print 'alors pal',range,bin_width
        num_bins=int((range[1]-range[0])/bin_width)
        #range=[0.0,1.2]
        
        #hista, bin_edgesa = np.histogram(sela[varname],bins=num_bins,range=range,weights=sela['X1_weight']*sela['Color_weight'])
        #histb, bin_edgesb = np.histogram(selb[varname],bins=num_bins,range=range,weights=selb['X1_weight']*selb['Color_weight'])
        hista, bin_edgesa = np.histogram(sela[varname],bins=num_bins,range=range)
        histb, bin_edgesb = np.histogram(selb[varname],bins=num_bins,range=range)
        #print hista
        #print histb
       
        #bin_center = (bin_edgesa[:-1] + bin_edgesa[1:]) / 2

        
        ratio=[]
        ratio_err=[]
        norm=[]
        norm_err=[]
        
        for a,b in zip(histb,hista):
            if b==0:
                ratio.append(a)
                ratio_err.append(0)
            else:
                #print 'hello',a,b
                effi=float(a)/float(b)
                ratio.append(effi)
                ratio_err.append(np.sqrt(1.-effi)*effi/np.sqrt(float(b)))
            eff=float(a)/float(np.sum(hista))
            norm.append(eff)
            norm_err.append(np.sqrt(1.-eff)*eff/np.sqrt(float(np.sum(hista))))
    
        return bin_edgesa[:-1], ratio, ratio_err,norm,norm_err


    def Plot_Effi_new(self,ls='-',color='k',tot_label=[],take_label=False,ax=None,selection_name=''):

        tab=self.efficiencies[selection_name]
        #print 'hello',tab                                                                                                                                                                                                                                                
        idx =(tab['X1']==self.x1_color[0])&(tab['Color']==self.x1_color[1])
        sel=tab[idx]
        sel.sort(order='z')
            #print('hello',season,len(sel),sel['effi'],sel['z'])                                                                                       
            #print 'label',str(season+1)                                                                                                               
        if take_label:
            tot_label.append(ax.errorbar(sel['z'],sel['effi'],yerr=sel['err_effi'],label='Season '+str(self.season+1),ls=ls,color=color))
        else:
            ax.errorbar(sel['z'],sel['effi'],yerr=sel['err_effi'],ls=ls,color=color)


    def Plot_Effi(self,ax,tab,fieldid,season,x1_color,ls,color,tot_label=[],take_label=False):

        #print 'hello',tab
        for (x1,c) in x1_color:
            #print 'there pal',x1,c
            idx = (tab['fieldid']==fieldid)&(tab['season']==season)&(tab['X1']==x1)&(tab['Color']==c)
            sel=tab[idx]
            sel.sort(order='z')
            #print('hello',season,len(sel),sel['effi'],sel['z'])
            #print 'label',str(season+1)
            if take_label:
                tot_label.append(ax.errorbar(sel['z'],sel['effi'],yerr=sel['err_effi'],label='Season '+str(season+1),ls=ls,color=color))
            else:
                ax.errorbar(sel['z'],sel['effi'],yerr=sel['err_effi'],ls=ls,color=color)
        
    def Plot_Effi_vs(self,axc,tabref,tab,varname,title):

        min_bin=np.min(tabref[varname])
        max_bin=np.max(tabref[varname])
        delta_bin=3.
        

        bin_center, ratio, ratio_err,norm,norm_err= self.Histo_ratio(tabref,tab,varname,min_bin,max_bin,delta_bin)

        #figc, axc = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        #figc.suptitle(title)
        axc.errorbar(bin_center,ratio,yerr=ratio_err)


    def Plot_Efficiency(self,data,dict_sel,data_dict,varname,min_bin,max_bin,delta_bin,obs=None):
        tot_label=[]
        figc, axc = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        for key, vals in data.items():
            self.Plot_Efficiency_Indiv(axc,vals,dict_sel[key],varname,tot_label,key,data_dict[key].colorfig,data_dict[key].season,min_bin,max_bin,delta_bin,obs)
        
        plt.gcf().savefig('effi.png')

    def Plot_Efficiency_Indiv(self,axc,sela,selb,varname,tot_label,ll,color,season,min_bin,max_bin,delta_bin,label='',obs=None):

        filtercolors = {'u':'b', 'g':'c', 'r':'g', 'i':'y', 'z':'r', 'y':'m'}

        fontsize=12
        bin_center, ratio, ratio_err,norm,norm_err= self.Histo_ratio(sela,selb,varname,min_bin,max_bin,delta_bin)
        
        ll='Y'+str(season+1)
        #axc.errorbar(bin_center,ratio, yerr=ratio_err,marker=self.ms[season], mfc=self.color[season], mec=self.color[season], ms=8, linestyle='-',color=color,label=ll)
        #tot_label.append(axc.errorbar(bin_center,ratio, yerr=ratio_err,marker=self.ms[season], mfc=color, mec=color, ms=8, linestyle='-',color=color,label=label))
        tot_label.append(axc.errorbar(bin_center,ratio, yerr=ratio_err,marker=self.ms[season], linestyle='-',label=label))

    def Plot_Sims(self,dict_ana):
        
        def Plot_Indiv(axbb,tab_resu):
            fontsize=10

            for (j,vals) in enumerate(['DayMax','z','X1','Color']):
                if j==0 or j==2:
                    k=0
                else:
                    k=1

                if vals != 'X1' and vals != 'Color':
                    axbb[j/2][k].hist(tab_resu[vals],bins=10,histtype='step')
                else:
                    axbb[j/2][k].hist(tab_resu[vals],weights=tab_resu[vals+'_weight'],bins=10,histtype='step')
            
                axbb[j/2][k].set_xlabel(r''+vals+'_sim',{'fontsize': fontsize})
                axbb[j/2][k].set_ylabel(r'Number of Entries',{'fontsize': fontsize})
                print vals,np.mean(tab_resu[vals]),np.std(tab_resu[vals])

        figb, axb = plt.subplots(ncols=2, nrows=2, figsize=(10,9))
        for key in dict_ana.keys():
            Plot_Indiv(axb,dict_ana[key])
 
    def Plot_Npoints(self,dict_ana):
        
        def Plot_Indiv_vs(axbb,tab_resu):
            fontsize=10

            for (j,vals) in enumerate(['g','r','i','z']):
                if j==0 or j==2:
                    k=0
                else:
                    k=1

                axbb[j/2][k].scatter(tab_resu['z'],tab_resu['N_bef_'+vals]+tab_resu['N_aft_'+vals])
            
                axbb[j/2][k].set_xlabel(r'z',{'fontsize': fontsize})
                axbb[j/2][k].set_ylabel(r'Npoints( '+vals+' band)',{'fontsize': fontsize})
               

        figb, axb = plt.subplots(ncols=2, nrows=2, figsize=(10,9))
        for key in dict_ana.keys():
            Plot_Indiv_vs(axb,dict_ana[key])

    """
    def Plot_Var(self,data_dict,varx):

        #figc, axc = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        tot_tab=None
        data={}
        for key,vals in self.data.items():
            print 'hello',key
            
            valscp=vals.copy()
            for ic,sel in enumerate(self.mysel):
                valscp=self.Select_Cut(valscp,sel[0],sel[1],sel[2])

            if tot_tab is None:
                tot_tab=valscp
            else:
                print valscp.dtype,len(valscp)
                if len(valscp) > 0:
                    print valscp
                    tot_tab=np.concatenate((tot_tab,valscp))
    
            data[data_dict[key].fieldid]=tot_tab

        for key, vals in data.items():
            print key
    """
