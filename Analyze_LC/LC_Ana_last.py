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
import operator
from scipy.interpolate import griddata
from Load_X1_C import Load_X1_Color
from matplotlib import colors as mcolors
import multiprocessing
import time
import os

class LC_Ana:
    def __init__(self,fieldname,fieldid,season,zrange,x1_color,main_dir='/sps/lsst/data/dev/pgris',data_dir='../../Fitted_Light_Curves_snsim',selection_name='',selection=[],simu_name='',Force_reprocess=False):

        self.fieldname=fieldname
        self.fieldid=fieldid
        self.season=season
        self.zrange=zrange
        self.x1_color=x1_color
        self.data_dir=main_dir+'/'+data_dir
        self.selection=selection[selection_name]
        #self.selection_name=selection_name

        self.effi_dir='Effi'+'_'+simu_name+'_'+selection_name
        self.Check_Dir(self.effi_dir)
        self.effi_name=self.effi_dir+'/Effi_'+self.fieldname+'_'+str(self.fieldid)+'_'+str(self.season)+'_'+str(self.x1_color[0])+'_'+str(self.x1_color[1])+'.npy'

        if not os.path.isfile(self.effi_name) or Force_reprocess is True:
            tot_data=self.Load_Data_Multiproc()
            self.efficiencies=self.Get_Efficiencies(tot_data) 
        else:
            self.efficiencies=np.load(self.effi_name)
            
        self.nsn_tot=0.
        self.err_tot=0.
        self.nsn_theo=0
        self.err_theo=0
        self.res_nsn=[]

        #self.Check_Zero_Efficiencies(tot_data)
        #self.efficiencies=self.Get_Efficiencies(tot_data)
        
        #self.Plot_Effi_new()

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

    """
    def Get_Selected(self,data_dict,selec):
        r=[]
        for key,vals in self.data.items():
            #print 'hello',key,len(vals),vals.dtype
            if vals is not None:
                valscp=vals.copy()
                for ic,sel in enumerate(selec):
                    valscp=self.Select_Cut(valscp,sel[0],sel[1],sel[2])
                    #print('after selection',len(valscp),len(vals),data_dict[key].fieldname,data_dict[key].fieldid,data_dict[key].z,data_dict[key].season)
                for val in valscp:    
                    r.append((data_dict[key].fieldname,data_dict[key].fieldid,data_dict[key].z,data_dict[key].season,data_dict[key].X1,data_dict[key].Color,np.sqrt(val['salt2.CovColorColor']),np.sqrt(val['salt2.CovX1X1']),vals['duration'][0],val['DayMax'],val['N_nights_m20_p_p30']))
    
        print('hello a',len(r))
        return np.rec.fromrecords(r,names=['fieldname','fieldid','z','season','X1','Color','sigma_color','sigma_X1','duration','DayMax','N_nights_m20_p_p30'])

    def Plot(self,data_dict,selec,z,x1_color=[(0.0,0.0)],colors=['ko','r*']):

        color=dict(zip([(-2.0,0.2),(2.0,-0.2),(0.0,0.0)],(['g','r','b'])))
        #color=dict(zip([0.1,0.6],(['g','r'])))
        figa, axa = plt.subplots(ncols=1, nrows=1)
        
        restot=[]
        for sel in selec:
            restot.append(self.Get_Selected(data_dict,sel))

        
        for i,res in enumerate(restot):
            #print(z,np.unique(res['z']),np.unique(res[['X1','Color']]))
            for (x1,c) in x1_color:
                #print('rr',x1,c)
                idd = (res['X1']==x1)&(res['Color']==c)&(np.abs(res['z']-z)<1.e-8)
                sel=res[idd]
                axa.plot(sel['DayMax'],sel['sigma_X1'],colors[i])
                #axa.plot(sel['DayMax'],sel['N_nights_m20_p_p30'],color[z]+'o')
            #axa.plot(vals['DayMax'],np.sqrt(vals['salt2.CovColorColor']),'k*')
        
        plt.show()

    """
    def Get_Efficiencies(self,tot_data):
        
        r=[]
        
        self.Check_Dir(self.effi_dir)

        sel_data=self.Select_all(tot_data,self.selection)
        for x1_color in np.unique(tot_data[['X1','Color']]):
            for z in np.unique(tot_data['z']):
                idxa=(tot_data['z']==z)&(tot_data[['X1','Color']]==x1_color)
                idxb=(sel_data['z']==z)&(sel_data[['X1','Color']]==x1_color)
            
                effi=float(len(sel_data[idxb]))/float(len(tot_data[idxa]))
                err_effi=np.sqrt(1.-effi)*effi/np.sqrt(float(len(tot_data[idxa])))
                r.append((self.fieldname,self.fieldid,self.season,z,x1_color[0],x1_color[1],effi,err_effi,tot_data[idxa]['duration'][0]))

                #print 'helli',len(r)
            res= np.rec.fromrecords(r,names=['fieldname','fieldid','season','z','X1','Color','effi','err_effi','duration'])
            np.save(self.effi_name,res)
            return res

    """
    def Plot_Efficiencies_Multiple(self,fieldid,data_dict,selec,x1_color=(0.0,0.0),eff_zero=[]):

        fig, axes = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        fig.suptitle('Field '+str(fieldid)+' - (x1,c) ='+str(x1_color))
        title='Field - '+str(fieldid)
        fontsize=12
        myls=['-','--']
        
        colors=dict(zip([i for i in range(10)],['blue','orange','green','red','black','brown','pink','gray','olive','cyan']))
        #colors=dict(zip([i for i in range(10)],mycolors))

        io=-1
        tot_label=[]
        zlim=dict(zip([-2.0,0.0,2.0],[0.75,1.1,1.1]))
        for key, vals in selec.items():
            io+=1
            res=self.Get_Efficiencies(data_dict,vals)
            idx = res['fieldid']==fieldid
            effi_f=res[idx]
            for season in np.unique(effi_f['season']):
                ijk=(effi_f['effi']<1.e-5)&(effi_f['season']==season)&(effi_f['z']<zlim[x1_color[0]])&(effi_f['X1']==x1_color[0])&(effi_f['Color']==x1_color[1])
                test_zero=effi_f[ijk]
                if len(test_zero) > 0:
                    print('Efficiency equal 0',season,test_zero['fieldid'],test_zero['z'],test_zero['effi'],x1_color)
                    for pval in test_zero:
                        eff_zero.append((pval['fieldid'],pval['season'],round(pval['z'],4),x1_color[0],x1_color[1]))
                if io ==0:
                    self.Plot_Effi(axes,effi_f,fieldid,season,x1_color=[x1_color],ls=myls[io],color=colors[season],tot_label=tot_label,take_label=True)
                else:
                    self.Plot_Effi(axes,effi_f,fieldid,season,x1_color=[x1_color],ls=myls[io],color=colors[season])
        labs = [l.get_label() for l in tot_label]
        axes.legend(tot_label, labs, ncol=5,loc='best',prop={'size':fontsize},frameon=False)
        axes.set_xlabel('z')
        axes.set_ylabel('Efficiency')
        #plt.show()
    """
    def Estimate_Efficiencies(self,data_dict,selec,selname,interp=False,x1_color=[(0.,0.)]):
        
        res=self.Get_Efficiencies(data_dict,selec)
        
        #print res

        #effi=self.Interpolate_Effi(res)
        effi=res
        x1_color=x1_color
        if interp:
           effi=self.Interpolate_Effi(res) 
           x1_color=[(-999.0,-999.0)]
        #x1_color=[(0.,0.)]
        #print effi[['fieldid','season']]
        fontsize=12
        cumul=True 
        myls=['-','--']
        colors=dict(zip([i for i in range(10)],['k','k','r','r','b','b','g','g','m','m']))
        for fieldid in np.unique(effi['fieldid']):
            fig, axes = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
            figb, axesb = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
            fig.suptitle('Field '+str(fieldid))
            title='Field - '+str(fieldid)
            self.nsn_tot=0.
            self.err_tot=0.
            self.nsn_theo=0
            self.err_theo=0
            self.res_nsn=[]
            idx = effi['fieldid']==fieldid
            effi_f=effi[idx]
            #print 'hello',effi_f.dtype,effi_f
            for season in np.unique(effi_f['season']):
                self.Plot_Effi(axes,effi_f,fieldid,season,x1_color,ls=myls[season%2],color=colors[season])
                self.Plot_N_SN(axesb,effi_f,effi_f['fieldname'][0],fieldid,season,cumul=cumul)
                
            axes.legend(loc='upper right',prop={'size':fontsize})
            axes.set_xlabel('z',{'fontsize': fontsize})
            axes.set_ylabel('Efficiency',{'fontsize': fontsize})

            axesb.legend(loc='upper left',prop={'size':fontsize})
            axesb.set_xlabel('z',{'fontsize': fontsize}) 
            if cumul:
                axesb.set_ylabel('Number of SN Ia < z')
            else:
                axesb.set_ylabel('Number of SN Ia per z bin')
            
            nsntot_str=str(int(self.nsn_tot))+'$\pm$'+str(int(np.power(self.err_tot,0.5)))
            nsntheo_str=str(int(self.nsn_theo))+'$\pm$'+str(int(np.power(self.err_theo,0.5)))
            title+=' - N$_{SN Ia}$ = '+nsntot_str +' / '+nsntheo_str
            axesb.set_title(title)
            fig.savefig('effi_'+str(fieldid)+'.png')
            figb.savefig('nsn_'+str(fieldid)+'.png')
            if len(self.res_nsn) > 0:
                rate_nsn=np.rec.fromrecords(self.res_nsn,names=['fieldname','fieldid','season','n_sn_detected','err_detected','n_sn_expected','err_expected'])
                #print 'allo',effi_f['fieldname']
                self.Check_Dir(outname)
                pkl_out=open(outname+'/N_SN_'+effi_f['fieldname'][0]+'_'+str(fieldid)+'.pkl','wb')
                pkl.dump(rate_nsn,pkl_out)
                pkl_out.close()
            
        #plt.show()

    def Check_Dir(self,dirout):
        if not os.path.isdir(dirout) :
            os.makedirs(dirout)

    def Interpolate_Effi(self,tab):

        restot=None
        
        for fieldid in np.unique(tab['fieldid']):
            idx = tab['fieldid']==fieldid
            tab_f=tab[idx]
            for season in np.unique(tab_f['season']):
                res=self.Interpolate_elem(tab_f,tab_f['fieldname'][0],fieldid,season)
                print 'interp',fieldid,season,len(res)
                if restot is None:
                    restot=res
                else:
                    restot=np.concatenate((restot,res))
        return restot

    def Interpolate_elem(self,tab,fieldname,fieldid,season):

        grid_x1_c=Load_X1_Color(-999.0,-999.0,1).tab

        """
        print 'ici',grid_x1_c.keys()
        for key in grid_x1_c.keys():
            print key, len(grid_x1_c[key])
        """
        idxa= (tab['fieldid']==fieldid)&(tab['season']==season)
        sel=tab[idxa]
        #print 'hhh',len(sel),fieldid,season
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
                """
                print 'hello',len(selb),z,grid
                print len(selb[['X1','Color']]),len(selb['Effi']), len(grid['x1']),len(grid['c'])
                print np.asarray(selb[['X1','Color']])
                print np.asarray(grid[['x1']]).shape
                
                gridx, gridy = np.meshgrid(grid['x1'],grid['c'])
                print (gridx,gridy)
                """
            #gridx, gridy = np.meshgrid(grid['x1'],grid['c'])
                grid_res=griddata(np.asarray(selb[['X1','Color']],dtype=np.float),np.asarray(selb['effi'],dtype=np.float),np.asarray(grid[['x1','c']], dtype=np.float), method='nearest')
                grid_res_err=griddata(np.asarray(selb[['X1','Color']],dtype=np.float),np.asarray(selb['err_effi'],dtype=np.float),np.asarray(grid[['x1','c']], dtype=np.float), method='nearest')
            #print grid_res
        
                effi_comb=np.sum([grid['weight_x1'][i]*grid['weight_c'][i]*grid_res[i] for i in range(len(grid_res))])
                effi_comb_err=np.sum([grid['weight_x1'][i]*grid['weight_c'][i]*grid_res_err[i] for i in range(len(grid_res_err))])
            #print effi_comb,effi_comb_err
                r.append((fieldname,fieldid,season,z,effi_comb,effi_comb_err,duration,-999.0,-999.0))
            else:
                r.append((fieldname,fieldid,season,z,0.,0.,duration,-999.0,-999.0))


            """
            #plt.imshow(grid_res.T, extent=np.meshgrid(grid['x1'],grid['c']), origin='lower')
            gridx, gridy = np.meshgrid(grid['x1'],grid['c']) 
            #plt.contourf(grid['x1'],grid['c'],grid_res,cmap=plt.cm.jet,s=5)
            plt.scatter(grid['x1'],grid['c'],c=grid_res,cmap=plt.cm.jet,s=100,marker='s')
            plt.colorbar() # draw colorbar
            plt.scatter(selb['X1'],selb['Color'],marker='o',c='b',s=5)
            plt.title('Nearest')
            plt.show()
            """
        return np.rec.fromrecords(r,names=['fieldname','fieldid','season','z','effi','err_effi','duration','X1','Color'])


    def Dump_In_File(self,filename,what,iseq):
        if not os.path.isfile(filename):
            what.write(filename, path='fit_lc_all_'+str(iseq), compression=True)
        else:
            what.write(filename, path='fit_lc_all_'+str(iseq), append=True, compression=True)

    def Load_Data(self,data_dict={}):

        print('loading data',len(data_dict))
        for key,val in data_dict.items():
            dirmeas=self.data_dir+'/'+val.thedir+'/'+val.fieldname+'/'+str(val.fieldid)+'/Season_'+str(val.season)
            what=dirmeas+'/z_'+str(val.z)+'/'+val.fieldname+'_'+str(val.fieldid)+'_'+str(val.z)+'_X1_'+str(val.X1)+'_C_'+str(val.Color)+'*.hdf5'
            files = glob.glob(what)
            print 'there man',len(files),dirmeas,what
            tot_data=None
            for fi in files:
                f = h5py.File(fi,'r')
                print('try loading',key,fi,len(f.keys()))
                n_lc=0
                for i,keyb in enumerate(f.keys()):
                    tab=Table.read(fi, path=keyb)
                    n_lc+=len(tab)
                    if tot_data is None:
                        tot_data=tab
                    else:
                        tot_data=vstack([tot_data,tab])
                print 'number of lc',n_lc
            print('loading',key,len(tot_data))
            self.data[key]=tot_data

    def Load_Data_new(self,data_dict={}):

        print('loading data',len(data_dict))
        tot_data=None
        for key,val in data_dict.items():
            dirmeas=self.data_dir+'/'+val.thedir+'/'+val.fieldname+'/'+str(val.fieldid)+'/Season_'+str(val.season)
            what=dirmeas+'/z_'+str(val.z)+'/'+val.fieldname+'_'+str(val.fieldid)+'_'+str(val.z)+'_X1_'+str(val.X1)+'_C_'+str(val.Color)+'*.hdf5'
            files = glob.glob(what)
            print 'there man',len(files),dirmeas,what
            #tot_data=None
            for fi in files:
                f = h5py.File(fi,'r')
                print('try loading',key,fi,len(f.keys()))
                n_lc=0
                for i,keyb in enumerate(f.keys()):
                    tab=Table.read(fi, path=keyb)
                    tab.add_column(Column([val.fieldid]*len(tab),name='fieldid'))
                    tab.add_column(Column([val.season]*len(tab),name='season'))
                    n_lc+=len(tab)
                    if tot_data is None:
                        tot_data=tab
                    else:
                        tot_data=vstack([tot_data,tab])
                print 'number of lc',n_lc
            print('loading',key,len(tot_data))
            #self.data[key]=tot_data
        return tot_data

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

    def Load_Data_Multiproc(self):


        dirmeas=self.data_dir+'/'+self.fieldname+'/'+str(self.fieldid)+'/Season_'+str(self.season)

        files=[]
        for z in self.zrange:
            what=dirmeas+'/z_'+str(z)+'/'+self.fieldname+'_'+str(self.fieldid)+'_'+str(z)+'_X1_'+str(self.x1_color[0])+'_C_'+str(self.x1_color[1])+'*.hdf5'
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

    def Load_Data_old(self, data_dict={}):

        print 'Data dir',self.data_dir
        for key,val in data_dict.items():
            dirmeas=self.data_dir+'/'+val.thedir+'/'+val.fieldname+'/'+str(val.fieldid)+'/Season_'+str(val.season)
            sum_file=dirmeas+'/'+val.fieldname+'_'+str(val.fieldid)+'_X1_'+str(val.X1)+'_C_'+str(val.Color)+'_all.hdf5'
            nfich=-1
            print 'sum file',sum_file
            if not os.path.exists(sum_file):
                tot_data=None
                inum=-1
                for nz,z in enumerate(val.z):
                    zstr=str(z)
                    
                    if z < 0:
                        zstr='*'
                    if val.X1 > -90. and val.Color > -90:
                        what=dirmeas+'/z_'+zstr+'/'+val.fieldname+'_'+str(val.fieldid)+'_'+zstr+'_X1_'+str(val.X1)+'_C_'+str(val.Color)+'*.hdf5'
                    else:
                        what=dirmeas+'/z_'+zstr+'/'+val.fieldname+'_'+str(val.fieldid)+'_'+zstr+'_X1_*_C_*.hdf5'
                    files = glob.glob(what)
                    print 'there man',len(files),dirmeas,what
                    for fi in files:
                        nfich+=1
                        f = h5py.File(fi,'r')
                        print 'loading',fi,len(f.keys())
                        n_lc=0
                        for i,keyb in enumerate(f.keys()):
                            tab=Table.read(fi, path=keyb)
                            n_lc+=len(tab)
                            if tot_data is None:
                                tot_data=tab
                            else:
                                tot_data=vstack([tot_data,tab])
                        print 'number of lc',n_lc
                        if nfich > 1 and nfich%3 == 0:
                            inum+=1
                            print 'Dumping',nfich
                            self.Dump_In_File(sum_file,tot_data,inum)
                            tot_data=None
                    if tot_data is not None:
                        inum+=1
                        self.Dump_In_File(sum_file,tot_data,inum)
                        tot_data=None
                   

            f = h5py.File(sum_file,'r')
            for i,keyb in enumerate(f.keys()):
                tab=Table.read(f, path=keyb)
                if not key in self.data.keys():
                    self.data[key]=tab
                else:
                    self.data[key]=vstack([self.data[key],tab])
            print 'nlc',len(self.data[key])

                #pkl_out = open(sum_file,'wb')
                
                #pkl.dump(self.data[key], pkl_out)
                
                #pkl_out.close()

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

    def Select_Data(self,data):

        sel=data.copy()
        logical_test={}
    
        #rint 'before',len(sel)
        """
        for band in 'grizy':
            logical_test[band]=np.logical_and(sel['N_bef_'+band]>=1,sel['N_aft_'+band]>=1)

        logical_and_g_r=np.logical_and(logical_test['g'],logical_test['r'])
        logical_and_r_i=np.logical_and(logical_test['r'],logical_test['i'])
        logical_and_i_z=np.logical_and(logical_test['i'],logical_test['z'])
        logical_and_z_y=np.logical_and(logical_test['z'],logical_test['y'])
        """            
        #print 'logical',len(sel)
        print sel.dtype
        idx=(sel['N_bef_all']>=4)&(sel['N_aft_all']>=10)
        idx&=(sel['status']=='go_fit')&(sel['fit_status']=='fit_ok')
        #idx&=(sel['phase_first']<=-5)&(sel['phase_last']>=20)
        
        """
        for val in sel:
            print val['phase_first'],(val['DayMax']-val['salt2.T0'])/(1.+val['z'])
        """
        #sel=sel[idx]
        
        #print 'all',len(sel)
        return sel[idx]

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


    def Plot_Effi_new(self,ls='-',color='k',tot_label=[],take_label=False,ax=None):

        tab=self.efficiencies
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


    def Plot_N_SN(self,axb,tab_tot,fieldname,fieldid,season,cumul=False):

        color='black'

        idx=(tab_tot['fieldid']==fieldid)&(tab_tot['season']==season)
        tab=tab_tot[idx].copy()
        tab.sort(order='z')
        #print 'ooo',tab['z'],type(tab)
        effi = interpolate.interp1d(tab['z'],tab['effi'])
        ratio_err=tab['err_effi']

        rate_name='Perrett'
        sn_rate=SN_Rate(rate=rate_name,duration=(tab['duration'][0]+0.0)/365.25)
        print 'attention',tab['duration']

        #zz,rate,err_rate,nsn,err_nsn=sn_rate(self.zmin,self.zmax,self.bin_z)                                                                                   
        zz,rate,err_rate,nsn,err_nsn=sn_rate(bins=tab['z'])

        print 'Nsn',fieldid,season,np.sum(nsn),zz,tab['z'],nsn,err_nsn,np.power(np.sum(np.power(err_nsn,2.)),0.5)
    
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
        self.res_nsn.append((fieldname,fieldid,int(season),np.sum(combi),err_N_sn,np.sum(nsn_season(zz)),np.sum(np.power(err_nsn_season(zz),2.))))

        nsn_str= str(int(np.sum(nsn_season(zz))))+'$\pm$'+str(int(np.power(np.sum(np.power(err_nsn_season(zz),2.)),0.5)))
        if season < 9:
            ll='Y'+str(season+1)+'   - N$_{SN Ia}$ = '+str(int(N_sn))+' $\pm$ '+str(int(err_N_sn))+' / '+nsn_str
        else:
            ll='Y'+str(season+1)+' - N$_{SN Ia}$ = '+str(int(N_sn))+' $\pm$ '+str(int(err_N_sn))+ ' / '+nsn_str

        if not cumul:
            axb.errorbar(zz,nsn_season(zz),yerr=err_nsn_season(zz),marker=self.ms[season], mfc=self.color[season], mec=self.color[season], ms=8, linestyle='--',color='black')
            axb.errorbar(zz,combi,yerr=yerr_combi,marker=self.ms[season], mfc=self.color[season], mec=self.color[season], ms=8, linestyle='-',color=color,label=ll)
        else:
            axb.errorbar(zz,np.cumsum(nsn_season(zz)),marker=self.ms[season], mfc=self.color[season], mec=self.color[season], ms=8, linestyle='--',color='black')
            axb.errorbar(zz,np.cumsum(combi),marker=self.ms[season], mfc=self.color[season], mec=self.color[season], ms=8, linestyle='-',color=color,label=ll)





    def Plot_N_SN_old(self,data,dict_sel,data_dict,varname,min_bin,max_bin,delta_bin,cumul):

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

        nsn_str= str(int(np.sum(nsn_season(zz))))+'$\pm$'+str(int(np.power(np.sum(np.power(err_nsn_season(zz),2.)),0.5)))
        if season < 9:
            ll='Y'+str(season+1)+'   - N$_{SN Ia}$ = '+str(int(N_sn))+' $\pm$ '+str(int(err_N_sn))+' / '+nsn_str
        else:
            ll='Y'+str(season+1)+' - N$_{SN Ia}$ = '+str(int(N_sn))+' $\pm$ '+str(int(err_N_sn))+ ' / '+nsn_str

        if not cumul:
            axb.errorbar(zz,nsn_season(zz),yerr=err_nsn_season(zz),marker=self.ms[season], mfc=self.color[season], mec=self.color[season], ms=8, linestyle='--',color='black')
            axb.errorbar(zz,combi,yerr=yerr_combi,marker=self.ms[season], mfc=self.color[season], mec=self.color[season], ms=8, linestyle='-',color=color,label=ll)
        else:
            axb.errorbar(zz,np.cumsum(nsn_season(zz)),marker=self.ms[season], mfc=self.color[season], mec=self.color[season], ms=8, linestyle='--',color='black')
            axb.errorbar(zz,np.cumsum(combi),marker=self.ms[season], mfc=self.color[season], mec=self.color[season], ms=8, linestyle='-',color=color,label=ll)
 
        

        axb.set_xlabel('z')
        if cumul:
            axb.set_ylabel('Number of SN Ia < z')
        else:
           axb.set_ylabel('Number of SN Ia per z bin') 
        axb.legend(loc='best',prop={'size':12})
        #axb.set_xlim(self.zmin,np.max(bin_center)+0.01)


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

    def Plot_debug(self,dict_ana):
        
        def Plot_Indiv_vs(axbb,tab_resu):
            fontsize=10

            idf = tab_resu['z']==0.2
            idfb = tab_resu['z']==0.225

            sela=tab_resu[idf]
            selb=tab_resu[idfb]

            for i,vals in enumerate(['X1','Color']):
                axbb[i].hist(sela[vals],weights=sela['X1_weight']*sela['Color_weight'],bins=10,histtype='step')
                axbb[i].hist(selb[vals],weights=selb['X1_weight']*selb['Color_weight'],bins=10,histtype='step')
                axbb[i].set_xlabel(r''+vals,{'fontsize': fontsize})
            #axbb[j/2][k].set_ylabel(r'Npoints( '+vals+' band)',{'fontsize': fontsize})
        
        def Plot_Indiv_vsb(axbb,tab_resu):
            fontsize=10

            idf = tab_resu['z']==0.2
            idfb = tab_resu['z']==0.225

            sela=tab_resu[idf]
            selb=tab_resu[idfb]

            print sela['DayMax'],len(sela['DayMax']),np.min(sela['DayMax']),np.max(sela['DayMax']),np.min(selb['DayMax']),np.max(selb['DayMax'])

            for val in sela['DayMax']:
                if val not in selb['DayMax']:
                    print 'problem here',val

            nbins_daymax=(np.max(selb['DayMax'])-np.min(selb['DayMax']))/0.2
            #bins=dict(zip(['DayMax','phase_last'],[nbins_daymax,10]))
            bins=dict(zip(['DayMax','phase_last'],[10,10]))
            for i,vals in enumerate(['DayMax','phase_last']):
                axbb[i].hist(sela[vals],bins=bins[vals],histtype='step')
                axbb[i].hist(selb[vals],bins=bins[vals],histtype='step')
                axbb[i].set_xlabel(r''+vals,{'fontsize': fontsize})
                    


        figb, axb = plt.subplots(ncols=2, nrows=1, figsize=(10,9))
        for key in dict_ana.keys():
            Plot_Indiv_vs(axb,dict_ana[key])
 
        figc, axc = plt.subplots(ncols=2, nrows=1, figsize=(10,9))
        for key in dict_ana.keys():
            Plot_Indiv_vsb(axc,dict_ana[key])

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
