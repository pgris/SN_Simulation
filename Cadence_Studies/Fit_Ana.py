import multiprocessing
from Fit_LC import *
import sncosmo
import cPickle as pkl

class Fit_Ana:
    def __init__(self, fname, telescope):
        self.Ana_LCs(fname,telescope)
        
    def Ana_LCs(self,fname,telescope):

        time_ref=time.time()

        transmission=telescope.throughputs
        for filtre in 'ugrizy':
            if telescope.airmass > 0:
                band=sncosmo.Bandpass(transmission.atmosphere[filtre].wavelen,transmission.atmosphere[filtre].sb, name='LSST::'+filtre,wave_unit=u.nm)
            else:
                band=sncosmo.Bandpass(transmission.system[filtre].wavelen,transmission.system[filtre].sb, name='LSST::'+filtre,wave_unit=u.nm) 
            sncosmo.registry.register(band, force=True)

        f = open(fname, "r")
        objs = []
        r=[]
        n_per_batch=10
        fit_lc=Fit_LC(telescope=telescope)

        while 1:
            try:
        #print 'trying'
                objs.append(pkl.load(f))
            except EOFError:
                break
        print(len(objs))
        for obj in objs:
        #print(len(obj),type(obj))
            inter=range(0,len(obj),n_per_batch)
            inter=np.append(inter,len(obj))
        
            for jo in range(len(inter)-1):
                ida=inter[jo]
                idb=inter[jo+1]
                result_queue = multiprocessing.Queue()
        
                for j in range(ida,idb):
                    p=multiprocessing.Process(name='Subprocess-'+str(j),target=self.Ana_LC,args=(obj[j],fit_lc,j,result_queue))
                    p.start()
 
                resultdict = {}

                for j in range(ida,idb):
                    resultdict.update(result_queue.get())

                for p in multiprocessing.active_children():
                    p.join()

                for j in range(ida,idb):
                    r+=resultdict[j][0]
                    names=resultdict[j][1]
            if len(r)%100 == 0:
                print('there',len(r),time.time()-time_ref)
            
        lc_ana=np.rec.fromrecords(r,names=names)
        outname=fname.replace('Light_Curves','SN_from_LC')
        fout = open(outname, "w")
        pkl.dump(lc_ana,fout)
        fout.close()
        #Plot(lc_ana)

    def Ana_LC(self,tab,fit_lc,j,output_q=None):

    
        idxc=tab['flux']/tab['fluxerr']>5.
        #fit this LC
        res,fitted_model,mbfit,fit_status=fit_lc(tab[idxc])

        ro, names=self.Fill_Infos(res,fitted_model,mbfit,fit_status)
        #print(res,fitted_model,mbfit,fit_status)

        nb,na=self.Bef_Aft(tab,tab.meta['DayMax'])
        nb5,na5=self.Bef_Aft(tab[idxc],tab.meta['DayMax'])
        phase_min,phase_max= self.phase(tab)
        phase_min_5,phase_max_5= self.phase(tab[idxc])

        #print('before',ro,names)
        ro+=[tab.meta['DayMax'],tab.meta['z'],self.snr(tab),self.snr(tab[idxc]),nb,na,nb5,na5,phase_min,phase_max,phase_min_5,phase_max_5]
        names+=['DayMax','z','snr_all','snr_5_all','Nbef_all','Naft_all','Nbef_5_all','Naft_5_all','phase_min_all','phase_max_all','phase_min_5_all','phase_max_5_all']

        for band in 'grizy':
            #ro=[]
            idx = tab['band']=='LSST::'+band
            idxb= (tab['band']=='LSST::'+band)&(tab['flux']/tab['fluxerr']>5.)
            #print(band,snr(tab[idx]),snr(tab[idxb]),Bef_Aft(tab[idx],tab.meta['DayMax']))
            nb,na=self.Bef_Aft(tab[idx],tab.meta['DayMax'])
            nb5,na5=self.Bef_Aft(tab[idxb],tab.meta['DayMax'])
            phase_min,phase_max= self.phase(tab[idx])
            phase_min_5,phase_max_5= self.phase(tab[idxb])
            #r.append((tab.meta['DayMax'],tab.meta['z'],band,snr(tab[idx]),snr(tab[idxb]),nb,na,nb5,na5,phase_min,phase_max,phase_min_5,phase_max_5))
            ro+=[self.snr(tab[idx]),self.snr(tab[idxb]),nb,na,nb5,na5,phase_min,phase_max,phase_min_5,phase_max_5]
            names+=['snr_'+band,'snr_5_'+band,'Nbef_'+band,'Naft_'+band,'Nbef_5_'+band,'Naft_5_'+band,'phase_min_'+band,'phase_max_'+band,'phase_min_5_'+band,'phase_max_5_'+band]
    
            #r.append((tab.meta['DayMax'],tab.meta['z'],'all',snr(tab),snr(tab[idxc]),nb,na,nb5,na5,phase_min,phase_max,phase_min_5,phase_max_5))

        rsb, nsb= self.Var_Science_Book(tab,tab.meta['DayMax'],tab.meta['z'])

        ro+=rsb
        names+=nsb

        r=[]
        r.append(tuple(ro))

        if output_q is not None:
            output_q.put({j : (r,names)})
        else:
            return r

    def Fill_Infos(self,res,fitted_model,mbfit,fit_status):

        r=[]
        names=[]
        vars_salt=['x1','x0','c']
        corresp=dict(zip(vars_salt,['X1','X0','Color']))
        corresp['t0']='t0'
        corresp['z']='z'
        r+=[fit_status]
        names+=['fit_status']
        if fit_status == 'fitok':
            r+=[mbfit]
            names+=['mbfit']

            corr={}
            for i,pal in enumerate(res['vparam_names']):
                corr[pal]=i
            for i,par in enumerate(fitted_model.param_names):
                r+=[fitted_model.parameters[i]]
                names+=['salt2.'+corresp[par]]
            if res['covariance'] is not None:
                for k,var in enumerate(vars_salt):
                    for varb in vars_salt[k:]:
                        r+=[res['covariance'][corr[var]][corr[varb]]]
                        names+=['salt2.Cov'+corresp[var]+corresp[varb]]
                r+=[res.errors['t0']**2]
                names+=['salt2.Covt0t0']


        else: 
            r+=[-999.]
            names+=['mbfit']
            thenames=['salt2.z', 'salt2.t0', 'salt2.X0', 'salt2.X1', 'salt2.Color', 'salt2.CovX1X1', 'salt2.CovX1X0', 'salt2.CovX1Color', 'salt2.CovX0X0', 'salt2.CovX0Color', 'salt2.CovColorColor', 'salt2.Covt0t0']
            r+=[-999.]*len(thenames)
            names+=thenames
       
        return r, names


    def Var_Science_Book(self,lc,T0,z):
       
        r=[]
        names=[]
 
        lc.sort('time')
        phases=(lc['time']-T0)/(1.+z)
        #print('phases',phases)
        idx = phases < -5.
        idxb = phases > 30.
        r+=[len(phases[idx]),len(phases[idxb])]
        names+=['N_Phase_m5','N_Phase_p30']
    

        diffa,seldiffa=self.Get_Diff(lc,phases,-20.,60.)
    

        r+=[len(seldiffa)-1]
        names+=['N_nights_m20_p_p30']

        diffb,seldiffb=self.Get_Diff(lc,phases,-5.,30.)
        nearpeak=-1
        if len(diffb) > 0:
            nearpeak=np.max(diffb)/(1.+z)
        r+=[nearpeak]
        names+=['Near_peak_gap']

        snr_refs=[10.,15.,20.]
        for band in 'grizy':
            idf = lc['band']=='LSST::'+band
            sel_idf= lc[idf]
            if len(sel_idf) > 0:
                #print('band',band,len(sel_idf))
                snr_max=np.max(sel_idf['flux']/sel_idf['fluxerr'])
                for valref in snr_refs:
                    nres=0
                    if snr_max > valref:
                        nres=1
                    r+=[nres]
                    names+=['N_'+band+'_snrmax_'+str(int(valref))]
                    #self.dict_quality['N_'+band+'_snrmax_'+str(int(valref))]=nres
        
            else:
                for valref in snr_refs:
                    r+=[0.]
                    names+=['N_'+band+'_snrmax_'+str(int(valref))]
       

        return r,names

    def Get_Diff(self,lc,phases,phase_min,phase_max):

        idxc=(phases > phase_min)&(phases < phase_max)
    
        sel=lc[idxc]
    
        diff_time=1.
    
        diff=[io-jo for jo,io in zip(sel['time'][:-1], sel['time'][1:])]

        seldiff=[i+1 for i in range(len(diff)) if diff[i]>= diff_time]
        seldiff=[0]+seldiff+[len(sel)]

        return diff,seldiff

    def phase(self,tab):
        if len(tab) > 0:
            phases=(tab['time']-tab.meta['DayMax'])/(1.+tab.meta['z'])
            return np.min(phases),np.max(phases)
        else:
            return -999.,-999. 

    def snr(self,sel):
        if len(sel) > 0:
            return np.sum(sel['flux'])/np.sqrt(np.sum(sel['fluxerr']*sel['fluxerr']))
        else:
            return 0.0

    def Bef_Aft(self,lc, T0):
        if len(lc) > 0:
            diff=lc['time']-T0
            idxa=diff <= 0
            idxb=diff > 0
            return len(lc[idxa]),len(lc[idxb])
        else:
            return 0.0,0.0
