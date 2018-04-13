from Fit_LC import *
import cPickle as pkl
import matplotlib.pyplot as plt
from astropy.table import vstack,Table,Column
import h5py
from astropy.table import Table


class Fit_Single_LC:
    def __init__(self, lc,telescope,inum,output_q=None):
        #time_begin=time.time()
        self.lc=lc
        self.telescope=telescope
        #print lc.dtype
        self.bands_rest = 'grizy'
        self.dict_quality={}
        self.outdict={}
        

        if lc is not None:

            idxc=lc['flux']/lc['fluxerr']>5.
            myfit=Fit_LC(z=lc.meta['z'],telescope=self.telescope,Plot=False)
        
            res,fitted_model,mbfit,fit_status=myfit(lc)
            ro, names=self.Fill_Infos(res,fitted_model,mbfit,fit_status)
            nb,na=self.Bef_Aft(lc,lc.meta['DayMax'])
            nb5,na5=self.Bef_Aft(lc[idxc],lc.meta['DayMax'])
            phase_min,phase_max= self.phase(lc)
            phase_min_5,phase_max_5= self.phase(lc[idxc])
            
            #print('before',ro,names)
            ro+=[self.snr(lc),self.snr(lc[idxc]),nb,na,nb5,na5,phase_min,phase_max,phase_min_5,phase_max_5]
            names+=['snr_all','snr_5_all','Nbef_all','Naft_all','Nbef_5_all','Naft_5_all','phase_min_all','phase_max_all','phase_min_5_all','phase_max_5_all']
            for band in 'grizy':
                #ro=[]
                idx = lc['band']=='LSST::'+band
                idxb= (lc['band']=='LSST::'+band)&(lc['flux']/lc['fluxerr']>5.)
                nb,na=self.Bef_Aft(lc[idx],lc.meta['DayMax'])
                nb5,na5=self.Bef_Aft(lc[idxb],lc.meta['DayMax'])
                phase_min,phase_max= self.phase(lc[idx])
                phase_min_5,phase_max_5= self.phase(lc[idxb])
                #r.append((lc.meta['DayMax'],lc.meta['z'],band,snr(lc[idx]),snr(lc[idxb]),nb,na,nb5,na5,phase_min,phase_max,phase_min_5,phase_max_5))
                ro+=[self.snr(lc[idx]),self.snr(lc[idxb]),nb,na,nb5,na5,phase_min,phase_max,phase_min_5,phase_max_5]
                names+=['snr_'+band,'snr_5_'+band,'Nbef_'+band,'Naft_'+band,'Nbef_5_'+band,'Naft_5_'+band,'phase_min_'+band,'phase_max_'+band,'phase_min_5_'+band,'phase_max_5_'+band]
    
                #r.append((lc.meta['DayMax'],lc.meta['z'],'all',snr(lc),snr(lc[idxc]),nb,na,nb5,na5,phase_min,phase_max,phase_min_5,phase_max_5))

            rsb, nsb= self.Var_Science_Book(lc,lc.meta['DayMax'],lc.meta['z'])

            ro+=rsb
            names+=nsb

            self.ro=ro
            self.names=names

        else:

            self.names=['fit_status', 'mbfit', 'salt2.z', 'salt2.t0', 'salt2.X0', 'salt2.X1', 'salt2.Color', 'salt2.CovX1X1', 'salt2.CovX1X0', 'salt2.CovX1Color', 'salt2.CovX0X0', 'salt2.CovX0Color', 'salt2.CovColorColor', 'salt2.Covt0t0', 'snr_all', 'snr_5_all', 'Nbef_all', 'Naft_all', 'Nbef_5_all', 'Naft_5_all', 'phase_min_all', 'phase_max_all', 'phase_min_5_all', 'phase_max_5_all', 'snr_g', 'snr_5_g', 'Nbef_g', 'Naft_g', 'Nbef_5_g', 'Naft_5_g', 'phase_min_g', 'phase_max_g', 'phase_min_5_g', 'phase_max_5_g', 'snr_r', 'snr_5_r', 'Nbef_r', 'Naft_r', 'Nbef_5_r', 'Naft_5_r', 'phase_min_r', 'phase_max_r', 'phase_min_5_r', 'phase_max_5_r', 'snr_i', 'snr_5_i', 'Nbef_i', 'Naft_i', 'Nbef_5_i', 'Naft_5_i', 'phase_min_i', 'phase_max_i', 'phase_min_5_i', 'phase_max_5_i', 'snr_z', 'snr_5_z', 'Nbef_z', 'Naft_z', 'Nbef_5_z', 'Naft_5_z', 'phase_min_z', 'phase_max_z', 'phase_min_5_z', 'phase_max_5_z', 'snr_y', 'snr_5_y', 'Nbef_y', 'Naft_y', 'Nbef_5_y', 'Naft_5_y', 'phase_min_y', 'phase_max_y', 'phase_min_5_y', 'phase_max_5_y', 'N_Phase_m5', 'N_Phase_p30', 'N_nights_m20_p_p30', 'Near_peak_gap', 'N_g_snrmax_10', 'N_g_snrmax_15', 'N_g_snrmax_20', 'N_r_snrmax_10', 'N_r_snrmax_15', 'N_r_snrmax_20', 'N_i_snrmax_10', 'N_i_snrmax_15', 'N_i_snrmax_20', 'N_z_snrmax_10', 'N_z_snrmax_15', 'N_z_snrmax_20', 'N_y_snrmax_10', 'N_y_snrmax_15', 'N_y_snrmax_20']
            self.ro=[-999.0]*len(self.names)

        if output_q is not None:
            output_q.put({inum : self.Summary()})  

    def Summary(self):

        #print(len(self.ro),len(self.names))
        #print(self.ro,self.names)
        t=Table([[ri] for ri in self.ro],names=tuple(self.names))
        for key in self.lc.meta.keys():
            t.add_column(Column([self.lc.meta[key]], name=key))

        return t
        
        
    def Get_Quality_LC(self,lc, T0, z):

        #estimate the number of LC points (5 sigma) before and after T0 - observer frame
        #time_begin=time.time()
        lc.sort('time')
        #self.dict_quality['SNR_tot']=5.*np.power(np.sum(np.power(lc['flux_e_sec']/lc['flux_5sigma_e_sec'],2.)),0.5)
                #print lc
        #print 'total elapse time filling table a',time.time()-time_begin
        #time_beginc=time.time()
        #calc=0
        lc_bef,lc_aft=self.Bef_Aft(lc,T0)
        for band in self.bands_rest:
            #time_beginb=time.time()
            idx=lc['band']=='LSST::'+band                   
            
            self.dict_quality['N_bef_'+band]=len(lc_bef[lc_bef['band']=='LSST::'+band])
            self.dict_quality['N_aft_'+band]=len(lc_aft[lc_aft['band']=='LSST::'+band])
            #print 'total elapse time filling table b1',time.time()-time_beginb

            #if len(lc[idx])>=1:
                
                #self.dict_quality['SNR_'+band]=5.*np.power(np.sum(np.power(lc['flux_e_sec'][idx]/lc['flux_5sigma_e_sec'][idx],2.)),0.5)
                #mean_cadence, rms_cadence=self.Get_cadence(lc[idx])
                #mean_cadence, rms_cadence=0,0
                #self.dict_quality['cadence_'+band]=mean_cadence
                #self.dict_quality['cadence_rms_'+band]=rms_cadence
                #self.dict_quality['m5sigma_'+band]=np.mean(lc[idx]['m5'])
                #self.dict_quality['m5sigma_rms_'+band]=np.std(lc[idx]['m5'])
            #calc+=time.time()-time_beginb
            #print 'total elapse time filling table b',time.time()-time_beginb
        #print 'total elapse time filling table c',time.time()-time_beginc
        self.dict_quality['N_bef_all']=int(np.sum([self.dict_quality['N_bef_'+band] for band in self.bands_rest]))
        self.dict_quality['N_aft_all']=int(np.sum([self.dict_quality['N_aft_'+band] for band in self.bands_rest]))
        #mean_cadence, rms_cadence=self.Get_cadence(lc)
        #mean_cadence, rms_cadence=0,0
        #self.dict_quality['cadence_all']=mean_cadence
        #self.dict_quality['cadence_rms_all']=rms_cadence
        
        self.dict_quality['phase_first']=(lc[0]['time']-T0)/(1.+z)
        self.dict_quality['phase_last']=(lc[-1]['time']-T0)/(1.+z)
           
        #print 'total elapse time filling table',time.time()-time_begin

    def Get_Sel_Science_Book(self,lc,T0,z):
        
        lc.sort('time')
        phases=(lc['time']-T0)/(1.+z)
        #print('phases',phases)
        idx = phases < -5.
        idxb = phases > 30.
        self.dict_quality['N_Phase_m5']=len(phases[idx])
        self.dict_quality['N_Phase_p30']=len(phases[idxb])
        #print('yep',T0,len(phases[idx]),len(phases[idxb]))


        diffa,seldiffa=self.Get_Diff(lc,phases,-20.,60.)
        diffb,seldiffb=self.Get_Diff(lc,phases,-5.,30.)

        self.dict_quality['N_nights_m20_p_p30']=len(seldiffa)-1
        #print('hello',len(diffb))
        nearpeak=-1
        if len(diffb) > 0:
            nearpeak=np.max(diffb)/(1.+z)
        self.dict_quality['Near_peak_gap']=nearpeak
        #print('yep',len(seldiffa)-1,nearpeak)

        #print('gap',np.max(diffb)/(1.+z))

        #print(lc['band'])
        snr_refs=[10.,15.,20.]
        for band in self.bands_rest:
            idf = lc['band']=='LSST::'+band
            sel_idf= lc[idf]
            if len(sel_idf) > 0:
            #print('band',band,len(sel_idf))
                snr_max=np.max(sel_idf['flux']/sel_idf['fluxerr'])
                for valref in snr_refs:
                    nres=0
                    if snr_max > valref:
                        nres=1
                    self.dict_quality['N_'+band+'_snrmax_'+str(int(valref))]=nres
        
            else:
                for valref in snr_refs:
                    self.dict_quality['N_'+band+'_snrmax_'+str(int(valref))]=0 


    def Get_Diff(self,lc,phases,phase_min,phase_max):

        idxc=(phases > phase_min)&(phases < phase_max)
        
        sel=lc[idxc]

        diff_time=1.
    
        diff=[io-jo for jo,io in zip(sel['time'][:-1], sel['time'][1:])]
        
        seldiff=[i+1 for i in range(len(diff)) if diff[i]>= diff_time]
        seldiff=[0]+seldiff+[len(sel)]

        return diff,seldiff


    def Get_nums(self, lc, T0):
        
        #idxa=(lc['time'] <= T0)&(lc['time'] > T0-20)
        #idxb=(lc['time']> T0)&(lc['time'] <= T0+40)
        idx=(lc['time'] > T0-20)&(lc['time'] <= T0+40)
        lc_sel=lc[idx]
        idxa=lc_sel['time']<=T0
        lc_bef=lc_sel[idxa]
        ntot=len(lc_sel)
        nbef=len(lc_bef)

        return nbef,ntot-nbef

    def Bef_Aft(self, lc, T0):
        """
        idxa=(lc['time'] <= T0)&(lc['time'] > T0-20)
        idxb=(lc['time']> T0)&(lc['time'] <= T0+40)
        return lc[idxa],lc[idxb]
        """
        if len(lc) > 0:
            diff=lc['time']-T0
            idxa=diff <= 0
            idxb=diff > 0
            return len(lc[idxa]),len(lc[idxb])
        else:
            return 0.0,0.0
    def Get_cadence(self, sel):

        if len(sel) == 0 or len(sel) == 1:
            return 0.0,0.0
        else:
            calc=[io-jo for jo,io in zip(sel['time'][:-1], sel['time'][1:])]  
            if len(calc)==1:
                return calc[0],0.0
            else:
                return np.mean(calc),np.std(calc)

    def Plot_LC(self,lc):
        
        figb, axb = plt.subplots(ncols=2, nrows=3, figsize=(10,9))
        
        x=dict(zip(['u','g','r','i','z','y'],[0,0,1,1,2,2]))
        y=dict(zip(['u','g','r','i','z','y'],[0,1,0,1,0,1]))   
        
        
        for band in self.bands_rest:
            idx = lc['band']=='LSST::'+band
            lc_sel=lc[idx]
            axb[x[band]][y[band]].errorbar(lc_sel['time'],lc_sel['flux'],yerr=lc_sel['fluxerr'])
            axb[x[band]][y[band]].errorbar([lc.meta['DayMax']]*2,[min(lc_sel['flux']),max(lc_sel['flux'])],color='k')
        
        plt.show()
        
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
