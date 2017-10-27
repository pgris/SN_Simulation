from Fit_LC import *
import cPickle as pkl
import matplotlib.pyplot as plt
from astropy.table import vstack,Table,Column

class Fit_Single_LC:
    def __init__(self, lc,telescope,inum,output_q):
        
        self.lc=lc
        self.telescope=telescope
        print lc.dtype
        self.bands_rest = 'grizy'
        self.dict_quality={}
        self.outdict={}
        """
        print lc.meta
        self.T0=lc.meta['DayMax'].value
        self.z=lc.meta['z'][0]
        """

        for band in self.bands_rest:
            self.dict_quality['N_bef_'+band]=0.
            self.dict_quality['N_aft_'+band]=0.
            self.dict_quality['cadence_'+band]=0.
            self.dict_quality['cadence_rms_'+band]=0.
            self.dict_quality['m5sigma_'+band]=0.
            self.dict_quality['m5sigma_rms_'+band]=0.
            self.dict_quality['SNR_'+band]=0.0
        self.dict_quality['N_bef_all']=0.
        self.dict_quality['N_aft_all']=0.
        self.dict_quality['phase_first']=0.
        self.dict_quality['phase_last']=0.
        self.dict_quality['cadence_all']=0.
        self.dict_quality['cadence_rms_all']=0.
        self.dict_quality['SNR_tot']=0.0

        if self.lc is not None:
            idx = self.lc['flux']>0.
            lc_sel=self.lc[idx]
            idxa = lc_sel['flux']/lc_sel['fluxerr']>5.
            lc_sel=lc_sel[idxa]
            if len(lc_sel)>0.:
                self.Get_Quality_LC(lc_sel,lc.meta['DayMax'][0],lc.meta['z'][0])
                #print self.dict_quality
                #self.Plot_LC(lc_sel)
        
        #print [val for key,val in self.dict_quality.items()]
        #print [key for key,val in self.dict_quality.items()]

                if (self.dict_quality['N_bef_all']+self.dict_quality['N_aft_all']) >= 5:
                    self.outdict['status']='go_fit'
                    self.outdict['fit_status']='unknow'
                    self.Fit_LC(lc_sel)
                    print self.outdict.keys()
                    self.Get_Quality_LC(lc_sel,self.outdict['sncosmo_fitted']['t0'],lc.meta['z'][0])
                else:
                    self.outdict['status']='no_obs'
                    self.outdict['fit_status']='unknow'
        else:
            self.outdict['status']='no_pha'
            self.outdict['fit_status']='unknow' 


        """
        tab_resu=Table(names=tuple([key for key,val in self.dict_quality.items()]),dtype=tuple(['f8']*len(self.dict_quality)),meta=lc.meta)
        tab_resu.add_row(tuple([val for key,val in self.dict_quality.items()]))

        print tab_resu
        print tab_resu.meta
        print self.outdict
        """

        output_q.put({inum : self.Summary()})

    def Get_Quality_LC(self,lc, T0, z):

        #estimate the number of LC points (5 sigma) before and after T0 - observer frame
       
        
        lc.sort('time')
        n_bef_tot=0
        n_aft_tot=0
        self.dict_quality['SNR_tot']=5.*np.power(np.sum(np.power(lc['flux_e_sec']/lc['flux_5sigma_e_sec'],2.)),0.5)
                #print lc
        n_bef_tot, n_aft_tot=self.Get_nums(lc,T0)
        for band in self.bands_rest:
            idx=lc['band']=='LSST::'+band
            n_bef, n_aft=self.Get_nums(lc[idx],T0)                   
            self.dict_quality['N_bef_'+band]=n_bef
            self.dict_quality['N_aft_'+band]=n_aft
            
            if len(lc[idx])>=1:
                self.dict_quality['SNR_'+band]=5.*np.power(np.sum(np.power(lc['flux_e_sec'][idx]/lc['flux_5sigma_e_sec'][idx],2.)),0.5)
                mean_cadence, rms_cadence=self.Get_cadence(lc[idx])
                self.dict_quality['cadence_'+band]=mean_cadence
                self.dict_quality['cadence_rms_'+band]=rms_cadence
                self.dict_quality['m5sigma_'+band]=np.mean(lc[idx]['m5'])
                self.dict_quality['m5sigma_rms_'+band]=np.std(lc[idx]['m5'])

        self.dict_quality['N_bef_all']=n_bef_tot 
        self.dict_quality['N_aft_all']=n_aft_tot
        mean_cadence, rms_cadence=self.Get_cadence(lc)
        self.dict_quality['cadence_all']=mean_cadence
        self.dict_quality['cadence_rms_all']=rms_cadence
        
        self.dict_quality['phase_first']=(lc[0]['time']-T0)/(1.+z)
        self.dict_quality['phase_last']=(lc[len(lc)-1]['time']-T0)/(1.+z)
           

    def Get_nums(self, lc, T0):
        
        idxa=np.logical_and(lc['time'] <= T0,lc['time'] > T0-20)
        idxb=np.logical_and(lc['time']> T0,lc['time'] <= T0+40)

        return float(len(lc[idxa])),float(len(lc[idxb]))

    def Get_cadence(self, sel):

        if len(sel) == 0:
            return 0.0,0.0
        else:
            if len(sel)==1:
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
            axb[x[band]][y[band]].errorbar([lc.meta['DayMax'][0]]*2,[min(lc_sel['flux']),max(lc_sel['flux'])],color='k')
        
        plt.show()

    def Fit_LC(self,lc):
        
        
        myfit=Fit_LC(z=lc.meta['z'][0],telescope=self.telescope,Plot=False)
        #print 'going to fit',len(self.lc),self.lc.dtype

        #self.outdict['m5sigma']=np.median(self.lc['m5'])
        self.outdict['airmass']=np.median(lc['airmass'])
                                          
        res,fitted_model,mbfit,covar_mb,fit_status=myfit(lc)

        #print 'hello',res,fitted_model,mbfit,fit_status
        if fit_status == 'ok':
            #print 'hello',res,fitted_model
            self.outdict['sncosmo_res']=res
        
            self.outdict['sncosmo_fitted']={}
            self.outdict['recalc']={}
            for i,par in enumerate(fitted_model.param_names):
                self.outdict['sncosmo_fitted'][par]=fitted_model.parameters[i]
                #self.outdict['fitted_model']=fitted_model
            for key, val in covar_mb.items():
                self.outdict['recalc'][key]=val
        
            #print 'allo',covar_mb
            self.outdict['mbfit']=mbfit
            self.outdict['fit_status']='fit_ok'

        if fit_status == 'crash':
            self.outdict['fit_status']='crashd'


    def Summary(self):
        
        resu=self.dict_quality

        #print 'yes',self.duration
            
        resu['status']=self.outdict['status']
        resu['fit_status']=self.outdict['fit_status']

        resu['salt2.T0']=-999.
        resu['salt2.X0']=-999.
        resu['salt2.X1']=-999.
        resu['salt2.Color']=-999.
        resu['salt2.CovX1X1']=-999.
        resu['salt2.CovX0X0']=-999.
        resu['salt2.Covmbmb']=-999.
        resu['salt2.CovColorColor']=-999.
        resu['salt2.CovX0X1']=-999.
        resu['salt2.CovColorX0']=-999.
        resu['salt2.CovColorX1']=-999.
        resu['salt2.CovColormb']=-999.
        resu['salt2.CovX1mb']=-999.

        resu['mbfit']=-999.

        if self.outdict['status']=='go_fit':

            if self.outdict['fit_status']=='fit_ok':
                corr={}
                for i,pal in enumerate(self.outdict['sncosmo_res']['vparam_names']):
                    corr[pal]=i
                #print 'hhe',i,pal
                resu['salt2.T0']=self.outdict['sncosmo_fitted']['t0']
                resu['salt2.X0']=self.outdict['sncosmo_fitted']['x0']
                resu['salt2.X1']=self.outdict['sncosmo_fitted']['x1']
                resu['salt2.Color']=self.outdict['sncosmo_fitted']['c']
                if self.outdict['sncosmo_res']['covariance'] is not None:
                    resu['salt2.CovX1X1']=self.outdict['sncosmo_res']['covariance'][corr['x1']][corr['x1']]
                    resu['salt2.CovColorColor']=self.outdict['sncosmo_res']['covariance'][corr['c']][corr['c']]
                    resu['salt2.CovX0X0']=self.outdict['sncosmo_res']['covariance'][corr['x0']][corr['x0']]
                    resu['salt2.CovX0X1']=self.outdict['sncosmo_res']['covariance'][corr['x0']][corr['x1']]
                    resu['salt2.CovColorX0']=self.outdict['sncosmo_res']['covariance'][corr['c']][corr['x0']]
                    resu['salt2.CovColorX1']=self.outdict['sncosmo_res']['covariance'][corr['c']][corr['x1']]
                    for val in ['salt2.Covmbmb','salt2.CovColormb','salt2.CovX1mb']:
                        resu[val]=self.outdict['recalc'][val]

                resu['mbfit']=self.outdict['mbfit']
                #print 'ohohoho',resu['mbfit']-resu['mbsim']
        
        #print 'hello',resu.values(),resu.keys()
        #restab=np.rec.fromrecords(tuple([res for res in resu.values()]),names=[key for key in resu.keys()])
        """
        print restab,restab.dtype.names
        for name in restab.dtype.names:
            print name,restab[name]
       
        print restab
        print 'hhhh',len(resu),self.outdict['observations']
        """
        #print 'test here'
        t=Table(meta=self.lc.meta)
        for key,val in resu.items():
            aa = Column([val], name=key)
            t.add_column(aa)

        #print t
        return t
