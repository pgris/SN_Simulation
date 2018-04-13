import sncosmo
import numpy as np
from astropy import (cosmology, units as u, constants as const)
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
from lsst.sims.photUtils import Bandpass,Sed
import os
from Telescope import *
from scipy import interpolate, integrate
from lsst.sims.photUtils import PhotometricParameters
from scipy.interpolate import InterpolatedUnivariateSpline as Spline1d
import time
from scipy.interpolate import griddata

class Fit_LC:
    def __init__(self,model='salt2-extended',version='1.0',z=0.1,telescope=None,Plot=False,bands='ugrizy'):

        self.Plot=Plot
        #transmission=telescope.throughputs
        self.z=z
        #bands=[b[-1:] for b in np.unique(select['band'])]
        self.model=model
        self.version=version

        #print 'hello',bands,telescope.airmass
        """
        for filtre in bands:
            if telescope.airmass > 0:
                band=sncosmo.Bandpass(transmission.atmosphere[filtre].wavelen,transmission.atmosphere[filtre].sb, name='LSST::'+filtre,wave_unit=u.nm)
            else:
                band=sncosmo.Bandpass(transmission.system[filtre].wavelen,transmission.system[filtre].sb, name='LSST::'+filtre,wave_unit=u.nm) 
            sncosmo.registry.register(band, force=True)
        """
        source=sncosmo.get_source(model,version=version)
        dust = sncosmo.OD94Dust()

          
        self.SN_fit_model=sncosmo.Model(source=source)
        self.SN_fit_model.set(z=self.z)
        #SN_fit_model.set_source_peakabsmag(peakAbsMagBesselB, 'bessellB', 'vega',cosmo=self.astropy_cosmo)

    def __call__(self,meas):

        #time_begin=time.time()
        t=meas
        select=t[np.where(np.logical_and(t['flux']/t['fluxerr']>5.,t['flux']>0.))]
        #print 'hello',select.meta
        #idx=select['band']!='LSST::u'
        #select=select[idx]

        #print 'I will fit',meas
        """
        if self.z > 0.35:
           idx=select['band']!='LSST::g'
           select=select[idx] 

        if self.z > 0.85:
           idx=select['band']!='LSST::r'
           select=select[idx]

        if self.z > 1.27:
           idx=select['band']!='LSST::i'
           select=select[idx]
        """
        #print 'what I have to fit',len(select)
        try:
            #print 'trying to fit',len(select)
            fit_status='notok'
            if len(select) > 0:
                res, fitted_model = sncosmo.fit_lc(select, self.SN_fit_model,['t0', 'x0', 'x1', 'c'],bounds={'z':(self.z-0.001, self.z+0.001)})
            #print 'total elapse time fit',time.time()-time_begin
            #self.sigma_c=res['errors']['c']
                mbfit=fitted_model._source.peakmag('bessellb','vega')
            #mbfit=0.
                params={}
            #print res
                for i,par in enumerate(fitted_model.param_names):
                    params[par]=fitted_model.parameters[i]
                #print(i,par)
            #snutils=SN_Utils()
            #mbfit_calc=snutils.mB(params)
            #print 'total elapse time mbfit',time.time()-time_begin,mbfit
                fit_status='fitok'
            else:
                res=None
                fitted_model=None
                mbfit=-1
                fit_status='nodat'
            """
            
            
            covar_mb=snutils.Covar(params,res['covariance'],res['vparam_names'])
            print covar_mb
            
            #print 'fitted',res['vparam_names']
            covar_mb={} 
            what=(params['z'],params['x0'],params['x1'],params['c'])
            print what
            Der=np.zeros(shape=(len(res['vparam_names']),1))
            ider=-1
            for i,key in enumerate(res['vparam_names']):
                ider+=1
                if key == 't0':
                    Der[ider]=0.
                else:
                    Der[ider]=griddata((deriv_mb['z'],deriv_mb['X0'],deriv_mb['X1'],deriv_mb['Color']), deriv_mb['dmb_d'+key],what , method='nearest')
                

            print 'hhh',res['covariance'],Der
            Prod=np.dot(res['covariance'],Der)

            #print 'prod',Prod
            for i,key in enumerate(res['vparam_names']):
                if key != 'c':
                    covar_mb['salt2.Cov'+key.upper()+'mb']=Prod[i,0]
                else:
                    covar_mb['salt2.CovColormb']=Prod[i,0] 

            covar_mb['salt2.Covmbmb']=np.asscalar(np.dot(Der.T,Prod))

            print covar_mb
            #print 'total elapse time covar',time.time()-time_begin
            #print 'after fit',mbfit,mbfit_calc,mbfit-mbfit_calc
            #snutils.Test()
            """
            if self.Plot:
                sncosmo.plot_lc(select, model=fitted_model,color='r',pulls=False,errors=res.errors) 
                plt.show()
            #print 'total elapse time',time.time()-time_begin
            #return res,fitted_model,mbfit,covar_mb,'ok'
            return res,fitted_model,mbfit,fit_status   

        except (RuntimeError, TypeError, NameError):
                """
                print select
                self.Plot_bands(select)
                plt.show()
                self.sigma_c=0.
                print 'crashed'
                """
                #print 'crashed'
                #self.Plot_bands(select)
                #plt.show()
                #return None,None,-1,None,'crash'
                return None,None,-1,'crash'

    @property
    def sigma_c(self):
        return self.sigma_c

    def Plot_bands(self,obj,time_name='time',flux_name='flux',errflux_name='fluxerr',filter_name='band',T0=0):
        
        figa, axa = plt.subplots(ncols=2, nrows=3, figsize=(10,9))
    
        for j,band in enumerate(['u','g','r','i','z','y']):
            if j<2:
                k=0
            if j>= 2 and j < 4:
                k=1
            if j>=4:
                k=2

            selobs=obj[np.where(obj[filter_name]=='LSST::'+band)]

            axa[k][j%2].errorbar(selobs[time_name],selobs[flux_name],yerr=selobs[errflux_name],fmt='.',ecolor='r',color='r')
           
        axa[k][j%2].set_xlim(T0-30,T0+50)

