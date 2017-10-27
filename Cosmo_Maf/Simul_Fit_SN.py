from SN_Object import SN_Object
from Parameters import parameters
from Throughputs import Throughputs
import numpy as np
import math

from lsst.sims.photUtils import SignalToNoise
from lsst.sims.photUtils import PhotometricParameters
import sncosmo
from astropy.table import vstack,Table
import astropy.units as u
import matplotlib.pyplot as plt
from lsst.sims.photUtils.EBV import EBVbase
from lsst.sims.photUtils import Bandpass,Sed
import time

class Simul_Fit_SN():

    def __init__(self,T0,c,x1,z,observations,model='salt2-extended',version='1.0',sn_type='Ia',ra=0.,dec=0.,syste=False):
        
        self.T0=T0
        self.c=c
        self.x1=x1
        self.z=z
        self.obs=observations
        self.model=model
        self.version=version
        self.sn_type=sn_type
        self.mbsim=-1.
        self.ra_field=ra
        self.dec_field=dec
        self.params=parameters()
        self.filterNames = ['u','g','r','i','z','y']
        self.syste=syste

        #This will be the data for sncosmo fitting
        self.table_for_fit={}
        self.table_for_fit['error_coadd_opsim'] = Table(names=('time','flux','fluxerr','band','zp','zpsys'), dtype=('f8', 'f8','f8','S7','f4','S4'))
        self.table_for_fit['error_coadd_through'] = Table(names=('time','flux','fluxerr','band','zp','zpsys'), dtype=('f8', 'f8','f8','S7','f4','S4'))
        self.table_LC=Table(names=('filter','expMJD','visitExpTime','FWHMeff','moon_frac','filtSkyBrightness','kAtm','airmass','fiveSigmaDepth','Nexp','e_per_sec','e_per_sec_err'), dtype=('S7','f8', 'f8','f8','f8', 'f8','f8','i8','f8','f8','f8','f8'))

        if syste:
            vars_a=['filter','expMJD','visitExpTime','FWHMeff','moon_frac','filtSkyBrightness','kAtm','airmass','fiveSigmaDepth','Nexp','e_per_sec','e_per_sec_err','mag','err_mag','fiveSigmaThrough']
            dtype_a=['S7','f8', 'f8','f8','f8', 'f8','f8','i8','f8','f8','f8','f8','f8','f8','f8']

            for i in range(1,6,1):
                vars_a.append('err_mag_plus_'+str(i))
                vars_a.append('err_mag_minus_'+str(i))
                vars_a.append('fiveSigmaThrough_plus_'+str(i))
                vars_a.append('fiveSigmaThrough_minus_'+str(i))
                dtype_a.append('f8')
                dtype_a.append('f8')
                dtype_a.append('f8')
                dtype_a.append('f8')

            vars_a.append('z')
            dtype_a.append('f8')

            self.table_LC_syste=Table(names=tuple(vars_a), dtype=tuple(dtype_a))

        self.transmission=Throughputs()

        # Register LSST band pass (system) in sncosmo

        for filtre in self.filterNames:
            band=sncosmo.Bandpass(self.transmission.lsst_system[filtre].wavelen, self.transmission.lsst_system[filtre].sb, name='LSST::'+filtre,wave_unit=u.nm)
            sncosmo.registry.register(band, force=True)

        #print 'hello',self.ra_field,self.dec_field
        self.SN=SN_Object(ra=np.rad2deg(self.ra_field),dec=np.rad2deg(self.dec_field),z=z,t0=T0,c=c,x1=x1,model=self.model,version=self.version,sn_type=self.sn_type)
        #print 'SN X0',self.SN.SN.get('x0')
        if self.SN.sn_type=='Ia':
            self.mbsim=self.SN.SN._source.peakmag('bessellb','vega')

        self.outdict={}
        self.outdict['t0']=self.T0
        self.outdict['c']=self.c
        self.outdict['x1']=self.x1
        self.outdict['z']=self.z
        self.outdict['ra']=self.ra_field
        self.outdict['dec']=self.dec_field
        self.outdict['status']='unkown'
        self.outdict['fit']=None
        self.outdict['mbsim']=-999.
        self.outdict['sn_type']=self.SN.sn_type
        self.outdict['sn_model']=self.SN.model
        self.outdict['sn_version']=self.SN.version
        self.outdict['mbsim']=self.mbsim
        self.outdict['x0']=self.SN.SN.get('x0')
        self.outdict['dL']=self.SN.lumidist
        
        m_begin_date=self.obs['mjd'].min()
        m_end_date=self.obs['mjd'].max()

        lowrange=-30
        highrange=50
        timelow = self.T0+ lowrange*(1+self.z)
        timehigh = self.T0 + highrange*(1+self.z)

        self.obs.sort('mjd')
        self.obs=self.obs[np.where(np.logical_and(self.obs['mjd']>=timelow,self.obs['mjd']<=timehigh))]

        self.Simulate_LC()
        """
        for band in self.filterNames:
            sel=self.table_LC[np.where(self.table_LC['filter']==band)]
            print sel
        """
        #self.Fit_LCs()
    def Simulate_LC(self):

        sed_SN=self.SN.get_SED(self.obs['mjd'])
            
        print 'yes'
        print 'dot',np.dot(sed_SN.wavelen,sed_SN.flambda)
        
    def Simulate_LC_last(self):

        
        sed_SN=self.SN.get_SED(self.obs['mjd'])
        
        #print 'SED loaded'
        for i in range(len(sed_SN.wavelen)):
            obs=self.obs[i]
            sed=Sed(wavelen=sed_SN.wavelen[i],flambda=sed_SN.flambda[i])
            visittime=obs['exptime']
            filtre=obs['band'][-1]
            photParams = PhotometricParameters(nexp=visittime/15.)
            e_per_sec = sed.calcADU(bandpass=self.transmission.lsst_atmos_aerosol[filtre], photParams=photParams) #number of ADU counts for expTime
            e_per_sec/=visittime/photParams.gain

            flux_SN=sed.calcFlux(bandpass=self.transmission.lsst_atmos_aerosol[filtre])

            FWHMeff=obs['FWHMeff']
            if flux_SN >0:
                 
                mag_SN=-2.5 * np.log10(flux_SN / 3631.0)
                m5_calc,snr_m5_through=self.Get_m5(filtre,mag_SN,obs['sky'],photParams,FWHMeff)
                
                self.table_LC.add_row(('LSST::'+filtre,obs['mjd'],visittime,FWHMeff,obs['moon_frac'],obs['sky'],obs['kAtm'],obs['airmass'],obs['m5sigmadepth'],obs['Nexp'],e_per_sec,e_per_sec/snr_m5_through))

        
    def Simulate_LC_old(self):

         for obs in self.obs:
             filtre=obs['band'].split('::')[1]
             
             #seeing=obs['rawSeeing']
             time_obs=obs['mjd']
             #m5_opsim=obs['fiveSigmaDepth']
             m5_opsim=obs['m5sigmadepth']
             sed_SN=self.SN.get_SED(time_obs)
             
             #self.transmission.Load_Atmosphere(obs['airmass'])
             #
             flux_SN=sed_SN.calcFlux(bandpass=self.transmission.lsst_atmos_aerosol[filtre])
             #flux_SN=sed_SN.calcFlux(bandpass=self.transmission.lsst_system[filtre])

             print 'hello flux',flux_SN
             #visittime=obs['visitExpTime']
             visittime=obs['exptime']

             photParams = PhotometricParameters(nexp=visittime/15.)
             #e_per_sec = sed_SN.calcADU(bandpass=self.transmission.lsst_system[filtre], photParams=photParams) #number of ADU counts for expTime
             e_per_sec = sed_SN.calcADU(bandpass=self.transmission.lsst_atmos_aerosol[filtre], photParams=photParams) #number of ADU counts for expTime
             e_per_sec/=visittime/photParams.gain
             
             """
             Filter_Wavelength_Correction = np.power(500.0 / self.params.filterWave[filtre], 0.3)
             Airmass_Correction = math.pow(obs['airmass'],0.6)
             FWHM_Sys = self.params.FWHM_Sys_Zenith * Airmass_Correction
             FWHM_Atm = seeing * Filter_Wavelength_Correction * Airmass_Correction
             finSeeing = self.params.scaleToNeff * math.sqrt(np.power(FWHM_Sys,2) + self.params.atmNeffFactor * np.power(FWHM_Atm,2))
             
             FWHMeff = SignalToNoise.FWHMgeom2FWHMeff(finSeeing)
             """
             #print 'alors pal',finSeeing,FWHMeff,obs['FWHMgeom'],obs['FWHMeff'],SignalToNoise.FWHMgeom2FWHMeff(obs['FWHMgeom'])
             FWHMeff=obs['FWHMeff']
             if flux_SN >0:
                  
                 mag_SN=-2.5 * np.log10(flux_SN / 3631.0)
             
                 #print 'hello',finSeeing,FWHMeff

                 m5_calc,snr_m5_through=self.Get_m5(filtre,mag_SN,obs['sky'],photParams,FWHMeff)
                 
                 m5_opsim+=1.25*np.log10(visittime/30.)
                 snr_m5_opsim,gamma_opsim=SignalToNoise.calcSNR_m5(mag_SN,self.transmission.lsst_atmos_aerosol[filtre],m5_opsim,photParams)

                 err_flux_SN_opsim=flux_SN/snr_m5_opsim
                 err_flux_SN_through=flux_SN/snr_m5_through

                 self.table_for_fit['error_coadd_opsim'].add_row((time_obs,flux_SN,err_flux_SN_opsim,'LSST::'+filtre,25,'ab'))
                 self.table_for_fit['error_coadd_through'].add_row((time_obs,flux_SN,err_flux_SN_through,'LSST::'+filtre,25,'ab'))
                 """
                 flatSed_err = Sed()
                 flatSed_err.setFlatSED(wavelen_min, wavelen_max, wavelen_step)
                 #flatSed_err.multiplyFluxNorm(err_flux_SN_through)
                 err_mag_SN=-2.5 * np.log10( err_flux_SN_through/ 3631.0)
                 flux0_err=np.power(10.,-0.4*err_mag_SN)
                 flatSed_err.multiplyFluxNorm(flux0_err)
                 e_per_sec_err = flatSed_err.calcADU(bandpass=self.transmission.lsst_atmos_aerosol[filtre], photParams=photParams) #number of ADU counts for expTime
                 e_per_sec_err/=obs['visitExpTime']/2.3
                 """
                 #print 'test',e_per_sec_err,e_per_sec/snr_m5_through
                 #print 'hello',filtre,obs['expMJD'],obs['visitExpTime'],obs['rawSeeing'],obs['moon_frac'],obs['filtSkyBrightness'],obs['kAtm'],obs['airmass'],obs['fiveSigmaDepth'],obs['Nexp'],e_per_sec,e_per_sec_err,flux_SN,err_flux_SN_through
                 self.table_LC.add_row(('LSST::'+filtre,obs['mjd'],visittime,FWHMeff,obs['moon_frac'],obs['sky'],obs['kAtm'],obs['airmass'],obs['m5sigmadepth'],obs['Nexp'],e_per_sec,e_per_sec/snr_m5_through))
                 if self.syste:
                     resu=['LSST::'+filtre,obs['mjd'],visittime,FWHMeff,obs['moon_frac'],obs['sky'],obs['kAtm'],obs['airmass'],m5_opsim,obs['Nexp'],e_per_sec,e_per_sec/snr_m5_through,mag_SN,mag_SN/snr_m5_through,m5_calc]
                    
                     for i in range(1,6,1):
                         m5_calc_plus,snr_m5_plus=self.Get_m5(filtre,mag_SN,obs['sky']+float(i)/10.,photParams,FWHMeff)
                         m5_calc_minus,snr_m5_minus=self.Get_m5(filtre,mag_SN,obs['sky']-float(i)/10.,photParams,FWHMeff)
                         resu.append(mag_SN/snr_m5_plus)
                         resu.append(mag_SN/snr_m5_minus)
                         resu.append(m5_calc_plus)
                         resu.append(m5_calc_minus)
                     resu.append(self.z)
                     self.table_LC_syste.add_row(tuple(resu))
             else:
                 self.table_LC.add_row(('LSST::'+filtre,obs['mjd'],visittime,FWHMeff,obs['moon_frac'],obs['sky'],obs['kAtm'],obs['airmass'],obs['m5sigmadepth'],obs['Nexp'],-1.,-1.))   


         self.table_LC.sort(['filter','expMJD'])

    def Fit_LCs(self):
        

        dict_fit={}
        for val in ['error_coadd_opsim','error_coadd_through']: 
            dict_fit[val]={}
            dict_fit[val]['sncosmo_fitted']={}
            dict_fit[val]['table_for_fit']=self.table_for_fit[val]
            res,fitted_model,mbfit,fit_status=self.Fit_LC(self.table_for_fit[val])
            #print res,fitted_model
            if res is not None:
                dict_fit[val]['sncosmo_res']=res
                for i,par in enumerate(fitted_model.param_names):
                    dict_fit[val]['sncosmo_fitted'][par]=fitted_model.parameters[i]
                dict_fit[val]['mbfit']=mbfit
            dict_fit[val]['fit_status']=fit_status

        self.outdict['fit']=dict_fit
        self.outdict['status']='try_fit'

    def Fit_LC(self,t):

        
        select=t[np.where(np.logical_and(t['flux']/t['fluxerr']>5.,t['flux']>0.))]
        #print 'fit',len(select)
 
        if len(select)>=5:
            #print 'data to be fitted',select
            try:
                
                
                #print 'before fit here'
                #print 'SN parameters',self.SN.SN
                
                #print 'fitting',select
                
                z_sim=self.SN.z
                #print 'hello z',z_sim
                #print 'fit it',val,time.time()-self.thetime
                res, fitted_model = sncosmo.fit_lc(select, self.SN.SN_fit_model,['z', 't0', 'x0', 'x1', 'c'],bounds={'z':(z_sim-0.01, z_sim+0.01)})
                #res, fitted_model = sncosmo.fit_lc(select, self.SN.SN,['t0', 'x0', 'x1', 'c'])
    

                #print 'after fit',res.keys()
                #print res.keys()
                
                """
                print 'after fit'
                print res['parameters'],res['errors']
                """

                mbfit=fitted_model._source.peakmag('bessellb','vega')
                #print 'oooo test',-2.5*np.log10(res['parameters'][2])+10.635,fitted_model.bandmag('bessellb','vega',res['parameters'][1]),mbsim,mbfit,mbsim-mbfit

                """
                sncosmo.plot_lc(t, model=fitted_model,color='k',pulls=False)
                
                plt.show()
                """
                #print 'fitted'
                return res,fitted_model,mbfit,'ok'
            
            except (RuntimeError, TypeError, NameError):
                #print 'crashed'
                return None,None,-1,'crash'

        else:
            return None,None,-1,'Noobs'

    def Get_m5(self,filtre,mag_SN,msky,photParams,FWHMeff):

        wavelen_min, wavelen_max, wavelen_step=self.transmission.lsst_system[filtre].getWavelenLimits(None,None,None)
                
        flatSed = Sed()
        flatSed.setFlatSED(wavelen_min, wavelen_max, wavelen_step)
        flux0=np.power(10.,-0.4*msky)
        flatSed.multiplyFluxNorm(flux0)
        m5_calc=SignalToNoise.calcM5(flatSed,self.transmission.lsst_atmos_aerosol[filtre],self.transmission.lsst_system[filtre],photParams=photParams,FWHMeff=FWHMeff)
        snr_m5_through,gamma_through=SignalToNoise.calcSNR_m5(mag_SN,self.transmission.lsst_atmos_aerosol[filtre],m5_calc,photParams)
        
        return m5_calc,snr_m5_through
                 
