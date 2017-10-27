from astropy.table import vstack,Table
import cPickle as pkl
from astropy.io import ascii
import numpy as np
#from Simul_Fit_SN import *
import matplotlib.pyplot as plt
import sncosmo
from astropy import (cosmology, units as u, constants as const)
from astropy.cosmology import FlatLambdaCDM
from Telescope import *
from lsst.sims.photUtils import Bandpass,Sed
from lsst.sims.photUtils import SignalToNoise
from lsst.sims.photUtils import PhotometricParameters
from lsst.sims.photUtils.EBV import EBVbase
from saunerie import instruments,salt2
import time

class Generate_LC:
    def __init__(self,parameters,fit=False,model='salt2-extended',version='1.0',telescope=None,airmass=-1):

        self.lc=[]
        self.m5={'u':23.61,'g':24.83,'r':24.35,'i':23.88,'z':23.30,'y':22.43}
        self.radeg=np.rad2deg(6.0979440)
        self.decdeg=np.rad2deg(-1.1051600)

        self.params=Table(names=('t0','c','x1','z','ra','dec','status','fit','sn_type','sn_model','sn_version','mbsim','x0','dL'),dtype=('f8','f8','f8','f8','f8','S8','S8','f8','S8','S8','S8','f8','f8','f8'))
       
        self.table_for_fit = Table(names=('time','flux','fluxerr','band','zp','zpsys'), dtype=('f8', 'f8','f8','S7','f4','S4'))
        self.model=model
        self.version=version

        self.peakAbsMagBesselB=-19.0906

        

        if self.model == 'salt2-extended':
            model_min=300.
            model_max=180000.
            wave_min=3000.
            wave_max=11501.

        if self.model=='salt2':
            model_min=3400.
            model_max=11501.
            wave_min=model_min
            wave_max=model_max

        self.wave= np.arange(wave_min,wave_max,1.)
        
        sn_type='Ia'
        source=sncosmo.get_source(self.model,version=self.version)
        #self.mycosmology=FlatLambdaCDM(H0=70, Om0=0.25)
        #astropy_cosmo=FlatLambdaCDM(H0= self.mycosmology.H0, Om0=self.mycosmology.Om0)

        dust = sncosmo.OD94Dust()

        #self.transmission=Throughputs(through_dir='FAKE_THROUGH',atmos_dir='FAKE_THROUGH',atmos=False,aerosol=False)
        
        self.transmission=telescope.throughputs

        """
        for band in self.bands:
            themax=np.max(self.transmission.lsst_system[band].sb)
            idx = self.transmission.lsst_system[band].sb >= 0.2*themax
            sel=self.transmission.lsst_system[band].wavelen[idx]
            print band,np.min(sel),np.max(sel)
        """
        self.airmass=airmass
        if self.airmass > 0:
            self.transmission.Load_Atmosphere(airmass)

        self.lc={}
        #print 'there we go',parameters,len(parameters),parameters.dtype

        self.param=parameters
        """
        SN=sncosmo.Model(source=source,effects=[dust, dust],
                         effect_names=['host', 'mw'],
                         effect_frames=['rest', 'obs'])
        """
        self.SN=sncosmo.Model(source=source)
        
        self.z=self.param['z']
        self.Cosmology()
        lumidist=self.astropy_cosmo.luminosity_distance(self.param['z']).value*1.e3
        X0 = self.X0_norm() / lumidist** 2
        #print 'before alpha beta',X0
        alpha=0.13
        beta=3.
        X0 *= np.power(10., 0.4*(alpha*self.param['X1'] -beta*self.param['Color']))
        #print 'llla',X0,alpha,beta,param['X1'],param['Color'],param['z'],lumidist
        self.X0=X0
        #SN=sncosmo.Model(source=source)
        self.SN.set(z=self.param['z'])
        self.SN.set(t0=self.param['DayMax'])
        self.SN.set(c=self.param['Color'])
        self.SN.set(x1=self.param['X1'])
        self.SN.set(x0=self.X0)
        #self.SED={}
        print 'sncosmo parameters',self.SN.param_names,self.SN.parameters
        #lsstmwebv = EBVbase()
        """
        ebvofMW = lsstmwebv.calculateEbv(
            equatorialCoordinates=np.array([[np.radians(self.radeg)], [np.radians(self.decdeg)]]))[0]
        SN.set(mwebv=ebvofMW)
        """
        #SN.set_source_peakabsmag(self.peakAbsMagBesselB, 'bessellB', 'vega',cosmo=self.astropy_cosmo)
        
    def __call__(self,mjds,filtre,expTime,out_q):

        fluxes=10.*self.SN.flux(mjds,self.wave)
        
        wavelength=self.wave/10.
        
        wavelength=np.repeat(wavelength[np.newaxis,:], len(fluxes), 0)
        SED_time = Sed(wavelen=wavelength, flambda=fluxes)
        
        r=[]
        visittime=expTime
        if self.airmass > 0.:
            trans=self.transmission.atmosphere[filtre]
        else:
            trans=self.transmission.system[filtre]

                
        photParams = PhotometricParameters(nexp=visittime/15.)
        itemindex = np.where(mjds==0.)

        #print 'hello',itemindex,mjds[itemindex]
        for i in range(len(SED_time.wavelen)):
                
            sed=Sed(wavelen=SED_time.wavelen[i],flambda=SED_time.flambda[i])
            
            """
            if filtre == 'r' and i==itemindex[0]:
                filtercolors = {'u':'b', 'g':'c', 'r':'g', 'i':'y', 'z':'r', 'y':'m'}
                plt.plot(SED_time.wavelen[i],SED_time.flambda[i],ls='-',label='z ='+str(self.param['z']))
                plt.ylabel('Flux density [$ergs/cm^2/s/nm$]')
                plt.xlabel('$\lambda$ [nm]')
                plt.legend(loc='upper right')
                plt.title('x0= '+str(self.X0)+' x1= '+str(self.param['X1'])+' c= '+str(self.param['Color']))
                #self.SED[self.param['z']]=(SED_time.wavelen[i],SED_time.flambda[i])
                
                if self.param['z']==0.5:
                    ax = plt.gca().twinx()
                    for i,band in enumerate(['u','g','r','i','z','y']):
                        ax.plot(self.transmission.system[filtre].wavelen,self.transmission.system[filtre].sb,linestyle='--',color=filtercolors[filtre])
            """
            flux_SN=sed.calcFlux(bandpass=self.transmission.system[filtre])
            if flux_SN >0:
                mag_SN=-2.5 * np.log10(flux_SN / 3631.0)  
                snr_m5_opsim,gamma_opsim=SignalToNoise.calcSNR_m5(mag_SN,trans,self.m5[filtre],photParams)
                err_flux_SN=flux_SN/snr_m5_opsim
                e_per_sec = sed.calcADU(bandpass=trans, photParams=photParams) #number of ADU counts for expTime
                    #e_per_sec = sed.calcADU(bandpass=self.transmission.lsst_atmos[filtre], photParams=photParams)
                e_per_sec/=visittime/photParams.gain
                #print 'ref',filtre,i,mjds[i],e_per_sec
                    #self.lc[filtre].append(e_per_sec)
                r.append((e_per_sec,mjds[i],flux_SN))
                
                self.table_for_fit.add_row((mjds[i],flux_SN,err_flux_SN,'LSST::'+filtre,25,'ab'))

            #print 'hello',r
        self.lc[filtre]=np.rec.fromrecords(r, names = ['flux','mjd','flux_SN'])

        #print 'yyyy',len(self.lc[filtre]),self.lc[filtre]['flux_SN']
       
        
        res=np.rec.fromrecords(r, names = ['flux','mjd','flux_SN'])
        out_q.put({filtre :(res,self.table_for_fit)})
        
        return res
        #plt.show()
    def Fit(self):

        #print self.table_for_fit
        t=self.table_for_fit
        select=t[np.where(np.logical_and(t['flux']/t['fluxerr']>5.,t['flux']>0.))]
        idx=select['band']!='LSST::u'
        select=select[idx]

        if self.z > 0.35:
           idx=select['band']!='LSST::g'
           select=select[idx] 


        for filtre in bands:
            band=sncosmo.Bandpass(self.transmission.atmosphere[filtre].wavelen, self.transmission.atmosphere[filtre].sb, name='LSST::'+filtre,wave_unit=u.nm)
            sncosmo.registry.register(band, force=True)
        
        source=sncosmo.get_source(self.model,version=self.version)
        dust = sncosmo.OD94Dust()
        """
        SN_fit_model=sncosmo.Model(source=source,effects=[dust, dust],
                              effect_names=['host', 'mw'],
                              effect_frames=['rest', 'obs'])
        """
        SN_fit_model=sncosmo.Model(source=source)
        SN_fit_model.set(z=self.z)
        SN_fit_model.set_source_peakabsmag(self.peakAbsMagBesselB, 'bessellB', 'vega',cosmo=self.astropy_cosmo)
        """
        lsstmwebv = EBVbase()
        ebvofMW = lsstmwebv.calculateEbv(
            equatorialCoordinates=np.array([[np.radians(self.radeg)], [np.radians(self.decdeg)]]))[0]
        
        SN_fit_model.set(mwebv=ebvofMW)
        """

        z_sim=self.z
        self.sigma_c=0.
        try:
            res, fitted_model = sncosmo.fit_lc(select, SN_fit_model,['t0', 'x0', 'x1', 'c'],bounds={'z':(z_sim-0.1, z_sim+0.1)})

            self.sigma_c=res['errors']['c']
            """
            print z_sim,self.sigma_c
            if z_sim>=0.409 and z_sim<0.411:
                sncosmo.plot_lc(select, model=fitted_model,color='k',pulls=True,errors=res.errors) 
                
                plt.show()
            """
        except (RuntimeError, TypeError, NameError):
                print 'crashed'
                for sel in select:
                    print sel
                print self.z
                self.Plot_bands(select)
                plt.show()
                print 'crashed'
  
    def Cosmology(self,H0=70,Om0=0.25):
    
        mycosmology=FlatLambdaCDM(H0=H0, Om0=Om0)
        self.astropy_cosmo=FlatLambdaCDM(H0= mycosmology.H0, Om0=mycosmology.Om0)
        
    def X0_norm(self):

        instrument = instruments.InstrumentModel("STANDARD")
        B = instrument.EffectiveFilterByBand("B")
        magsys = instruments.MagSys('VEGA')
        zp = magsys.ZeroPoint(B)

        #print 'zeropoint for b-band',zp
        flux_at_10pc = np.power(10., -0.4 * (self.peakAbsMagBesselB-zp))
        fs = salt2.load_filters(['STANDARD::B'])
        mc = salt2.ModelComponents('salt2.npz')
        s = salt2.SALT2([0.], ['STANDARD::B'], mc, fs, z=0.)
        raw_model_norm = s()[0]
        # (10pc)^2 in kpc^2
        #print 'alors man',flux_at_10pc * 1.E-4 / raw_model_norm
        return flux_at_10pc * 1.E-4 / raw_model_norm

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
       
