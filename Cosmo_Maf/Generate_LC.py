from astropy.table import vstack,Table,Column
import cPickle as pkl
#from astropy.io import ascii
import numpy as np
#from Simul_Fit_SN import *
#import matplotlib.pyplot as plt
import sncosmo
from astropy import (cosmology, units as u, constants as const)
from astropy.cosmology import FlatLambdaCDM
from Telescope import *
from lsst.sims.photUtils import Bandpass,Sed
from lsst.sims.photUtils import SignalToNoise
from lsst.sims.photUtils import PhotometricParameters
from lsst.sims.photUtils.EBV import EBVbase
#from saunerie import instruments,salt2
import time
import os
from scipy import interpolate, integrate

class Generate_LC:
    def __init__(self,parameters,fit=False,model='salt2-extended',version='1.0',telescope=None,ra=6.0979440,dec=-1.1051600,airmass=1.2,X0=None,dL=None):

        #time_begin=time.time()

        self.lc=[]
        self.m5={'u':23.61,'g':24.83,'r':24.35,'i':23.88,'z':23.30,'y':22.43}
        self.radeg=np.rad2deg(ra)
        self.decdeg=np.rad2deg(dec)

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
        
        self.sn_type='Ia'

        source=sncosmo.get_source(self.model,version=self.version)
        #self.mycosmology=FlatLambdaCDM(H0=70, Om0=0.25)
        #astropy_cosmo=FlatLambdaCDM(H0= self.mycosmology.H0, Om0=self.mycosmology.Om0)

        dust = sncosmo.OD94Dust()

        #self.transmission=Throughputs(through_dir='FAKE_THROUGH',atmos_dir='FAKE_THROUGH',atmos=False,aerosol=False)
        
        self.transmission=telescope.throughputs
        #self.transmission.Load_Atmosphere(airmass)
        #print 'total elapse time init a',time.time()-time_begin
        self.telescope=telescope

        """
        for band in self.bands:
            themax=np.max(self.transmission.lsst_system[band].sb)
            idx = self.transmission.lsst_system[band].sb >= 0.2*themax
            sel=self.transmission.lsst_system[band].wavelen[idx]
            print band,np.min(sel),np.max(sel)
        """
        """
        self.airmass=airmass
        if self.airmass > 0:
            self.transmission.Load_Atmosphere(airmass)
        """
        #self.lc={}
        #print 'there we go',parameters,len(parameters),parameters.dtype

        self.param=parameters
        """
        SN=sncosmo.Model(source=source,effects=[dust, dust],
                         effect_names=['host', 'mw'],
                         effect_frames=['rest', 'obs'])
        """
        self.SN=sncosmo.Model(source=source)
        
        self.z=self.param['z']
        if X0 is None:
            self.Cosmology()
            lumidist=self.astropy_cosmo.luminosity_distance(self.param['z']).value*1.e3
            self.dL=lumidist
        #X0_snsim = self.X0_norm_snsim() / lumidist** 2
            X0= self.X0_norm()/ lumidist** 2
            #print 'before alpha beta',X0
            alpha=0.13
            beta=3.
            X0 *= np.power(10., 0.4*(alpha*self.param['X1'] -beta*self.param['Color']))
            self.X0=X0

        else:
            self.X0=X0
            self.dL=dL
        #print 'llla',X0,alpha,beta,param['X1'],param['Color'],param['z'],lumidist
        #self.X0=X0
        #print 'hello x0',X0,lumidist
        #SN=sncosmo.Model(source=source)
        self.SN.set(z=self.param['z'])
        self.SN.set(t0=self.param['DayMax'])
        self.SN.set(c=self.param['Color'])
        self.SN.set(x1=self.param['X1'])
        self.SN.set(x0=self.X0)
    
        #print 'total elapse time init b',time.time()-time_begin

        #self.SED={}
        #print 'sncosmo parameters',self.SN.param_names,self.SN.parameters
        lsstmwebv = EBVbase()
        
        ebvofMW = lsstmwebv.calculateEbv(
            equatorialCoordinates=np.array([[np.radians(self.radeg)], [np.radians(self.decdeg)]]))[0]
        #self.SN.set(mwebv=ebvofMW)
        
        #self.SN.set_source_peakabsmag(self.peakAbsMagBesselB, 'bessellB', 'vega',cosmo=self.astropy_cosmo)
        if self.sn_type=='Ia':
            self.mbsim=self.SN._source.peakmag('bessellb','vega')
            #print 'hello',self.mbsim


        #self.metadata=np.rec.fromrecords([(self.X0,self.param['X1'],self.param['Color'],self.param['DayMax'],self.param['z'],self.mbsim)], names = ['X0','X1','Color','DayMax','z','mbsim'])
        
        #self.metadata=Table(rows=[(self.X0,self.param['X1'],self.param['Color'],self.param['DayMax'],self.param['z'],self.mbsim)], names=('X0','X1','Color','DayMax','z','mbsim'), meta={'name': 'metadata'},dtype=tuple(['f8']*6))
        #print 'total elapse time init',time.time()-time_begin

    def __call__(self,obs,out_q=None):

        #print 'hello here',len(mjds),filtre
        #time_begin=time.time()
        
        """
        mjds=obs['mjd']
        airmass=obs['airmass']
        m5=obs['m5sigmadepth']
        filtre=obs['band']
        expTime=obs['exptime']
        """

        fluxes=10.*self.SN.flux(obs['mjd'],self.wave)
        
        wavelength=self.wave/10.
        
        wavelength=np.repeat(wavelength[np.newaxis,:], len(fluxes), 0)
        SED_time = Sed(wavelen=wavelength, flambda=fluxes)
        #print 'total elapse time seds',time.time()-time_begin

        """
        time_begin=time.time()
        for i in range(len(SED_time.wavelen)):
            

            photParams=PhotometricParameters(nexp=expTime[i]/15.)
            sed=Sed(wavelen=SED_time.wavelen[i],flambda=SED_time.flambda[i])
            
            self.transmission.Load_Atmosphere(airmass[i])
            trans=self.transmission.atmosphere[filtre]
            
            flux_SN=sed.calcFlux(bandpass=trans)
            
            if flux_SN >0:
                mag_SN=-2.5 * np.log10(flux_SN / 3631.0)  
                snr_m5_opsim,gamma_opsim=SignalToNoise.calcSNR_m5(mag_SN,trans,m5[i],photParams)
                err_flux_SN=flux_SN/snr_m5_opsim
                e_per_sec = sed.calcADU(bandpass=trans, photParams=photParams) #number of ADU counts for expTime
                    #e_per_sec = sed.calcADU(bandpass=self.transmission.lsst_atmos[filtre], photParams=photParams)
                #print 'alors',e_per_sec,expTime[i],photParams.gain
                e_per_sec/=expTime[i]/photParams.gain
                #print 'alors b',e_per_sec
                #print 'ref',filtre,i,mjds[i],e_per_sec
                    #self.lc[filtre].append(e_per_sec)
                r.append((e_per_sec,mjds[i],flux_SN))
                
                #self.table_for_fit.add_row((mjds[i],flux_SN,err_flux_SN,'LSST::'+filtre,25,'ab'))
                #self.table_for_fit.add_row((mjds[i],mag_SN,mag_SN/snr_m5_opsim,'LSST::'+filtre,25,'ab'))
                self.table_obs.add_row((mjds[i],flux_SN,err_flux_SN,'LSST::'+filtre,2.5*np.log10(3631),'ab',airmass[i],m5[i],expTime[i],e_per_sec,self.telescope.mag_to_flux(m5[i],filtre)))
                #print 'there we go',filtre,e_per_sec,self.telescope.mag_to_flux(m5[i],filtre)
            
        #print 'total elapse time GEN LC',time.time()-time_begin
        """
        #time_begin=time.time()
        fluxes=[]
        transes=[]
        seds=[Sed(wavelen=SED_time.wavelen[i],flambda=SED_time.flambda[i]) for i in range(len(SED_time.wavelen))]
        
        #photParams=[]
        #time_begin=time.time()
        """
        for i in range(len(SED_time.wavelen)):
            
            #photParams=PhotometricParameters(nexp=expTime[i]/15.)
            #sed=Sed(wavelen=SED_time.wavelen[i],flambda=SED_time.flambda[i])
            
            self.transmission.Load_Atmosphere(airmass[i])
            trans=self.transmission.atmosphere[filtre]
            
            #flux_SN=sed.calcFlux(bandpass=trans)
            
            #fluxes.append(flux_SN)
            transes.append(trans)
            #seds.append(sed)
        """
        
        transes=[self.transmission.atmosphere[obs['band'][i][-1]] for i in range(len(SED_time.wavelen))]
        #print 'total elapse time sed',time.time()-time_begin,len(transes)    
        fluxes=[seds[i].calcFlux(bandpass=transes[i]) for i in range(len(SED_time.wavelen))]
        #print 'total elapse time sed',time.time()-time_begin
        #print 'before',len(fluxes),len(seds),len(transes),len(m5),len(expTime),len(mjds),self.param['X1'],self.param['Color']
        #print fluxes,mjds
        #time_begin=time.time()

        table_obs=Table(obs)
        snr_m5=[]
        e_per_sec_list=[]
        for i in range(len(SED_time.wavelen)):
            photParams=PhotometricParameters(nexp=table_obs['exptime'][i]/15.)
            flux_SN=fluxes[i]
            if flux_SN > 0:
                trans=self.transmission.atmosphere[table_obs['band'][i][-1]]
                mag_SN=-2.5 * np.log10(flux_SN / 3631.0)  
                snr_m5_opsim,gamma_opsim=SignalToNoise.calcSNR_m5(mag_SN,trans,table_obs['m5sigmadepth'][i],photParams)
                err_flux_SN=flux_SN/snr_m5_opsim
                e_per_sec = seds[i].calcADU(bandpass=trans, photParams=photParams) #number of ADU counts for expTime
                    #e_per_sec = sed.calcADU(bandpass=self.transmission.lsst_atmos[filtre], photParams=photParams)
                #print 'alors',e_per_sec,expTime[i],photParams.gain
                e_per_sec/=table_obs['exptime'][i]/photParams.gain
                #print('hello',mag_SN,flux_SN,e_per_sec,snr_m5_opsim)
                snr_m5.append(snr_m5_opsim)
                e_per_sec_list.append(e_per_sec)
                #print fluxes,mags
            else:
                snr_m5.append(1)
                e_per_sec_list.append(1)

        #print('passed')
        
        table_obs.add_column(Column(fluxes, name='flux'))
        table_obs.add_column(Column(snr_m5, name='snr_m5'))
        table_obs.add_column(Column(e_per_sec_list, name='flux_e'))
        idx = table_obs['flux'] >= 0.
        table_obs=table_obs[idx]
        """
        mags=-2.5 * np.log10(table_obs['flux'] / 3631.0)
        def snr(bands,mags,sky,seeing,expTime,ron):
            pixel_scale=0.2
            mag_sky=dict(zip('ugrizy',[22.95,22.24,21.20,20.47,19.60,18.63]))
            FWHMeff = {'u':0.92, 'g':0.87, 'r':0.83, 'i':0.80, 'z':0.78, 'y':0.76}
            neff=np.array(2.266*(seeing/pixel_scale)**2)
            flux_e=self.telescope.mag_to_flux(mags,[band[-1] for band in bands])*expTime
            #flux_sky=self.telescope.mag_to_flux([mag_sky[band[-1]] for band in bands],[band[-1] for band in bands])*expTime
            flux_sky=self.telescope.mag_to_flux(sky,[band[-1] for band in bands])*expTime
            #res=flux_e/np.sqrt((flux_sky*pixel_scale**2+ron**2)*neff)
            res=flux_e/np.sqrt(flux_e+(flux_sky*pixel_scale**2+ron**2)*neff)
            #for i in range(len(mags)):
            #print mags[i],flux_e[i],res[i]
            return res,flux_e/expTime

        snrs,flux_e=snr(table_obs['band'],mags,table_obs['sky'],table_obs['seeing'],table_obs['exptime'],5.)
        print('ici',flux_e,snrs)
        """
        table_obs.add_column(Column(table_obs['flux']/table_obs['snr_m5'], name='fluxerr'))
        #table_obs.add_column(Column(flux_e, name='flux_e'))
        #table_obs.add_column(Column(flux_e/snrs, name='flux_e_err'))
        table_obs.add_column(Column(table_obs['flux_e']/table_obs['snr_m5'], name='flux_e_err'))
        table_obs.add_column(Column([2.5*np.log10(3631)]*len(table_obs),name='zp'))
        table_obs.add_column(Column(['ab']*len(table_obs),name='zpsys'))
        #table_obs.add_column(Column(self.telescope.mag_to_flux(mags,[band[-1] for band in table_obs['band']]),name='flux_e_sec'))
        #table_obs.add_column(Column(self.telescope.mag_to_flux(obs['m5sigmadepth'],[band[-1] for band in table_obs['band']]),'flux_5sigma_e_sec'))
        table_obs['band'][:]=np.array([b.replace('LSSTPG','LSST') for b in table_obs['band']])

        table_obs.rename_column('mjd','time')
        table_obs.rename_column('m5sigmadepth','m5')
         
        for colname in ['exptime','rawSeeing','seeing','moon_frac','sky','kAtm','airmass','m5','Nexp','Ra','Dec']:
            if colname in table_obs.colnames:
                table_obs.remove_column(colname)

        #print(table_obs[['flux_e','snr_m5']])
        """
        table_obs.remove_column('exptime')
        table_obs.remove_column('rawSeeing')
        table_obs.remove_column('seeing')
        table_obs.remove_column('moon_frac')
        table_obs.remove_column('sky')
        table_obs.remove_column('kAtm')
        table_obs.remove_column('airmass')
        table_obs.remove_column('m5')
        table_obs.remove_column('Nexp')
        table_obs.remove_column('Ra')
        table_obs.remove_column('Dec')
        """
        #table_obs.remove_column('flux_e_sec')
        #table_obs.remove_column('flux_5sigma_e_sec')
    

        #print len(table_obs)
        #print table_obs
        #print len(table_obs),table_obs.dtype
        #print np.array(np.sort(table_obs,order='band'))
        """
        idx=fluxes > 0.
        fluxes=fluxes[idx]
        seds=np.array(seds)[idx]
        transes=np.array(transes)[idx]
        m5=m5[idx]
        expTime=expTime[idx]
        photParams=[PhotometricParameters(nexp=expTime[i]/15.) for i in range(len(expTime))]
        airmass=airmass[idx]
        mjds=mjds[idx]
        filtre=filtre[idx]

        #print 'allors',len(fluxes),len(seds),len(transes),len(m5),len(expTime),len(photParams),len(mjds)
        #print fluxes,mjds
        mags=-2.5 * np.log10(fluxes / 3631.0)  
        gamma = np.asarray([SignalToNoise.calcGamma(transes[i],m5[i],photParams[i]) for i in range(len(mags))])
        x=np.power(10.,0.4*(mags-m5))
        snr_m5_gamma=1./np.asarray(np.sqrt((0.04-gamma)*x+gamma*(x**2)))
        #snr_m5_gamma_orig=[SignalToNoise.calcSNR_m5(mags[i],transes[i],m5[i],photParams[i]) for i in range(len(mags))]
        #print 'hhh',snr_m5_gamma,snr_m5_gamma_orig
        #snr_m5_opsim=[SignalToNoise.calcSNR_m5(mags[i],transes[i],m5[i],photParams[i]) for i in range(len(mags))]
        err_fluxes=fluxes/snr_m5_gamma
        #err_fluxes_orig=[fluxes[i]/snr_m5_gamma_orig[i][0] for i in range(len(fluxes))]
        
        e_per_sec = [seds[i].calcADU(bandpass=transes[i], photParams=photParams[i]) for i in range(len(transes))] #number of ADU counts for expTime
        
        e_per_sec=[e_per_sec[i]/(expTime[i]/photParams[i].gain) for i in range(len(e_per_sec))]
        #mags_new=[self.telescope.flux_to_mag(e_per_sec[i],filtre[i][-1])[0] for i in range(len(e_per_sec))]

        #print 10**(-0.4*(mags-mags_new))
        print filtre[:]
        print 'test',e_per_sec,self.telescope.mag_to_flux(mags,[band[-1] for band in filtre]),e_per_sec/self.telescope.mag_to_flux(mags,[band[-1] for band in filtre])

        

        #snrs=snr(filtre,mags,obs['sky'],np.array([0.8]*len(mags)),obs['exptime'],5.)
        print snr_m5_gamma,snrs,np.median(snr_m5_gamma/snrs),np.mean(snr_m5_gamma/snrs),np.min(snr_m5_gamma/snrs),np.max(snr_m5_gamma/snrs)

        table_obs=Table()
        table_obs.add_column(Column(mjds, name='time'))
        table_obs.add_column(Column(fluxes, name='flux'))
        table_obs.add_column(Column(err_fluxes, name='fluxerr'))
        table_obs.add_column(Column(['LSST::'+filtre[i][-1] for i in range(len(filtre))], name='band'))
        table_obs.add_column(Column([2.5*np.log10(3631)]*len(mjds),name='zp'))
        table_obs.add_column(Column(['ab']*len(mjds),name='zpsys'))
        table_obs.add_column(Column(airmass,name='airmass'))
        table_obs.add_column(Column(m5,name='m5'))
        table_obs.add_column(Column(expTime,name='expTime'))
        table_obs.add_column(Column(e_per_sec,name='flux_e_sec'))
        table_obs.add_column(Column(self.telescope.mag_to_flux(m5,[filtre[i][-1] for i in range(len(filtre))]),'flux_5sigma_e_sec'))
        #print 'total elapse time GEN LC b',time.time()-time_begin
        #print table_obs
        """
        if out_q is not None:
            out_q.put({filtre[0][-1] : table_obs}) 
 
        return table_obs
        #plt.show()
        
  
    def Cosmology(self,H0=70,Om0=0.25):
    
        mycosmology=FlatLambdaCDM(H0=H0, Om0=Om0)
        self.astropy_cosmo=FlatLambdaCDM(H0= mycosmology.H0, Om0=mycosmology.Om0)
        
    def X0_norm_snsim(self):

        instrument = instruments.InstrumentModel("STANDARD")
        B = instrument.EffectiveFilterByBand("B")
        magsys = instruments.MagSys('VEGA')
        zp = magsys.ZeroPoint(B)

        #print 'snsim zp',zp
        #print 'zeropoint for b-band',zp
        flux_at_10pc = np.power(10., -0.4 * (self.peakAbsMagBesselB-zp))
        fs = salt2.load_filters(['STANDARD::B'])
        mc = salt2.ModelComponents('salt2.npz')
        s = salt2.SALT2([0.], ['STANDARD::B'], mc, fs, z=0.)
        raw_model_norm = s()[0]
        #print 'flux snsim',raw_model_norm
        # (10pc)^2 in kpc^2
        #print 'alors man',flux_at_10pc * 1.E-4 / raw_model_norm
        return flux_at_10pc * 1.E-4 / raw_model_norm

    def X0_norm(self):

        name='STANDARD'
        band='B'
        thedir='../Cosmo_Maf'
        #thedir='.'

        os.environ[name] = thedir+'/Instruments/Landolt'

        trans_standard=Throughputs(through_dir='STANDARD',telescope_files=[],filter_files=['sb_-41A.dat'],atmos=False,aerosol=False,filterlist=('A'),wave_min=3559,wave_max=5559)
     
        mag, spectrum_file =self.Get_Mag(thedir+'/MagSys/VegaBD17-2008-11-28.dat',name,band)

        #print 'alors mag',mag, spectrum_file
        sed=Sed()

        sed.readSED_fnu(filename=thedir+'/'+spectrum_file)
        CLIGHT_A_s  = 2.99792458e18         # [A/s]
        HPLANCK = 6.62606896e-27
        sedb=Sed(wavelen=sed.wavelen,flambda=sed.wavelen*sed.fnu/(CLIGHT_A_s * HPLANCK))
        #print 'alors man',sed.wavelen*sed.fnu/(CLIGHT_A_s * HPLANCK)
        flux=self.calcInteg(bandpass=trans_standard.system['A'],signal=sedb.flambda,wavelen=sedb.wavelen)

        zp=2.5*np.log10(flux)+mag
        flux_at_10pc = np.power(10., -0.4 * (self.peakAbsMagBesselB-zp))
        
        #print 'zp',zp,flux_at_10pc

        source=sncosmo.get_source(self.model,version=self.version)
        SN=sncosmo.Model(source=source)

        SN.set(z=0.)
        SN.set(t0=0)
        SN.set(c=self.param['Color'])
        SN.set(x1=self.param['X1'])
        SN.set(x0=1)

        fluxes=10.*SN.flux(0.,self.wave)
        
        wavelength=self.wave/10.
        SED_time = Sed(wavelen=wavelength, flambda=fluxes)

        expTime=30.
        photParams = PhotometricParameters(nexp=expTime/15.)
        trans=Bandpass(wavelen=trans_standard.system['A'].wavelen/10., sb=trans_standard.system['A'].sb)
        e_per_sec = SED_time.calcADU(bandpass=trans, photParams=photParams) #number of ADU counts for expTime
                    #e_per_sec = sed.calcADU(bandpass=self.transmission.lsst_atmos[filtre], photParams=photParams)
        e_per_sec/=expTime/photParams.gain*photParams.effarea
        #print 'hello',e_per_sec
        """
        SN.set(c=self.param['Color'])
        SN.set(x1=self.param['X1'])
        """

        
        #print 'My zp',zp,flux
        return flux_at_10pc * 1.E-4 /e_per_sec

    def Get_Mag(self,filename,name,band):
        
        sfile=open(filename,'rb')
        spectrum_file='unknown'
        for line in sfile.readlines():
            if 'SPECTRUM' in line:
                spectrum_file=line.split(' ')[1].strip()
            if name in line and band in line:
                return float(line.split(' ')[2]),spectrum_file

        sfile.close()


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
       

    def calcInteg(self, bandpass, signal,wavelen):
        
        #print 'oyer',len(signal),len(wavelen),len(bandpass.sb),np.min(wavelen),np.max(wavelen),np.min(bandpass.wavelen),np.max(bandpass.wavelen)
        

        fa = interpolate.interp1d(bandpass.wavelen,bandpass.sb)
        fb = interpolate.interp1d(wavelen,signal)
        
        min_wave=np.max([np.min(bandpass.wavelen),np.min(wavelen)])
        max_wave=np.min([np.max(bandpass.wavelen),np.max(wavelen)])
        
        #print 'before integrand',min_wave,max_wave
        
        wavelength_integration_step=5
        waves=np.arange(min_wave,max_wave,wavelength_integration_step)
        
        integrand=fa(waves) *fb(waves) 
        #print 'rr',len(f2(wavelen)),len(wavelen),len(integrand)
        
        range_inf=min_wave
        range_sup=max_wave
        n_steps = int((range_sup-range_inf) / wavelength_integration_step)

        x = np.core.function_base.linspace(range_inf, range_sup, n_steps)
        #print len(waves),len(x)
        return integrate.simps(integrand,x=waves) 
