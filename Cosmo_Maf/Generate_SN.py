from astropy.table import vstack,Table
import cPickle as pkl
from astropy.io import ascii
import numpy as np
#from Simul_Fit_SN import *
import matplotlib.pyplot as plt
import sncosmo
from astropy import (cosmology, units as u, constants as const)
from astropy.cosmology import FlatLambdaCDM
from Throughputs import Throughputs
from lsst.sims.photUtils import Bandpass,Sed
from lsst.sims.photUtils import SignalToNoise
from lsst.sims.photUtils import PhotometricParameters
from saunerie import instruments,salt2

class Generate_SN:
    def __init__(self,logobs,parameters,model,version):


        self.mjd_name='mjd'
        self.FWHMeff='seeing'
        peakAbsMagBesselB=-19.0906
        alpha=0.13
        beta=3.
        
        mycosmology=FlatLambdaCDM(H0=70, Om0=0.25)
        self.astropy_cosmo=FlatLambdaCDM(H0= mycosmology.H0, Om0=mycosmology.Om0)
        instrument = instruments.InstrumentModel("STANDARD")
        B = instrument.EffectiveFilterByBand("B")
        magsys = instruments.MagSys('VEGA')
        zp = magsys.ZeroPoint(B)

        print 'zeropoint for b-band',zp
        flux_at_10pc = np.power(10., -0.4 * (peakAbsMagBesselB-zp))
        fs = salt2.load_filters(['STANDARD::B'])
        mc = salt2.ModelComponents('salt2.npz')
        s = salt2.SALT2([0.], ['STANDARD::B'], mc, fs, z=0.)
        raw_model_norm = s()[0]
        # (10pc)^2 in kpc^2
        self.X0_norm = flux_at_10pc * 1.E-4 / raw_model_norm

        table_obs=ascii.read(logobs,fast_reader=False)

        colnames=[]
        sfile=open(logobs,'r')
        for line in sfile.readlines():
            if line.count('#'):
                colnames.append(line[1:].split(':')[0].strip())

        for i,val in enumerate(table_obs.colnames):
            #if val.count('#') >= 1:
            table_obs.rename_column(val, colnames[i])

        self.lc=[]

        print table_obs

        self.params=Table(names=('t0','c','x1','z','ra','dec','status','fit','sn_type','sn_model','sn_version','mbsim','x0','dL'),dtype=('f8','f8','f8','f8','f8','S8','S8','f8','S8','S8','S8','f8','f8','f8'))

        self.transmission=Throughputs()
       
        dust = sncosmo.OD94Dust()

        model=model
        version=version

        if model == 'salt2-extended':
            model_min=300.
            model_max=180000.
            wave_min=3000.
            wave_max=11501.

        if model=='salt2':
            model_min=2000.
            model_max=9200.0
            wave_min=model_min
            wave_max=model_max

        wave= np.arange(wave_min,wave_max,1.)
        
        sn_type='Ia'
        source=sncosmo.get_source(model,version=version)
        ra_field=0.
        dec_field=0.
        
        table_obs.sort(self.mjd_name)
       
        for iv,param in enumerate(parameters):
            #mysn=Simul_Fit_SN(param['DayMax'],param['Color'],param['X1'],param['z'],table_obs,ra=param['ra'],dec=param['dec'])
            table_LC=Table(names=('idSN','filter','expMJD','visitExpTime','FWHMeff','moon_frac','filtSkyBrightness','kAtm','airmass','fiveSigmaDepth','Nexp','e_per_sec','e_per_sec_err'), dtype=('i8','S7','f8', 'f8','f8','f8', 'f8','f8','i8','f8','f8','f8','f8'))


            lumidist=self.astropy_cosmo.luminosity_distance(param['z']).value*1.e3
       
            X0 = self.X0_norm / lumidist** 2
            X0 *= np.power(10., 0.4*(alpha*param['X1'] -beta*param['Color']))

            print 'X0 val',X0

            SN=sncosmo.Model(source=source,effects=[dust, dust],
                              effect_names=['host', 'mw'],
                              effect_frames=['rest', 'obs'])
            SN.set(z=param['z'])
            SN.set(t0=param['DayMax'])
            SN.set(c=param['Color'])
            SN.set(x1=param['X1'])
            SN.set(x0=X0)
            
            #SN.set_source_peakabsmag(peakAbsMagBesselB, 'bessellB', 'vega',cosmo=mycosmology)
            fluxes=10.*SN.flux(table_obs[self.mjd_name],wave)
            wavelength=wave/10. 
            wavelength=np.repeat(wavelength[np.newaxis,:], len(fluxes), 0)
            SED_time = Sed(wavelen=wavelength, flambda=fluxes)
           
            self.params.add_row((param['DayMax'],param['Color'],param['X1'],param['z'],ra_field,dec_field,'unkown',None,sn_type,model,version,SN._source.peakmag('bessellb','vega'),SN.get('x0'),lumidist))
            
            for i in range(len(SED_time.wavelen)):
                obs=table_obs[i]
                sed=Sed(wavelen=SED_time.wavelen[i],flambda=SED_time.flambda[i])
                visittime=obs['exptime']
                filtre=obs['band'][-1]
                photParams = PhotometricParameters(nexp=visittime/15.)
                e_per_sec = sed.calcADU(bandpass=self.transmission.lsst_atmos_aerosol[filtre], photParams=photParams) #number of ADU counts for expTime
                e_per_sec/=visittime/photParams.gain
                
                flux_SN=sed.calcFlux(bandpass=self.transmission.lsst_atmos_aerosol[filtre])

                FWHMeff=obs[self.FWHMeff]
                if flux_SN >0:
                    #print 'hello flux',flux_SN,fluxes[i]
                    mag_SN=-2.5 * np.log10(flux_SN / 3631.0)
                    m5_calc,snr_m5_through=self.Get_m5(filtre,mag_SN,obs['sky'],photParams,FWHMeff)
                    
                    table_LC.add_row((iv,'LSST::'+filtre,obs[self.mjd_name],visittime,FWHMeff,obs['moon_frac'],obs['sky'],obs['kAtm'],obs['airmass'],obs['m5sigmadepth'],obs['Nexp'],e_per_sec,e_per_sec/snr_m5_through))

            self.lc.append(table_LC)
            #self.params.append(outdict)
            #self.lc.append(mysn.table_LC)
            #self.params.append(mysn.outdict)
            #break

    def Get_m5(self,filtre,mag_SN,msky,photParams,FWHMeff):
        
        wavelen_min, wavelen_max, wavelen_step=self.transmission.lsst_system[filtre].getWavelenLimits(None,None,None)
                
        flatSed = Sed()
        flatSed.setFlatSED(wavelen_min, wavelen_max, wavelen_step)
        flux0=np.power(10.,-0.4*msky)
        flatSed.multiplyFluxNorm(flux0)
        m5_calc=SignalToNoise.calcM5(flatSed,self.transmission.lsst_atmos_aerosol[filtre],self.transmission.lsst_system[filtre],photParams=photParams,FWHMeff=FWHMeff)
        snr_m5_through,gamma_through=SignalToNoise.calcSNR_m5(mag_SN,self.transmission.lsst_atmos_aerosol[filtre],m5_calc,photParams)
        
        return m5_calc,snr_m5_through
