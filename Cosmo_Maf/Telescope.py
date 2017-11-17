from lsst.sims.photUtils import SignalToNoise
from lsst.sims.photUtils import PhotometricParameters
from Throughputs import Throughputs
from lsst.sims.photUtils import Bandpass,Sed
import numpy as np
from astropy.table import vstack,Table
import matplotlib.pyplot as plt
from lsst.sims.photUtils import Bandpass,Sed
from Throughputs import Throughputs
from lsst.sims.photUtils import PhotometricParameters
from Parameters import parameters
import math

class Telescope(object):
    def __init__(self, **kwargs):
        
        self.filters=['u','g','r','i','z','y'] 
        self.Diameter=6.5 #m
        self.DeltaT=30 #s
        self.platescale=0.2 #arsec
        self.gain=2.3
        self.pixel_size=0.2
        self.pixel_area=self.pixel_size**2
        visittime=30.
        self.expTime=dict(zip(self.filters,[30 for i in range(len(self.filters))]))
        self.photParams = PhotometricParameters(nexp=visittime/15.)
       
        self.Cte=3631.*np.pi*6.6*6.5*30/(100.*4*6.626)

        #self.photParams = PhotometricParameters()

        params=['mag_sky','m5','FWHMeff','Tb','Sigmab','zp','counts_zp','Skyb','flux_sky']
        self.data={}
        for par in params:
            self.data[par]={}
        
        self.data['FWHMeff']={'u': 0.92 ,'g': 0.87 ,'r':0.83,'i':0.80,'z':0.78,'y':0.76}
        #self.data['mbsky']={'u': 22.92 ,'g': 22.27 ,'r':21.20,'i':20.47,'z':19.59,'y':18.63}
         
        self.atmos=True
        self.aerosol=True
        self.airmass=1.2
        
        if 'atmos' in kwargs.keys():
            self.atmos=kwargs['atmos']
        
        if 'aerosol' in kwargs.keys():
            self.aerosol=kwargs['aerosol']

        if 'airmass' in kwargs.keys():
            self.airmass=kwargs['airmass']
       

        self.throughputs=Throughputs(atmos=self.atmos,aerosol=self.aerosol)
        self.throughputs.Load_Atmosphere(self.airmass)


        self.Inputs()
        self.Sky()
        self.ZP()
        
        
    @property
    def FWHMeff(self):
        return self.data['FWHMeff']

    @property
    def mag_sky(self):
        return self.data['mag_sky']

    @property
    def m5(self):
        return self.data['m5']

    @property
    def Tb(self):
        return self.data['Tb']

    @property
    def Sigmab(self):
        return self.data['Sigmab']
    @property
    def zp(self):
        return self.data['zp']

    @property
    def ADU_zp(self):
        return self.data['counts_zp']

    @property
    def flux_sky(self):
        return self.data['flux_sky']

    def Inputs(self):
        
        for filtre in self.filters:    
            myup=self.throughputs.darksky.calcInteg(self.throughputs.system[filtre])
            self.data['Tb'][filtre]=self.Calc_Integ(self.throughputs.atmosphere[filtre])
            self.data['Sigmab'][filtre]=self.Calc_Integ(self.throughputs.system[filtre])
            self.data['mag_sky'][filtre]=-2.5*np.log10(myup/(3631.*self.Sigmab[filtre]))

    def Sky(self):

        for filtre in self.filters:
            self.Calc_m5(filtre)

    def Calc_m5(self,filtre):

        filtre_trans=self.throughputs.system[filtre]
        wavelen_min, wavelen_max, wavelen_step=filtre_trans.getWavelenLimits(None,None,None)
                    
        bandpass=Bandpass(wavelen=filtre_trans.wavelen, sb=filtre_trans.sb)

        flatSedb = Sed()
        flatSedb.setFlatSED(wavelen_min, wavelen_max, wavelen_step)
        flux0b=np.power(10.,-0.4*self.mag_sky[filtre])
        flatSedb.multiplyFluxNorm(flux0b)
        if self.atmos:
            self.data['m5'][filtre]=SignalToNoise.calcM5(flatSedb,self.throughputs.atmosphere[filtre],self.throughputs.system[filtre],photParams=self.photParams,FWHMeff=self.FWHMeff[filtre])
            adu_int= flatSedb.calcADU(bandpass=self.throughputs.atmosphere[filtre], photParams=self.photParams)
            self.data['flux_sky'][filtre]=adu_int*self.pixel_area/self.expTime[filtre]/self.gain
        else:
            self.data['m5'][filtre]=SignalToNoise.calcM5(flatSedb,self.throughputs.system[filtre],self.throughputs.system[filtre],photParams=self.photParams,FWHMeff=self.FWHMeff[filtre])
            adu_int= flatSedb.calcADU(bandpass=self.throughputs.system[filtre], photParams=self.photParams)
            self.data['flux_sky'][filtre]=adu_int*self.pixel_area/self.expTime[filtre]/self.gain   


    def ZP(self):
        
        for filtre in self.filters:
            self.ZP_filtre(filtre)
        print 'zeropoints',self.data['zp'],self.data['counts_zp']
        #self.data['zp']=dict(zip(['u','g','r','i','z','y'],[27.03,28.53,28.27,27.91,27.49,26.78]))

    def ZP_filtre(self,filtre):

        self.data['Skyb'][filtre]=self.Cte*np.power(self.Diameter/6.5,2.)*np.power(self.DeltaT/30.,2.)*np.power(self.platescale,2.)*np.power(10.,0.4*(25.-self.mag_sky[filtre]))*self.Sigmab[filtre]
            
        Zb=181.8*np.power(self.Diameter/6.5,2.)*self.Tb[filtre]
        mbZ=25.+2.5*np.log10(Zb) 
        filtre_trans=self.throughputs.system[filtre]
        wavelen_min, wavelen_max, wavelen_step=filtre_trans.getWavelenLimits(None,None,None)
        bandpass=Bandpass(wavelen=filtre_trans.wavelen, sb=filtre_trans.sb)
        flatSed = Sed()
        flatSed.setFlatSED(wavelen_min, wavelen_max, wavelen_step)
        flux0=np.power(10.,-0.4*mbZ)
        flatSed.multiplyFluxNorm(flux0)
        counts = flatSed.calcADU(bandpass, photParams=self.photParams) #number of counts for exptime
        self.data['zp'][filtre]=mbZ
        #print 'hello',counts/self.photParams.exptime
        self.data['counts_zp'][filtre]=counts/self.photParams.exptime

    def Calc_Integ(self,bandpass):
        resu=0.
        dlam=0
        for i,wave in enumerate(bandpass.wavelen):
            if i < len(bandpass.wavelen)-1:
                dlam=bandpass.wavelen[i+1]-wave
                resu+=dlam*bandpass.sb[i]/wave
            #resu+=dlam*bandpass.sb[i]

        return resu
        
    def flux_to_mag(self, flux, band, zp=None):
        if zp is None:
            zp = self.zero_points(band)
        #print 'zp',zp,band
        m = -2.5 * np.log10(flux) + zp
        return m

    def mag_to_flux(self, mag, band, zp=None):
        if zp is None:
            zp = self.zero_points(band)
        return np.power(10., -0.4 * (mag-zp))


    def zero_points(self, band):
        return np.asarray([self.zp[b] for b in band])
