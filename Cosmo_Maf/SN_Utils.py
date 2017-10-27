from scipy.interpolate import InterpolatedUnivariateSpline as Spline1d
from scipy import interpolate, integrate
import numpy as np

class SN_Utils:

    def __init__(self):
        self.splB=None

    def Load(self):

         #from F. Mondon 2017/10/20
        # wavelength limits for salt2 model
        wl_min_sal = 3000
        wl_max_sal = 7000
        thedir='../Cosmo_Maf'
    
    #interpolation of TB and Trest
        filt2 = np.genfromtxt(thedir+'/snfit_data/Instruments/SNLS3-Landolt-model/sb-shifted.dat')
        filt2=np.genfromtxt(thedir+'/Instruments/Landolt/sb_-41A.dat')
        wlen = filt2[:,0]
        tran = filt2[:,1]
        self.splB = Spline1d(wlen, tran, k=1,ext = 1)

    #interpolation of ref spectrum
        data = np.genfromtxt(thedir+'/snfit_data/MagSys/bd_17d4708_stisnic_002.ascii')
        data = np.genfromtxt(thedir+'/MagSys/bd_17d4708_stisnic_002.ascii')
        dispersion = data[:,0]
        flux_density = data[:,1]
        self.splref = Spline1d(dispersion, flux_density, k=1,ext = 1)

  
    #interpolation of the spectrum model
        template_0 = np.genfromtxt(thedir+'/snfit_data/salt2-4/salt2_template_0.dat')    
        template_1 = np.genfromtxt(thedir+'/snfit_data/salt2-4/salt2_template_1.dat')
#    salt2source=sncosmo.SALT2Source('/users/divers/lsst/mondon/hubblefit/sncosmo_jla/salt2-4')
    
        
        wlM0 = []
        M0 = []
        for i in range(len(template_0[:,0])):
            if template_0[:,0][i] == 0.0:
                wlM0.append(template_0[:,1][i]) 
                M0.append(template_0[:,2][i])
        self.splM0 = Spline1d(wlM0, M0, k=1,ext = 1)

        wlM1 = []
        M1 = []
        for i in range(len(template_1[:,0])):
            if template_1[:,0][i] == 0.0:
                wlM1.append(template_1[:,1][i]) 
                M1.append(template_1[:,2][i])
        self.splM1 = Spline1d(wlM1, M1, k=1,ext = 1)

    #computation of the integral
        dt = 100000
        self.xs = np.linspace(float(wl_min_sal), float(wl_max_sal), dt)
        self.dxs = (float(wl_max_sal-wl_min_sal)/(dt-1))


    def mB(self,params):
    
        if self.splB is None:
            self.Load()

#    I1=np.sum((splM0(xs)*10**-12+res.parameters[3]*splM1(xs)*10**-12)*(10**(-0.4*salt2source.colorlaw(xs)*res.parameters[4]))*xs*splB(xs)*dxs)
        I1 = np.sum((self.splM0(self.xs)*10**-12+params['x1']*self.splM1(self.xs)*10**-12)*(10**(-0.4*self.color_law_salt2(self.xs)*params['c']))*self.xs*self.splB(self.xs)*self.dxs)    
        I2 = np.sum(self.splref(self.xs)*self.xs*self.splB(self.xs)*self.dxs)
        #print I1, I2   
    
    #computation of mb
        mref = 9.907
        mb = -2.5*np.log10(params['x0']*(I1/I2))+mref
        
        return mb
    
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

    def Get_Mag(self,filename,name,band):
        
        sfile=open(filename,'rb')
        spectrum_file='unknown'
        for line in sfile.readlines():
            if 'SPECTRUM' in line:
                spectrum_file=line.split(' ')[1].strip()
            if name in line and band in line:
                return float(line.split(' ')[2]),spectrum_file

        sfile.close()

    def color_law_salt2(self,wl):
        B_wl = 4302.57
        V_wl = 5428.55
        l = (wl-B_wl)/(V_wl-B_wl)
        l_lo = (2800.-B_wl)/(V_wl-B_wl)
        l_hi = (7000.-B_wl)/(V_wl-B_wl)
        a = -0.504294
        b = 0.787691
        c = -0.461715
        d = 0.0815619
        cst = 1-(a+b+c+d)
        cl = []
        for i in range (len(l)):
            if l[i] > l_hi:
                cl.append(-(cst*l_hi+l_hi**2*a+l_hi**3*b+l_hi**4*c+l_hi**5*d+(cst+2*l_hi*a+3*l_hi**2*b+4*l_hi**3*c+5*l_hi**4*d)*(l[i]-l_hi)))
            if l[i] < l_lo:
                cl.append(-(cst*l_lo+l_lo**2*a+l_lo**3*b+l_lo**4*c+l_lo**5*d+(cst+2*l_lo*a+3*l_lo**2*b+4*l_lo**3*c+5*l_lo**4*d)*(l[i]-l_lo))) 
            if l[i]>= l_lo and l[i]<= l_hi:
                cl.append(-(cst*l[i]+l[i]**2*a+l[i]**3*b+l[i]**4*c+l[i]**5*d)) 
        return np.array(cl)

    def Covar(self,params,covar,vparam_names):

        res={}
        h=1.e-6
        Der=np.zeros(shape=(len(vparam_names),1))

        #print params
        par_var=params.copy()
        ider=-1
        for i,key in enumerate(vparam_names):
            par_var[key]+=h
            ider+=1
            Der[ider]=(self.mB(par_var)-self.mB(params))/h
            par_var[key]-=h

        Prod=np.dot(covar,Der)

        for i,key in enumerate(vparam_names):
            if key != 'c':
                res['salt2.Cov'+key.upper()+'mb']=Prod[i,0]
            else:
               res['salt2.CovColormb']=Prod[i,0] 

        res['salt2.Covmbmb']=np.asscalar(np.dot(Der.T,Prod))

        return res
        

    def Test(self):


        print 'in Test'
        """
        Salt2Model
        BEGIN_OF_FITPARAMS Salt2Model
        DayMax 53690.0336018 0.105513809169 
        Redshift 0.1178 0 F 
        Color -0.0664131339433 0.0234330339301 
        X0 0.00030732251016 8.89813428854e-06 
        X1 -0.0208012409076 0.160846457522 
        CovColorColor 0.00054910707917 -1 
        CovColorDayMax 0.00040528682468 -1 
        CovColorX0 -1.68238293879e-07 -1 
        CovColorX1 0.00114702847231 -1 
        CovDayMaxDayMax 0.0111331639253 -1 
        CovDayMaxX0 -2.94345317778e-07 -1 
        CovDayMaxX1 0.0131008809199 -1 
        CovX0X0 7.91767938168e-11 -1 
        CovX0X1 -7.23852420336e-07 -1 
        CovX1X1 0.0258715828973 
        """

        salt2_res={}
        salt2_res['DayMax']=53690.0336018
        salt2_res['Color']=-0.0664131339433
        salt2_res['X0']=0.00030732251016
        salt2_res['X1']=-0.0208012409076
        salt2_res['CovColorColor']=0.00054910707917
        salt2_res['CovColorDayMax']=0.00040528682468
        salt2_res['CovColorX0']=-1.68238293879e-07
        salt2_res['CovColorX1']=0.00114702847231
        salt2_res['CovDayMaxDayMax']=0.0111331639253
        salt2_res['CovDayMaxX0']=-2.94345317778e-07
        salt2_res['CovDayMaxX1']=0.0131008809199
        salt2_res['CovX0X0']=7.91767938168e-11
        salt2_res['CovX0X1']=-7.23852420336e-07
        salt2_res['CovX1X1']=0.0258715828973 
        #salt2_res['']=
        vparam_names=['t0','c','x0','x1']
        covar=np.zeros(shape=(len(vparam_names),len(vparam_names)))

        covar[0,1]=salt2_res['CovColorDayMax']
        covar[0,2]=salt2_res['CovDayMaxX0']
        covar[0,3]=salt2_res['CovDayMaxX1']
        
        covar[1,2]=salt2_res['CovColorX0']
        covar[1,3]=salt2_res['CovColorX1']
        
        covar[2,3]=salt2_res['CovX0X1']
        
        covar=covar+covar.T

        covar[0,0]=salt2_res['CovDayMaxDayMax']
        covar[1,1]=salt2_res['CovColorColor']
        covar[2,2]=salt2_res['CovX0X0']
        covar[3,3]=salt2_res['CovX1X1']

        
        #print covar
        
        
        
        params={}
        params['t0']=salt2_res['DayMax']
        params['c']=salt2_res['Color']
        params['x0']=salt2_res['X0']
        params['x1']=salt2_res['X1']

        self.Covar(params,covar,vparam_names)
