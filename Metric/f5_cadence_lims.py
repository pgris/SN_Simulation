from Generate_LC import Generate_LC
from Telescope import *

from croaks import NTuple
from pycosmo.cosmo import CosmoLambda
from saunerie import salt2, snsim, psf

model_components = None
cosmo = CosmoLambda()

def init_lcmodel(bands, filename='salt2.npz'):
    """
    Utility function to load a SALT2 light curve model. 
    The model components are cached. 

    This should be a function of the snsim module.
    """
    global model_components
    if model_components is None:
        print 'we have to reload model_components'
        if filename is None:
            model_components = salt2.ModelComponents.load_from_saltpath()
        else:
            model_components = salt2.ModelComponents(filename)
    fs = salt2.load_filters(np.unique(bands))
    lcmodel = snsim.SALT2LcModel(fs, model_components)
    return lcmodel

def get_pulse_shapes(bands, 
                     restframe_phase_range = (-20., 40.),
                     z=1.1, X1=-2., Color=0.2, DayMax=0., 
                     filename='salt2.npz',
                     cosmo=cosmo, plot=False):
    """
    Call the SALT model for the fiducial SN specified in argument, on
    a grid of regularly spaced mjd (observer frame), corresponding to
    ``restframe phase range``
    
    .. note : slow because no 
    """
    pmin, pmax = restframe_phase_range
    mjd_min = np.floor(pmin * (1.+z) + DayMax)
    mjd_max = np.ceil(pmax * (1.+z) + DayMax)
    mjd = np.arange(mjd_min, mjd_max, 1.)
    
    sn = np.rec.fromrecords([(z, cosmo.dL(z), X1, Color, DayMax)], 
                            names=['z', 'dL', 'X1', 'Color', 'DayMax'])    
    
    ret = {}
    lcmodel = init_lcmodel(bands, filename=filename)
    for bn in bands:
        b = [bn] * len(mjd)
        #        lcmodel = init_lcmodel(b)
        ret[bn] = lcmodel(sn, mjd, b)
        
    return mjd, ret


def f5_cadence_lims(SNR=dict(zip(['LSSTPG::' + b for b in "rizy"], 
                                 [25., 25., 25., 15.])),
                    zs=[0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1],
                    X1=-2., Color=0.2):
    #                    bands=["LSSTPG::" + b for b in "grizy"]):
    """
    like the functin above, but just returns a dict 
    with the cadence limits for the requested SNR
    """
    bands = SNR.keys()
    lims = {}
    for z in zs:
        mjd, shapes = get_pulse_shapes(z=z, X1=X1, Color=Color, bands=bands)
        lims[z] = {}
        for b in bands:
            Li2 = np.sqrt(np.sum(shapes[b]**2))
            lim = 5. * Li2 / SNR[b]
            print 'there we go',z,b,lim
            lims[z][b] = lim
    return lims

def f5_cadence_lims_sncosmo(SNR=dict(zip([b for b in "griz"], 
                                 [30., 40., 30., 20.])),
                    zs=[0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1],
                    X1=-2., Color=0.2,
                            expTime=dict(zip([b for b in "griz"], 
                                               [30., 30., 30., 30.])),telescope=None):

    bands=SNR.keys()
    restframe_phase_range = (-20., 40.)	
    pmin, pmax = restframe_phase_range
    r=[]
    lims={}
    for z in zs:
        r.append((z,-2.0,0.2,0.))

    #params=np.rec.fromrecords(r, names = ['z', 'X1', 'Color', 'DayMax'])
    params=np.rec.array(r, dtype=[('z', 'f8'),('X1', 'f8'), ('Color', 'f8'),('DayMax','f8')])
    #print 'before',params['z']
    cadence=1.
    
    #m5_lim={'u':23.61,'g':24.83,'r':24.35,'i':23.88,'z':23.30,'y':22.43}

    #this m5_lim values correspond to an exposure time of 30s...needs to be corrected for DDF

    obs_param=parameters()
    m5_lim=obs_param.m5
    for key, vals in expTime.items():
        m5_lim[key]=m5_lim[key]+1.25*np.log10(vals/30.)

    for param in params:
        rb=[]
        z=param['z']
        lims[z]={}
        mjd_min = np.floor(pmin * (1.+z) + param['DayMax'])
        mjd_max = np.ceil(pmax * (1.+z) + param['DayMax'])
        mjd = np.arange(mjd_min, mjd_max, cadence)
        lc_sncosmo=Generate_LC(param,telescope=telescope)
        #lc_sncosmo=Generate_LC(param,mjd,expTime).lc
        #airmass_val=np.repeat(airmass, len(mjd), 0)
       
        for mjd_val in mjd:
            for band in bands:
                rb.append((mjd_val,band,obs_param.seeing[band],m5_lim[band],expTime[band],obs_param.msky[band],))
        table_obs=np.rec.fromrecords(rb,names=['mjd','band','seeing','m5sigmadepth','exptime','sky'])
        lc=lc_sncosmo(table_obs)
       
        print('redshift',z)
        for filtre in bands:
           
            #res=resultdict[filtre][0]
            idx=lc['band']==filtre
            res=lc[idx]
            #print('bou',filtre,res['flux_e']/res['flux_e_err'])
            #idx = res['flux_e']/res['flux_e_err'] > 5.
            #res=res[idx]
            Li2=np.sqrt(np.sum(np.square(res['flux_e'])))
            lim = 5. * Li2 / SNR[filtre]
            #print(filtre,lim,Li2,SNR[filtre])
            lims[z][filtre]=lim

    return lims 

atmos=True
instr_snsim=psf.find('LSSTPG')
instr=Telescope(atmos=atmos,aerosol=False,airmass=1.2)
zs=[0.1, 0.2, 0.3, 0.4, 0.5]
SNR=dict(zip(['LSSTPG::' + b for b in "griz"],
             [30., 40., 30., 20.]))
SNR_sncosmo=dict(zip([b for b in "griz"],
                     [30., 40., 30., 20.]))

bands=SNR.keys()

lims_snsim = f5_cadence_lims(zs=zs, SNR=SNR)
lims_sncosmo = f5_cadence_lims_sncosmo(zs=zs, SNR=SNR_sncosmo,telescope=instr)

print('types',type(lims_snsim),type(lims_sncosmo))
np.save('f5_cadence_lims_snsim.npy',lims_snsim)
np.save('f5_cadence_lims_sncosmo.npy',lims_sncosmo)
