import numpy as np
import h5py
from astropy.table import Table,vstack
import sncosmo
import matplotlib.pyplot as plt
from Telescope import *
from astropy import (cosmology, units as u, constants as const)

#bands='ugrizy'
bands='g'

def Read_Fitted(fi):
    tot_data=None
    
    f = h5py.File(fi,'r')
    for i,keyb in enumerate(f.keys()):
        tab=Table.read(fi, path=keyb)
        if tot_data is None:
            tot_data=tab
        else:
            tot_data=vstack([tot_data,tab])

    return tot_data

def Get_Index(list_lc, simu_name, what):

    for i,lca in enumerate(list_lc[simu_name]):
        #print(lca.meta)
        if lca.meta[what] == val[what]:
            return i

    return -1

def Read_LC(fi):
    
    f = h5py.File(fi,'r')
    nlc=len(f.keys())

    r=[]
    for j in range(nlc):
        r.append(Table.read(fi, path='lc_'+str(j)))

    return r

def Plot_LC(lc,sn):

     dust = sncosmo.OD94Dust()
     fitted_model=sncosmo.Model(source='salt2-extended', effects=[dust, dust],
                                effect_names=['host', 'mw'],
                                effect_frames=['rest', 'obs'])
     fitted_model.set(z=sn['z'])
     fitted_model.set(t0=sn['salt2.T0'])
     fitted_model.set(x0=sn['salt2.X0'])
     fitted_model.set(x1=sn['salt2.X1'])
     fitted_model.set(c=sn['salt2.Color']) 
     
     errors={}
     errors['t0']=np.sqrt(sn['salt2.CovT0T0'])
     errors['x0']=np.sqrt(sn['salt2.CovX0X0'])
     errors['x1']=np.sqrt(sn['salt2.CovX1X1'])
     errors['c']=np.sqrt(sn['salt2.CovColorColor'])
     
     print 'Phases : first',(lc['time'][0]-sn['salt2.T0'])/(1.+sn['z']),'last',(lc['time'][-1]-sn['salt2.T0'])/(1.+sn['z'])
     """
     res, fitted_modelb = sncosmo.fit_lc(lc, fitted_model,['t0', 'x0', 'x1', 'c'],bounds={'z':(sn['z']-0.001, sn['z']+0.001)})
     
     print 'ooo',res.errors
     print 'bbb',errors
     """
     sncosmo.plot_lc(lc, model=fitted_model,pulls=True)
     
def Get_Ratios(lca,lcb,t0=0,z=0):

    #print(lca.dtype)
    #print(lca['time'],lcb['time'],lca['band'])
    r=[]
    for band in bands:
        idxa = lca['band']=='LSST::'+band
        sela= lca[idxa]
        idxb = lcb['band']=='LSST::'+band
        selb= lcb[idxb]
        #print(band,sela['time'],selb['time'])
        for i in range(len(sela)):
            
            idxc = np.abs(selb['time']-sela[i]['time'])<1.e-6
            selc=selb[idxc]
            if len(selc) ==1:
                print(i,sela[i]['time'],np.asscalar(selc['time'])-sela[i]['time'])
                #print('hello',len(sela[i]),len(selc),len(sela))
                r.append((band,sela[i]['time'],sela[i]['flux']/np.asscalar(selc['flux']),sela[i]['fluxerr']/np.asscalar(selc['fluxerr']),(sela[i]['time']-t0)/(1.+z)))
        
    return r


fieldname='DD'
fieldid=744
season=2
x1=0.0
color=0.0
z=0.375

simu_name=['sncosmo','snsim']

dirfiles_fitted='/sps/lsst/data/dev/pgris/Fitted_Light_Curves_'
dirfiles_lc='/sps/lsst/data/dev/pgris/Light_Curves_'
prefix='_0_5_m20.0_60.0'

dirname=fieldname+'/'+str(fieldid)+'/Season_'+str(season)+'/z_'+str(z)
fichname=fieldname+'_'+str(fieldid)+'_'+str(z)+'_X1_'+str(x1)+'_C_'+str(color)+'.hdf5'

lc_fitted={}
list_lc={}

telescope=Telescope(airmass=np.median(1.2))

transmission=telescope.throughputs

print(transmission.atmosphere.keys())
for filtre in 'ugrizy':
    band=sncosmo.Bandpass(transmission.atmosphere[filtre].wavelen,transmission.atmosphere[filtre].sb, name='LSST::'+filtre,wave_unit=u.nm)
    sncosmo.registry.register(band)



for simu in simu_name:
    lc_fitted[simu]=Read_Fitted(dirfiles_fitted+simu+prefix+'/'+dirname+'/'+fichname)
    list_lc[simu]=Read_LC(dirfiles_lc+simu+prefix+'/'+dirname+'/'+fichname)

for simu in simu_name:
    print(simu,len(lc_fitted[simu]),len(list_lc[simu]))
 
r=[]
for val in lc_fitted['sncosmo']:
    idx = lc_fitted['snsim']['DayMax']==val['DayMax']
    valb= lc_fitted['snsim'][idx]
    print(val,valb)
    #now try to grab light curves ...
    ia=Get_Index(list_lc,'sncosmo','DayMax')
    ib=Get_Index(list_lc,'snsim','DayMax')
    print(ia,ib,val.dtype)
    
   
    r+=Get_Ratios(list_lc['sncosmo'][ia],list_lc['snsim'][ia],t0=val['DayMax'],z=val['z'])
    #Plot_LC(list_lc['sncosmo'][ia],val)
    #Plot_LC(list_lc['snsim'][ia],valb)

    #break
res=np.rec.fromrecords(r,names=['band','time','flux_ratio','flux_err_ratio','phase'])

#print(res)


for band in bands:
    fig, ax = plt.subplots(1,2)
    idf= res['band']==band
    sel=res[idf]
    """
    ax[0].hist(sel['flux_ratio'],histtype='step',bins=20)
    ax[1].hist(sel['flux_err_ratio'],histtype='step',bins=20)
   
    ax[0].plot(sel['phase'],sel['flux_ratio'],'k.')
    ax[1].plot(sel['phase'],sel['flux_err_ratio'],'k.')
    """
    ax[0].plot(sel['time'],sel['flux_ratio'],'k.')
    ax[1].plot(sel['time'],sel['flux_err_ratio'],'k.')
    

plt.show()
