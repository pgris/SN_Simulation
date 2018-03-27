import numpy as np
import glob
from Observations import *
import pylab as plt

#simuname='_feature_baseline_10yrs'
simuname='_alt_sched_rolling'
#thedir='/sps/lsst/users/gris/Files_from_OpSim/OpSimLogs'
thedir='/pbs/throng/lsst/users/gris/SN_Simulation/Analyze_OpSim/OpSimLogs'

thedir+=simuname

files=glob.glob(thedir+'/WFD/Year_2/N_healpix_64_coadd/*.txt')

#print(files)


r=[]
num=-1
for fi in files:
    spl=fi.split('/')[-1]
    fieldid=spl.split('_')[-1].split('.')[0]
    
    myobs=Observations(fieldid=int(fieldid), filename=fi,season_length=365.)
    if len(myobs.seasons) >= 1:
        season=myobs.seasons[0]
        if len(season) >= 100:
            #print fieldid,len(season)
            num+=1
            rcad=[]
            for band in 'ugrizy':
                idx = season['band']=='LSSTPG::'+band
                sel=season[idx]
                diff=sel['mjd'][1:]-sel['mjd'][:-1]
                if len(sel) >= 2:
                    #print fieldid,band,len(season),np.median(diff)
                    #r.append((fieldid,band,np.median(diff)))
                    rcad.append(np.median(diff))
                else:
                    rcad.append(0.)
            #print(rcad,len(rcad))
            if len(rcad) != 6:
                print 'Pb hre',len(rcad)
            rtot=[fieldid,band]
            rtot+=rcad
            r.append(tuple(rtot))

            #if num > 100:
             #   break

names=['fieldid','band']+['cadence_'+band for band in 'ugrizy']
res=np.rec.fromrecords(r,names=names)
print 'hello',len(res),res.dtype

idx=True
for band in 'rizy':
    idx&= (res['cadence_'+band]<3.1)&(res['cadence_'+band]>2.9)

sel=res[idx]
print len(sel),sel[['fieldid','cadence_u','cadence_g','cadence_r','cadence_z','cadence_y']]

for fieldid in np.unique(res['fieldid']):
    idd=res['fieldid']==fieldid
    sell=res[idd]
    #print([np.median(sell['cadence_'+band]) for band in 'ugrizy'])
    """
    for band in 'ugrizy':
        ixb=sell['band']==band
        selb=sell[ixb]
        figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(7,5))
        figa.suptitle(band + ' band')
        axa.hist(selb['cadence'],bins=20)
        axa.set_xlabel('Cadence [day-1]')
        axa.set_ylabel('Number of Entries')
    break
    """    
plt.show()

                             

