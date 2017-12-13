import numpy as np
from Observations import *
import pylab as plt

thedir='OpSimLogs'
fieldname='WFD'
fieldid=310

name='Observations_'+fieldname+'_'+str(fieldid)+'.txt'
myobs=Observations(fieldid=fieldid, filename=thedir+'/'+fieldname+'/'+name)

print len(myobs.seasons)

for seas in range(1,2):
    obs_tot=myobs.seasons[seas]
    idx = obs_tot['band']=='LSSTPG::i'
    obs=obs_tot[idx]
    diff=obs['mjd'][1:]-obs['mjd'][:-1]
    print np.min(diff),np.max(diff),np.max(obs['mjd'])-np.min(obs['mjd']),np.median(diff)
    plt.plot(obs['mjd'],obs['airmass'],'bo')

plt.show()
