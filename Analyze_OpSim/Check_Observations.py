import numpy as np
from Observations import *
import pylab as plt
import heapq
import pickle as pkl

thedir='/sps/lsst/users/gris/Files_from_OpSim/OpSimLogs'
fieldname='WFD'
fieldid=586

name='Observations_'+fieldname+'_'+str(fieldid)+'.txt'
myobs=Observations(fieldid=fieldid, filename=thedir+'/'+fieldname+'/'+name)

"""
dir_pkl='/sps/lsst/data/dev/pgris/sims_operations/SN_LC/Obs_minion_1016'

name_pkl='Observations_'+fieldname+'_'+str(fieldid)+'.pkl'

pkl_file = open(dir_pkl+'/'+name_pkl,'rb')
thedict=pkl.load(pkl_file)
obser_orig=thedict['dataSlice']
"""
#print len(myobs.seasons),obser_orig.dtype

for seas in range(10):
    obs_tot=myobs.seasons[seas]
    for band in 'g':
        idx = obs_tot['band']=='LSSTPG::'+band
        obs=obs_tot[idx]
        obs=obs_tot
        diff=obs['mjd'][1:]-obs['mjd'][:-1]
        largest=heapq.nlargest(2,diff)
        idxa = diff == largest[0]
        idxb = diff == largest[1]
        print band,np.min(diff),np.max(diff),np.max(obs['mjd'])-np.min(obs['mjd']),np.median(diff),np.argmax(diff),obs['mjd'][np.argmax(diff)],heapq.nlargest(2,diff),obs['mjd'][idxa],obs['mjd'][idxb]
        #idx = obser_orig['filter']==band
        #sel = obser_orig[idx]
        plt.scatter(obs['mjd'],obs['airmass'],s=80,facecolors='none', edgecolors='k')
        #plt.plot(sel['expMJD'],sel['airmass'],'r*')

plt.show()
