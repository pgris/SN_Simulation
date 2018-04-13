import cPickle as pkl
import glob
from astropy.table import Table, vstack
import matplotlib.pyplot as plt
import numpy as np
from collections import OrderedDict
import h5py

def Get_LC(files):
    r=[]
    for fi in files:
        f = h5py.File(fi,'r')
        for i,key in enumerate(f.keys()):
            lc=Table.read(fi, path=key)
            meta=lc.meta
            r.append((meta['DayMax'],meta['z'],meta['X1'],meta['Color']))
            print('Fitting',i,len(lc))   
            print(lc.meta)
    return np.rec.fromrecords(r,names=['DayMax','z','X1','Color'])
    
def Get_Fitted(files):
    res=None
    for fi in files:
        f = h5py.File(fi,'r')
        print('loading',fi)
        for i,key in enumerate(f.keys()):
            lc=Table.read(fi, path=key)
            if res is None:
                res=lc
            else:
                res=vstack([res,lc])
    return res

fieldname='DD'
fieldid=744
season=0
X1=-999.
Color=-999.
dirmeas='/sps/lsst/data/dev/pgris/Light_Curves_sncosmo_random_m20.0_60.0/'+fieldname+'/'+str(fieldid)+'/Season_'+str(season)
dirmeasb=dirmeas.replace('Light_Curves','Fitted_Light_Curves')
#thefile='Prod_LC/DD_Obs/290/Season_0/DD_290_0.26_0.3_X1_0.0_C_0.0_0.pkl'

#pkl_file = open(thefile,'rb')

#resu=pkl.load(pkl_file)

files = glob.glob(dirmeas+'/*/'+fieldname+'_'+str(fieldid)+'*_X1_*_C_*.hdf5')
filesb = glob.glob(dirmeasb+'/*/'+fieldname+'_'+str(fieldid)+'*_X1_*_C_*.hdf5')

metadata=Get_LC(files)

fitted_LC=Get_Fitted(filesb)


"""
print 'hello',len(metadata),set(metadata['DayMax'])
idx = metadata['DayMax']==metadata['DayMax'][0]
print len(metadata[idx])

idx = metadata['DayMax']==metadata['DayMax'][1]
print len(metadata[idx])
"""


figb, axb = plt.subplots(ncols=2, nrows=2, figsize=(10,9))

axb[0][0].hist(metadata['DayMax'])
axb[0][1].hist(metadata['z'])
axb[1][0].hist(metadata['X1'])
axb[1][1].hist(metadata['Color'])

figc, axc = plt.subplots(ncols=2, nrows=2, figsize=(10,9))

axc[0][0].hist(fitted_LC['DayMax'])
axc[0][1].hist(fitted_LC['z'])
axc[1][0].hist(fitted_LC['X1'])
axc[1][1].hist(fitted_LC['Color'])

plt.show()
