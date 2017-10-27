import cPickle as pkl
import glob
from astropy.table import Table, vstack
import matplotlib.pyplot as plt
import numpy as np
from collections import OrderedDict

fieldname='DD'
fieldid=290
season=0
X1=-999.
Color=-999.
dirmeas='Prod_LC/'+fieldname+'/'+str(fieldid)+'/Season_'+str(season)

#thefile='Prod_LC/DD_Obs/290/Season_0/DD_290_0.26_0.3_X1_0.0_C_0.0_0.pkl'

#pkl_file = open(thefile,'rb')

#resu=pkl.load(pkl_file)

files = glob.glob(dirmeas+'/'+fieldname+'_'+str(fieldid)+'*_X1_'+str(X1)+'_C_'+str(Color)+'*.pkl')


metadata=None
r=[]
for fi in files:
    pkl_file = open(fi,'rb')
    resu=pkl.load(pkl_file)
    
    #print 'go',len(resu)
    for val in resu:
        ro=[]
        names=[]
        #print val.meta
        for i, (key, value) in enumerate(val.meta.iteritems()):
            ro.append(value[0])
            names.append(key)
            #print 'hello',value[0]
        r.append(tuple(ro))
        #print r,names
    

metadata=np.rec.fromrecords(r,names=names)


figb, axb = plt.subplots(ncols=2, nrows=2, figsize=(10,9))

axb[0][0].hist(metadata['DayMax'])
axb[0][1].hist(metadata['X1'])
axb[1][1].hist(metadata['Color'])

plt.show()
