import numpy as np
import cPickle as pkl
from astropy.table import vstack
import matplotlib.pyplot as plt
import h5py
import glob
from astropy.table import Table,vstack,Column

data=None
for season in range(10):
	filename='Sel_Selb_10_2/Sel_DD_744_Season_'+str(season)+'.pkl'

	loaded=pkl.load(open(filename,'rb'))
	if data is None:
		data=loaded
	else:
		data=vstack([data,loaded])

""" this is for JLA
names=['name','z','zhel','dz','mbfit','salt2.Covmbmb','salt2.X1','salt2.CovX1X1','salt2.Color','salt2.CovColorColor','3rdvar','d3rdvar','tmax','dtmax','salt2.CovX1mb','salt2.CovColormb','salt2.CovColorX1','set','ra','dec','biascor']
jla=np.loadtxt('jla_lcparams.txt',dtype={'names': tuple(names),'formats': tuple(['S15']+[np.float]*(len(names)-1))})
jla=np.sort(jla,order='z')
plt.hist(jla['salt2.CovColorColor'],histtype='step',color='k')
"""

#load hdf5 data
tot_data=None
fieldname='DD'
fieldid=744
season=3
dirmeas='/sps/lsst/data/dev/pgris/Fitted_Light_Curves_sncosmo_0_5_m20.0_60.0/'+fieldname+'/'+str(fieldid)+'/Season_'+str(season)
what=dirmeas+'/z_*/'+fieldname+'_'+str(fieldid)+'_*_X1_0.0_C_0.0.hdf5'
print('thter',what)
files= glob.glob(what)
for fi in files:
            f = h5py.File(fi,'r')
            print('try loading',fi,len(f.keys()))
            n_lc=0
            for i,keyb in enumerate(f.keys()):
                tab=Table.read(fi, path=keyb)
                if tot_data is None:
                    tot_data=tab
                else:
                    tot_data=vstack([tot_data,tab])

fig, axes = plt.subplots(figsize=(10,9))
axes.hist(np.sqrt(data['salt2.CovColorColor']),histtype='step',color='r',range=[0.,0.2])
axes.hist(np.sqrt(tot_data['salt2.CovColorColor']),histtype='step',color='k',range=[0.,0.2])


figb, axesb = plt.subplots(figsize=(10,9))
axesb.hist(np.sqrt(data['salt2.CovX1X1']),histtype='step',color='r',range=[0.,1.])
axesb.hist(np.sqrt(tot_data['salt2.CovX1X1']),histtype='step',color='k',range=[0.,1.])

figc, axesc = plt.subplots(figsize=(10,9))
axesc.plot(data['X1'],data['Color'],'k*')
axesc.plot(tot_data['X1'],tot_data['Color'],'ro')

figd, axesd = plt.subplots(figsize=(10,9))
#axesd.hist(data['X1'],histtype='step',color='r',range=[-2.,2.],bins=20)
axesd.plot(data['X1'],np.sqrt(data['salt2.CovX1X1']),'ko')

fige, axese = plt.subplots(figsize=(10,9))
axese.hist(data['Color'],histtype='step',color='r',range=[-0.5,0.5])

plt.show()
