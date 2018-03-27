import h5py
from astropy.table import Table

thefile='/sps/lsst/data/dev/pgris/Light_Curves_snsim_0_2_b/DD/290/Season_0/DD_290_0.05_X1_-999.0_C_-999.0_190_200.hdf5'
thefileb='/sps/lsst/data/dev/pgris/Light_Curves_sncosmo_0_2_b/DD/290/Season_0/DD_290_0.05_X1_-999.0_C_-999.0_190_200.hdf5'

f = h5py.File(thefile,'r')
fb = h5py.File(thefileb,'r')

print len(f.keys()),len(fb.keys())

for key in f.keys():
    tab=Table.read(thefile, path=key)
    tabb=Table.read(thefileb, path=key)
    print key,len(tab),len(tabb),len(tab.meta),len(tabb.meta)
    
