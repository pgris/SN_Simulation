import h5py
from astropy.table import Table

thefile='/sps/lsst/data/dev/pgris/Fitted_test/DD/290/Season_0/DD_290_0.1_X1_2.0_C_-0.2_0_684.hdf5'

f = h5py.File(thefile,'r')

print f.keys(),len(f.keys())

for i,key in enumerate(f.keys()):
    tab=Table.read(thefile, path=key)
    print tab
    if i > 10:
        break
