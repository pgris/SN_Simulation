from astropy.table import Table,vstack,Column
import h5py

def read_file(fi):
    tot_data=None
    f = h5py.File(fi,'r')
    print('try loading',fi,len(f.keys()))
    n_lc=0
    for i,keyb in enumerate(f.keys()):
        tab=Table.read(fi, path=keyb)
        n_lc+=len(tab)
        if tot_data is None:
            tot_data=tab
        else:
            tot_data=vstack([tot_data,tab])
        print('number of lc',n_lc)
    return tot_data

fia='/sps/lsst/data/dev/pgris/Fitted_Light_Curves_sncosmo_0_5_m20.0_60.0_covmb/DD/2786/Season_7/z_0.1/DD_2786_0.1_X1_0.0_C_0.0.hdf5'

fib='/sps/lsst/data/dev/pgris/Fitted_Light_Curves_sncosmo_0_5_m20.0_60.0/DD/2786/Season_7/z_0.1/DD_2786_0.1_X1_0.0_C_0.0.hdf5'

tota=read_file(fia)
totb=read_file(fib)

print(len(tota.colnames),len(totb.colnames),tota['salt2.Covmbmb'])
