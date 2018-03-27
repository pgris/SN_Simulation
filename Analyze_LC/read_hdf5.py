import h5py
from astropy.table import Table
import numpy as np
thedir='/sps/lsst/data/dev/pgris/Fitted_Light_Curves_sncosmo_0_5_m20.0_60.0'

z=1.075
X1=2.0
Color=-0.2

thefile=thedir+'/DD/290/Season_1/z_'+str(z)+'/DD_290_'+str(z)+'_X1_'+str(X1)+'_C_'+str(Color)+'.hdf5'

f = h5py.File(thefile,'r')

print f.keys(),len(f.keys())
thresh=10
for i,key in enumerate(f.keys()):
    tab=Table.read(thefile, path=key)
    for val in tab:
        print(val['N_Phase_m5'],val['N_Phase_p30'],val['N_nights_m20_p_p30'],val['Near_peak_gap'],val['N_g_snrmax_'+str(thresh)]+val['N_r_snrmax_'+str(thresh)]+val['N_i_snrmax_'+str(thresh)]+val['N_z_snrmax_'+str(thresh)]+val['N_y_snrmax_'+str(thresh)],val['fit_status'])
    if i > 10:
        break
