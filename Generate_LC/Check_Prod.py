import os
import numpy as np

from optparse import OptionParser

parser = OptionParser()

parser.add_option("--simul_name", type="string", default="sncosmo", help="filter [%default]")
parser.add_option("--min_rf_phase", type="float", default=-20.0, help="filter [%default]")
parser.add_option("--max_rf_phase", type="float", default=60., help="filter [%default]")
parser.add_option("--paramfile", type="string", default="", help="filter [%default]")
parser.add_option("--T0random", type="string", default="no", help="filter [%default]")
opts, args = parser.parse_args()

min_rf_phase=opts.min_rf_phase
max_rf_phase=opts.max_rf_phase

simul_name=opts.simul_name
T0random=opts.T0random

params=np.load(opts.paramfile)

print params

outdir='/sps/lsst/data/dev/pgris'
r=[]
torep=None
for val in params:
    filename=outdir
    if T0random == 'no':
        filename+='/'+val['dirout']+'_'+str(val['T0step']).replace('.','_')+'_'+str(min_rf_phase).replace('-','m')+'_'+str(max_rf_phase)
    else:
        filename+='/'+val['dirout']+'_random_'+str(min_rf_phase).replace('-','m')+'_'+str(max_rf_phase)
    filename+='/'+val['fieldname']+'/'+str(val['fieldid'])
    filename+='/Season_'+str(val['season'])
    filename+='/z_'+str(val['z'])
    if T0random == 'no' :
        filename+='/'+val['fieldname']+'_'+str(val['fieldid'])+'_'+str(val['z'])+'_X1_'+str(val['x1'])+'_C_'+str(val['c'])+'.hdf5'
    else:
        filename+='/'+val['fieldname']+'_'+str(val['fieldid'])+'_'+str(val['z'])+'_X1_'+str(val['x1'])+'_C_'+str(val['c'])+'_DayMax_'+str(val['DayMax'])+'.hdf5'
    if not os.path.isfile(filename):
        print('problem here',filename)
        r.append([val[name] for name in params.dtype.names])

if len(r) > 0:
    resu=np.rec.fromrecords(r,names=params.dtype.names)
    
    print(len(params),len(resu))
    np.save('torepro_'+simul_name+'.npy',resu)
else:
    print('Production ok.')
