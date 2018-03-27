import os
import numpy as np
from optparse import OptionParser

parser = OptionParser()
parser.add_option("--paramfile", type="string", default='', help="filter [%default]")
parser.add_option("--simul_name", type="string", default='', help="filter [%default]")
parser.add_option("--min_rf_phase", type="float", default=-20.0, help="filter [%default]")
parser.add_option("--max_rf_phase", type="float", default=60., help="filter [%default]")

opts, args = parser.parse_args()

paramfile=opts.paramfile
simul_name=opts.simul_name
min_rf_phase=-20.
max_rf_phase=60.


params=np.load(paramfile)

print params

outdir='/sps/lsst/data/dev/pgris'
r=[]
torep=None
for val in params:
    filename=outdir
    if val['T0step'] > 0.:
        filename+='/'+val['dirout'].replace('Light_','Fitted_Light_')+'_'+str(val['T0step']).replace('.','_')+'_'+str(min_rf_phase).replace('-','m')+'_'+str(max_rf_phase)
    else:
        filename+='/'+val['dirout'].replace('Light_','Fitted_Light_')+'_random_'+str(min_rf_phase).replace('-','m')+'_'+str(max_rf_phase)
    filename+='/'+val['fieldname']+'/'+str(val['fieldid'])
    filename+='/Season_'+str(val['season'])
    filename+='/z_'+str(val['z'])
    if val['T0step'] > 0:
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
    print('Production ok')
