import os
import numpy as np

#field_ids=[290,744,1427,2412,2786]
#field_ids=[1427,2412,2786]
field_ids=[290]

field_names=['DD']*len(field_ids)
fields=dict(zip(field_ids,field_names))

"""
params_sn=[(2.0,-0.2),(-2.0,0.2),(0.0,0.0)] #n_per_batch=10, z=-1
n_per_batch=10
zvals=[-1]
"""

params_sn=[(-999.0,-999.0)]
n_per_batch=2
zvals=[0.01]
zvals+=[z for z in np.arange(0.1,1.,0.1)]


print zvals
for fieldid, fieldname in fields.items():
    for season in range(0,1):
    #for season in [0]:
        for params in params_sn:
            for simu in ['snsim']:
                for z in zvals:
                    cmd='python multiple_batch.py --fieldname '+fieldname+' --fieldid '+str(fieldid)+' --z '+str(round(z,2))+' --stretch '+str(params[0])+' --color '+str(params[1])+' --dirmeas Light_Curves_'+simu+' --season '+str(season)+' --simulator '+simu+' --dirout Fitted_Light_Curves_'+simu+' --n_per_batch '+str(n_per_batch)+' --multiproc no'
                    print cmd
                    os.system(cmd)
