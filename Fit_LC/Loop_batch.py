import os
import numpy as np

#field_ids=[290,744,1427,2412,2786]
#field_ids=[1427,2412,2786]
field_ids=[290]

field_names=['DD']*len(field_ids)
fields=dict(zip(field_ids,field_names))
zstep=0.025

"""
#params_sn=[(2.0,-0.2),(-2.0,0.2),(0.0,0.0)] #n_per_batch=10, z=-1
params_sn=[(0.0,0.0)]
n_per_batch=2
zvals=[-1]

"""
params_sn=[(-999.0,-999.0)]
#params_sn=[(-1.2,-0.06),(-1.2,0.08),(1.0,-0.04)]
n_per_batch=10
zvals=[0.01]
zvals+=[z for z in np.arange(zstep,1.1+zstep,zstep)]
#zvals=[0.8]
#zvals=np.arange(0.1,0.21,0.01)

T0step=[0.3,0.5,1.0]
T0step=[0.5]
#X1Cval='_b'
X1Cval=''
config_pro=[]
for val in T0step:
    config_pro.append('sncosmo_'+str(val).replace('.','_')+X1Cval)
#X1Cval='' 

print zvals
for fieldid, fieldname in fields.items():
    for season in range(0,1):
    #for season in [0]:
        for params in params_sn:
            #for simu in ['sncosmo_test_'+str(T0step).replace('.','_')+X1Cval,'snsim_test']:
            for simu in config_pro:
                for z in zvals:
                    cmd='python multiple_batch.py --fieldname '+fieldname+' --fieldid '+str(fieldid)+' --z '+str(z)+' --stretch '+str(params[0])+' --color '+str(params[1])+' --dirmeas Light_Curves_'+simu+' --season '+str(season)+' --simulator '+simu+' --dirout Fitted_Light_Curves_'+simu+' --n_per_batch '+str(n_per_batch)+' --multiproc yes'
                    print cmd
                    os.system(cmd)
