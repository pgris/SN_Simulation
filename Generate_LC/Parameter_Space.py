import numpy as np


simul_name='sncosmo'

r=[]
weight_x1=1.
weight_c=1.
#extent='_feature_baseline_10yrs'
extent='_feature_baseline_10yrs'
extent_b='Year_2/N_healpix_64_coadd'
extent=''
#dirout='Light_Curves_sncosmo'+extent+'/'+extent_b
#fieldnames=['WFD']
dirout='Light_Curves_'+simul_name
fieldnames=['DD']
Opsimlog_dir='/sps/lsst/users/gris/Files_from_OpSim/OpSimLogs'+extent+'/'+fieldnames[0]
#Opsimlog_dir='/sps/lsst/users/gris/Files_from_OpSim/OpSimLogs/'+fieldnames[0]
fieldids=[290,744,1427,2412,2786]
#fieldids=[100,101,102,103]
#fieldids=[101]
#fieldids=[6084]
seasons=range(0,10)
#seasons=[0]
T0steps=[0.5]
x1_c_vals= [(0.0,0.0),(2.0,-0.2),(-2.0,0.2)]
zmax=1.4
zmin=0.0
zmin_0=0.01
zstep=0.025
zvals=[zmin_0]
zvals+=[zz for zz in np.arange(zmin+zstep,zmax+2.*zstep,zstep)]

#print zvals

for fieldname in fieldnames:
    for fieldid in fieldids:
        for season in seasons:
            for T0step in T0steps:
                for (x1,c) in x1_c_vals:
                    for z in zvals:
                        r.append((fieldname,fieldid,season,z,T0step,x1,c,weight_x1,weight_c,fieldname,dirout,'Ia',Opsimlog_dir,-1.))

params=np.rec.fromrecords(r,names=['fieldname','fieldid','season','z','T0step','x1','c','weight_x1','weight_c','dirmeas','dirout','sntype','Opsimlog','DayMax'])

print params

np.save('prod_lc_'+simul_name+extent+'.npy',params)
