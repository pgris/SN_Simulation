import numpy as np
from Observations import *
from SN_Rate import *
import matplotlib.pyplot as plt
from Load_X1_C import Load_X1_Color

def Load_File(filename):
        res=np.loadtxt(filename, dtype={'names': ('x1', 'c', 'weight_x1','weight_c','weight_tot'),'formats': ('f8', 'f8', 'f8','f8','f8')})
        res.sort(order='weight_tot')
        res[:]=res[::-1]
    
        return res

def Get_duration(fieldname,fieldid,season,OpSim_Logs_dir):

    filename=OpSim_Logs_dir+'/Observations_'+fieldname+'_'+str(fieldid)+'.txt'
    myobs=Observations(fieldid=fieldid, filename=filename)
    myseason=myobs.seasons[season]
    iddx=myseason['band']!='LSSTPG::u'
    mysel=myseason[iddx]

    min_season=np.min(mysel['mjd'])
    max_season=np.max(mysel['mjd'])

    return min_season,max_season

def Get_params(fieldname,fieldid,season,Opsimlog_dir,modif_x1c,dirout):

    min_season,max_season=Get_duration(fieldname,fieldid,season,Opsimlog_dir)
    duration=max_season-min_season
    print(duration)
    rate_name='Perrett'
    sn_rate=SN_Rate(rate=rate_name,duration=duration/365.25)
    zz,rate,err_rate,nsn,err_nsn=sn_rate(bins=np.arange(0.01,1.5,0.0001))
#print(zz,nsn,err_nsn,np.sum(nsn))
    r=[]
    for i in range(int(np.sum(nsn))):
        z=np.random.choice(zz,1,p=nsn/np.sum(nsn))[0]
    #zval.append(z)
        T0 = np.random.uniform(min_season,max_season)
        ztype='high_z'
        if z < 0.1:
            ztype='low_z'
        x1=np.random.choice(modif_x1c[ztype]['x1'],1,p=modif_x1c[ztype]['weight_x1']/np.sum(modif_x1c[ztype]['weight_x1']))[0]
        color=np.random.choice(modif_x1c[ztype]['c'],1,p=modif_x1c[ztype]['weight_c']/np.sum(modif_x1c[ztype]['weight_c']))[0]
        r.append((fieldname,fieldid,season,z,-1.,x1,color,1.,1.,fieldname,dirout,'Ia',Opsimlog_dir,round(T0,3)))
    return r

fieldname='DD'
Opsimlog_dir='/sps/lsst/users/gris/Files_from_OpSim/OpSimLogs/'+fieldname

grid_x1_c=Load_X1_Color(-999.0,-999.0,1).tab

modif_x1c={}
for ztype in ['low_z','high_z']:
    r=[]
    sel_x1c=grid_x1_c[ztype]
    for val in sel_x1c:
        r.append((np.asscalar(val['x1']),np.asscalar(val['weight_x1']),np.asscalar(val['c']),np.asscalar(val['weight_c'])))
    #print(r)
    modif_x1c[ztype]=np.rec.fromrecords(r,names=['x1','weight_x1','c','weight_c'])
for ztype in ['low_z','high_z']:
    modif_x1c[ztype]=Load_File('Dist_X1_Color_jla_'+ztype+'.txt')

simul_name='sncosmo'
dirout='Light_Curves_'+simul_name
r=[]
fieldid=744
for season in range(1):
    r+=Get_params(fieldname,fieldid,season,Opsimlog_dir,modif_x1c,dirout)

print(len(r))

params=np.rec.fromrecords(r,names=['fieldname','fieldid','season','z','T0step','x1','c','weight_x1','weight_c','dirmeas','dirout','sntype','Opsimlog','DayMax'])

print(params['z'])
np.save('prod_lc_'+simul_name+'_random_'+str(fieldid)+'_test.npy',params)


figb, axb = plt.subplots(2, 2, figsize=(10,9))
idx = params['z'] < 0.1
idxb = params['z'] >= 0.1
axb[0][0].hist(params[idx]['x1'],histtype='step',color='k')
axb[0][1].hist(params[idx]['c'],histtype='step',color='k')
axb[1][0].hist(params[idxb]['x1'],histtype='step',color='k')
axb[1][1].hist(params[idxb]['c'],histtype='step',color='k')
plt.show()


"""
figb, axb = plt.subplots(ncols=2, nrows=2, figsize=(10,9))

axb[0][0].hist(params['z'],histtype='step')
axb[0][1].hist(params['DayMax'],histtype='step')
axb[1][0].hist(params['x1'],histtype='step')
axb[1][1].hist(params['c'],histtype='step')

plt.show()
"""
