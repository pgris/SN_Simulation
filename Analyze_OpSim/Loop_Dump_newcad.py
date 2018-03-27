import os
import glob
import multiprocessing
import numpy as np
import time
from optparse import OptionParser

def Process(fieldname,fieldid,season,simu_name,coadd,time_coadd,j,output_q):
    cmd='python Dump_OpSim_in_File.py --fieldname '+fieldname+' --fieldid '+str(fieldid)+' --season '+str(season)+' --simu '+simu_name+' --coadd '+coadd+' --time_coadd '+str(time_coadd)
    #print(cmd)
    os.system(cmd)
 
    if output_q is not None:
        output_q.put({j : j})

def Multiproc(n_per_batch,tab,coadd,time_coadd,target,args):
    n_per_batch=n_per_batch

    inter=range(0,len(tab),n_per_batch)
    inter=np.append(inter,len(tab))

    restot=[]
    time_begin=time.time()
#print range(len(tab)-1)
    for jo in range(len(inter)-1):
        ida=inter[jo]
        idb=inter[jo+1]
        result_queue = multiprocessing.Queue()

        if (jo%10) == 0:
            print('Processing',jo,time.time()-time_begin)

        for j in range(ida,idb):
        #print val
            myargs=[tab[var][j] for var in args]
            myargs+=[coadd,time_coadd]
            myargs+=[j,result_queue]
            
            p=multiprocessing.Process(name='Subprocess-'+str(j),target=target,args=tuple(myargs))
        #restot=Get_coadd(val,key,season)
            p.start()
 
        resultdict = {}
    
        for j in range(ida,idb):
            resultdict.update(result_queue.get())

        #print('joining')
        for p in multiprocessing.active_children():
            p.join()

    """
        for j in range(ida,idb):
            restot+=resultdict[j]

    return restot
    """

parser = OptionParser()
parser.add_option("--fieldname", type="string", default="WFD", help="filter [%default]")
parser.add_option("--season", type="int", default=2, help="filter [%default]")
parser.add_option("--simu", type="string", default="feature_baseline_10yrs", help="filter [%default]")
parser.add_option("--coadd", type="string", default="no", help="filter [%default]")
parser.add_option("--time_coadd", type="float", default=24.*3600., help="filter [%default]")

opts, args = parser.parse_args()

#suffix='feature_baseline_10yrs'
simu_name=opts.simu
#suffix='feature_rolling_half_mask_10yrs'
#suffix='feature_rolling_twoThird_10yrs'
#suffix='alt_sched'
#suffix='alt_sched_rolling'
season=opts.season
fieldname=opts.fieldname
coadd=opts.coadd
time_coadd=opts.time_coadd # 1 day
#field_radec=[(9.45,-44.),(35.708333,-4.75),(53.125,-28.1),(150.1,2.18194444444)]


#for season in range(1,11):
field_radec=[]
fields=[]
thedir='/sps/lsst/data/dev/pgris/OpSimFiles_'+simu_name+'/Season_'+str(season)+'/N_healpix_64'
files=glob.glob(thedir+'/Observations_*.npy')
print(thedir+'/Observations*_'+fieldname+'_*.npy')

r=[]
for fi in files:
    splita=fi.split('/')[-1]
    splitb=splita.split('_')
    #ra=float(splitb[-2])
    fieldid=int(splitb[-1].split('.npy')[0])
    #field_radec.append((ra,dec))
    fields.append(fieldid)
    #print('Number of Fields',len(field_radec))
    r.append((fieldname,fieldid,season,simu_name))

fieldids=np.rec.fromrecords(r,names=['fieldname','fieldid','season','simu_name'])

Multiproc(8,fieldids,coadd,time_coadd,Process,['fieldname','fieldid','season','simu_name'])
"""
for val in fields:
    #cmd='python Dump_OpSim_in_File.py --fieldRA '+str(val[0])+' --fieldDec '+str(val[1])+' --fieldname '+fieldname+' --simu '+suffix
    cmd='python Dump_OpSim_in_File.py --fieldname WFD --fieldid '+str(val)+' --season '+str(season)+' --simu '+suffix
    print(cmd)
    os.system(cmd)
"""        
