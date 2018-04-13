import numpy as np
from Observations import *
import pylab as plt
from mpl_toolkits.basemap import Basemap
import pickle as pkl
import multiprocessing
import glob

def Get_coadd(filt,fieldid,season,coadd,j,output_q=None):

    resu=None

    #print('alors',len(filt))
    if len(filt) >= 2:
        
        filt.sort(order='mjd')
    
        diff_time=0.
        if coadd:
            diff_time=1.
    
        diff=[io-jo for jo,io in zip(filt['mjd'][:-1], filt['mjd'][1:])]
        
        sel=[i+1 for i in range(len(diff)) if diff[i]>= diff_time]
        sel=[0]+sel+[len(filt)]
        r=[]

        var_tot=['band','Ra','Dec']
        vars_mean=['mjd','rawSeeing','moon_frac','sky','airmass','m5sigmadepth']
    
        for i in range(len(sel)-1):
            ia=sel[i]
            ib=sel[i+1]
       
            restot=dict(zip([var for var in var_tot],[filt[var][ia] for var in var_tot]))

            theslice=filt[ia:ib]
            theslice.sort(order='mjd')
            res=dict(zip([var for var in vars_mean],[np.median(theslice[var]) for var in vars_mean]))

            res['exptime']=np.sum(filt['exptime'][ia:ib])
        #res['finSeeing']=Get_mean_finSeeing(theslice)
        #res['kAtm']=params.kAtm[filt['band'][ia]]               
            res['Nvisits']=np.sum(theslice['Nexp'])
            if len(theslice)> 1:
                res['cadence']=np.mean(theslice['mjd'][1:]-theslice['mjd'][:-1])
            else:
                res['cadence']=-999.
            res['Difftime']=filt['mjd'][ib-1]-filt['mjd'][ia]
            res['mjd']=filt['mjd'][ia]
            restot['fieldid']=int(fieldid)
            restot['season']=season
            #print('hello season',season)
        #print 'aka', theslice['band'],[po[-1] for po in theslice['band']]
       
            restot.update(res)
      
            btot=''
            for po in theslice['band']:
                btot+=po[-1]
       
            restot['band']='LSSTPG::'+btot
            #print(theslice['band'])
            for bband in 'ugrizy':
                idx=theslice['band']=='LSSTPG::'+bband
                ssel=theslice[idx]
                if len(ssel) > 0:
                    restot['Nvisits_'+bband]=np.sum(ssel['Nexp'])
                    #print('hello',bband,season,np.sum(ssel['Nexp']))
                    if len(ssel)== 1:
                        restot['median_m5sigmadepth_'+bband]=np.asscalar(ssel['m5sigmadepth'])
                        restot['mjd_'+bband]=np.asscalar(ssel['mjd'])
                        
                    else:
                        restot['median_m5sigmadepth_'+bband]=np.median(ssel['m5sigmadepth'])
                        restot['mjd_'+bband]=np.median(ssel['mjd'])
                else:
                    restot['median_m5sigmadepth_'+bband]=-999.
                    restot['mjd_'+bband]=-999.
                    restot['Nvisits_'+bband]=0

                #print('alors',bband,restot['median_m5_'+bband],len(ssel))

            r.append(tuple([restot[key] for key in restot.keys()]))
            dtypes=[]
            for key in restot.keys():
                if key != 'band':
                    if key != 'fieldid':
                        dtypes.append((key,'f8'))
                    else:
                        dtypes.append((key,'i8')) 
                else:
                    dtypes.append((key,'S20'))
        
        
   
        resu=np.rec.array(r,dtype=dtypes)
        #print(resu['Nvisits_g'],np.sum(resu['Nvisits_g']))

    #print resu.dtype
    #print resu['visitExpTime'],len(resu)
    #print(resu)
    if output_q is not None:
        output_q.put({j : resu})

    return resu


def Cadence(obs):

    oob=obs['mjd'][1:]-obs['mjd'][:-1]

    return np.median(oob)

def Concat(tot_tab,res):
    val=None
    if tot_tab is not None:
        val=tot_tab.copy()

    if tot_tab is None:
        val=res
    else:
        if res is not None:
            val=np.concatenate((val,res))
            
    return val

#def Multi(valb,tot_tab,coadd,season,fichname):
def Multi(valb,coadd,season,fichname):
    
    res=None
    pkl_file=open(fichname+'_Season_'+str(season)+'.pkl','wb')

    ivals=range(0,len(valb),10)
    ivals=np.append(ivals,len(valb))

    for i in range(len(ivals)-1):
        
        if i > 1 and i%100==0:
            print("dumping",i,len(res))
            pkl.dump(res, pkl_file)
            res=None

        result_queue = multiprocessing.Queue()
        ia=ivals[i]
        ib=ivals[i+1]

        for j in range(ia,ib):
        #print val
            #print('processing',valb[j][1],valb[j][0])
            p=multiprocessing.Process(name='Subprocess-'+str(season),target=Get_coadd,args=(valb[j][1],valb[j][0],season,coadd,j,result_queue))
            #print('multiproc',j,ia,ib,ivals)
        #restot=Get_coadd(val,key,season)
            p.start()
    
        resultdict = {}

        for j in range(ia,ib):
            resultdict.update(result_queue.get())
        

        for p in multiprocessing.active_children():
            p.join()

        for key, vl in resultdict.items():
            res=Concat(res,vl)

    if res is not None:
        print('Final Dumping',len(res))
        pkl.dump(res, pkl_file)

    pkl_file.close()

def Process(fieldname,fieldids,simu_name,thedir,ifirst,ilast,inum,season,coadd=True):

    myobs={}
    mylist={}
    season_loop=[season]
    if season == -1:
        season_loop=range(10)
    print('season loop',season_loop)
    for seasonl in season_loop:
        mylist[seasonl]=[]
    for i,fieldid in enumerate(fieldids[ifirst:ilast]):
        name='Observations_'+fieldname+'_'+str(fieldid)+'.txt'
        obs=Observations(fieldid=fieldid, filename=thedir+'/'+name,season_length=50.)
        #print('Nseasons',len(obs.seasons))
        myobs[fieldid]={}
        for seasonl in season_loop:
            #myobs[fieldid][season]=obs.seasons[season]
            if seasonl in range(len(obs.seasons)):
                mylist[seasonl].append((fieldid,obs.seasons[seasonl]))
    #print(len(fieldids),ifirst,ilast,len(myobs))

    bands='ugrizy'

    tot_tab=None
    tab_all=None
    coadd=coadd
    suffend=''
    if coadd:
        suffend='_all'
    
    ifile=0
    """
    mylist=[]
    for key,valb in myobs.items():
        print(key,valb.keys())
        if season in valb.keys():
            mylist.append((key,valb[season]))
    """
    #print('hello',len(mylist))
    #inum=1
    name='Cadence'+extent+'_'+fieldname+suffend
    #tot_tab=Multi(mylist,tot_tab,coadd,season,name)
    for seasonl in season_loop:
        print('processing',seasonl)
        Multi(mylist[seasonl],coadd,seasonl,name)
    """
    for key,valb in myobs.items():
        ifile+=1
        tot_tab=Multi(valb,key,tot_tab,coadd,season)
        if ifile > 1 and ifile%100==0:
            print "dumping",ifile
            pkl.dump(tot_tab, pkl_file)
            tot_tab=None

    if tot_tab is not None:
        pkl.dump(tot_tab, pkl_file)
    pkl_file.close()
    """
season=-1
fieldname='WFD'
extent='_feature_baseline_10yrs'
extent='_feature_rolling_half_mask_10yrs'
extent='_feature_rolling_twoThird_10yrs'
extent='_alt_sched'
extent='_alt_sched_rolling'
extent=''
"""
"""
#extent='minion_1016_'
"""
thefile='fieldIDs'+extent+'_'+fieldname+'.txt'

fieldids=np.loadtxt(thefile,dtype={'names': ('name','fieldid'),'formats': ('S8','i4')})
"""
#thedir='OpSimLogs'+extent+'/'+fieldname+'/Year_'+str(season)+'/N_healpix_64'
thedir='/sps/lsst/users/gris/Files_from_OpSim/OpSimLogs'+extent+'/'+fieldname

files=glob.glob(thedir+'/Observations_*.txt')
#files=glob.glob(thedir+'/Observations_*_586.txt')

conv=dict(zip(['LSSTPG::'+band for band in 'ugrizy'],[1,2,3,4,5,6]))
class_map = lambda x: conv[x]
"""
print fieldids
for i,fieldid in enumerate(np.atleast_1d(fieldids)):
    name='Observations_'+fieldname+'_'+str(fieldid[1])+'.txt'
    obs=Observations(fieldid=fieldid[1], filename=thedir+'/'+name)
    myobs[fieldid[1]]={}
    for season in range(len(obs.seasons)):
        myobs[fieldid[1]][season]=obs.seasons[season]
    #break
    #if i>500: 
        #break
"""
fieldids=[]
for fi in files:
    spl=fi.split('/')[-1]
    splb=spl.split('_')[-1]
    fieldids.append(int(splb.split('.')[0]))

#print(fieldids)

ivals=range(0,len(fieldids),5000)
ivals=np.append(ivals,len(fieldids))

#ivals=range(0,10)
#print('hello',ivals)


for iv in range(len(ivals)-1):
    Process(fieldname,fieldids,extent,thedir,ivals[iv],ivals[iv+1],iv+1,season,coadd=True)

     
"""
for band in bands:
idx = val['band'] == 'LSSTPG::'+band
sel= val[idx]

res=Get_coadd(sel,key,season)
tot_tab=Concat(tot_tab,res)
"""
          
        #break

#print tot_tab

"""
pkl_file = open('Cadence_'+fieldname+'.pkl','wb')
pkl.dump(tot_tab, pkl_file)
pkl_file.close()

"""



#plt.show()
