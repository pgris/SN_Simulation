import numpy as np
from Observations import *
import pylab as plt
from mpl_toolkits.basemap import Basemap
import pickle as pkl
import multiprocessing
import glob

def Get_coadd(filt,fieldid,season,coadd,output_q=None):

    if len(filt) == 0:
        return None

    #print type(filt)
    filt.sort(order='mjd')
    diff_time=0.
    if coadd:
        diff_time=1.


    #print filt['mjd']
    """
    for val in filt:
        print val
    print np.mean(filt['expMJD'][:10])
    """
    #print filt['expMJD']
    
    diff=[io-jo for jo,io in zip(filt['mjd'][:-1], filt['mjd'][1:])]
    
    sel=[i+1 for i in range(len(diff)) if diff[i]>= diff_time]
    sel=[0]+sel+[len(filt)]
    r=[]

    #print sel

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
        #print 'aka', theslice['band'],[po[-1] for po in theslice['band']]
       
        restot.update(res)
        
        
        btot=''
        for po in theslice['band']:
            btot+=po[-1]
        #restot['band']=''.join(sorted(btot))
        
        if season == 6 and fieldid==1427:
            if btot.count('rg') and btot.count('zy'):
                print 'alors?',btot,theslice['mjd']
        
            #btot='rgizy'
        restot['band']='LSSTPG::'+btot

        """
        b=[restot[key] for key in restot.keys() if key.count('band')==0]
        b+=[restot['band']]
        r.append(tuple(b))
        """
        r.append(tuple([restot[key] for key in restot.keys()]))
        dtypes=[]
        for key in restot.keys():
            if key != 'band':
                dtypes.append((key,'f8'))
            else:
                dtypes.append((key,'S20'))
        
        
    #print dtypes,len(dtypes),len(r)
    #print r
    #resu=np.rec.fromrecords(r,names=[key for key in restot.keys()],dtype=['f8']*(len(b)-1)+['S15'])   
    resu=np.rec.array(r,dtype=dtypes)
    #print resu.dtype
    #print resu['visitExpTime'],len(resu)
    if output_q is not None:
        output_q.put({season : resu})

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

def Multi(valb,key,tot_tab,coadd):
    result_queue = multiprocessing.Queue()
    if tot_tab is not None:
        res=tot_tab.copy()
    else:
        res=None
    for season, val in valb.items():
        #print val
        p=multiprocessing.Process(name='Subprocess-'+str(season),target=Get_coadd,args=(val,key,season,coadd,result_queue))
        #restot=Get_coadd(val,key,season)
        p.start()
    
    resultdict = {}
    for season, val in valb.items():
        resultdict.update(result_queue.get())
        

    for p in multiprocessing.active_children():
        p.join()

    for key, vl in resultdict.items():
        res=Concat(res,vl)

    return res

myobs={}

#season=2
fieldname='DD'

conv=dict(zip(['LSSTPG::'+band for band in 'ugrizy'],[1,2,3,4,5,6]))
class_map = lambda x: conv[x]

suffix='_feature_baseline_10yrs'

field_radec=[]
thedir='OpSimLogs'+suffix+'/'+fieldname
files=glob.glob(thedir+'/*_'+fieldname+'_*.txt')
for fi in files:
    splita=fi.split('/')[-1]
    splitb=splita.split('_')
    ra=float(splitb[-2])
    dec=float(splitb[-1].split('.txt')[0])
    field_radec.append((ra,dec))

print 'found',len(field_radec)
for i,val in enumerate(field_radec):
    fieldid_str=str(val[0])+'_'+str(val[1])
    name='Observations_'+fieldname+'_'+fieldid_str+'.txt'
    fieldid=100+i
    obs=Observations(fieldid=fieldid,filename=thedir+'/'+name)
    myobs[fieldid]={}
    for season in range(len(obs.seasons)):
        myobs[fieldid][season]=obs.seasons[season]
    #break
    #if i>500: 
        #break

"""
for fieldid in myobs.keys():
    print fieldid,Cadence(myobs[fieldid])
    #plt.plot(myobs[fieldid]['mjd']*24.*60.,[fieldid]*len(myobs[fieldid]['mjd']))
    #plt.plot(myobs[fieldid]['mjd'],map(class_map,myobs[fieldid]['band']))
    plt.hist(myobs[fieldid]['mjd'][1:]-myobs[fieldid]['mjd'][:-1],histtype='step',bins=35)
"""

bands='ugrizy'

tot_tab=None
tab_all=None
#bands='g'
coadd=False
suffend=''
if coadd:
    suffend='_all'
pkl_file = open('Cadence'+suffix+'_'+fieldname+suffend+'.pkl','wb')
ifile=0
for key,valb in myobs.items():
    print 'hello',key,len(valb)
    ifile+=1
    tot_tab=Multi(valb,key,tot_tab,coadd)
    if ifile > 1 and ifile%100==0:
        print "dumping",ifile
        pkl.dump(tot_tab, pkl_file)
        tot_tab=None

if tot_tab is not None:
    pkl.dump(tot_tab, pkl_file)
    #idx = tot_tab['season']==1.
    #print tot_tab[idx]
        #print tab_all
        #break
    #break
pkl_file.close()
     
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
