import numpy as np
from Observations import *
import pylab as plt
from mpl_toolkits.basemap import Basemap
import pickle as pkl

def Get_coadd(filt,fieldid,season):

    if len(filt) == 0:
        return None

    #print type(filt)
    filt.sort(order='mjd')

    #print filt['mjd']
    """
    for val in filt:
        print val
    print np.mean(filt['expMJD'][:10])
    """
    #print filt['expMJD']
    
    diff=[io-jo for jo,io in zip(filt['mjd'][:-1], filt['mjd'][1:])]
    
    sel=[i+1 for i in range(len(diff)) if diff[i]>= 1.]
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

        res=dict(zip([var for var in vars_mean],[np.mean(theslice[var]) for var in vars_mean]))

        res['exptime']=np.sum(filt['exptime'][ia:ib])
        #res['finSeeing']=Get_mean_finSeeing(theslice)
        #res['kAtm']=params.kAtm[filt['band'][ia]]               
        res['Nvisits']=len(theslice)
        res['cadence']=np.median(theslice['mjd'][1:]-theslice['mjd'][:-1])
        res['Difftime']=filt['mjd'][ib-1]-filt['mjd'][ia]
        res['mjd']=filt['mjd'][ia]
        restot['fieldid']=fieldid
        restot['season']=season
        restot.update(res)
       
        r.append(tuple([restot[key] for key in restot.keys()]))
       
    resu=np.rec.fromrecords(r,names=[key for key in restot.keys()])   

    #print resu['visitExpTime'],len(resu)
    
    return resu


def Cadence(obs):

    oob=obs['mjd'][1:]-obs['mjd'][:-1]

    return np.median(oob)

thefile='fieldIDs_minion_1016_WFD.txt'

fieldids=np.loadtxt(thefile,dtype={'names': ('name','fieldid'),'formats': ('S8','i4')})

myobs={}

season=2
fieldname='WFD'
thedir='OpSimLogs/'+fieldname

conv=dict(zip(['LSSTPG::'+band for band in 'ugrizy'],[1,2,3,4,5,6]))
class_map = lambda x: conv[x]

print fieldids,len(fieldids)
for i,fieldid in enumerate(fieldids):
    name='Observations_'+fieldname+'_'+str(fieldid[1])+'.txt'
    obs=Observations(fieldid=fieldid[1], filename=thedir+'/'+name)
    myobs[fieldid[1]]={}
    for season in range(len(obs.seasons)):
        myobs[fieldid[1]][season]=obs.seasons[season]

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


for key,valb in myobs.items():
    for season, val in valb.items():
        for band in bands:
            idx = val['band'] == 'LSSTPG::'+band
            sel= val[idx]
        #print len(sel)
            res=Get_coadd(sel,key,season)
    #print res[['mjd','Difftime','cadence','Nvisits']],np.median(res['mjd'][1:]-res['mjd'][:-1]),np.median(res['Nvisits'])
            if tot_tab is None:
                tot_tab=res
            else:
                if res is not None:
                    tot_tab=np.concatenate((tot_tab,res))

#print tot_tab

pkl_file = open('Cadence_'+fieldname+'.pkl','wb')
pkl.dump(tot_tab, pkl_file)
pkl_file.close()





plt.show()
