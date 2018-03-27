import numpy as np
from Observations import *
import pylab as plt
#import pickle as pkl
import glob
import multiprocessing
import time
from mpl_toolkits.basemap import Basemap
from optparse import OptionParser

parallels = np.arange(-90.,90,30.)
meridians = np.arange(0.,360.,60.)

def Plot_Cadence(tab,what,season):
    
    """
    radec=np.unique((tab_all[['Ra','Dec']]))
    
    norm=[]
    
    for val in radec:
        idd = (tab['Ra']==val[0])&(tab['Dec']==val[1])
        sel=tab[idd]
        norm.append(np.median(sel['median_cadence']))
       

    lons = np.rad2deg([radec[i][0] for i in range(len(radec))])
    lats = np.rad2deg([radec[i][1] for i in range(len(radec))])
    """
    idx = tab['Nvisits'] > 150.
    print('Nfields',len(tab[idx]))
    lons=tab['RA']
    lats=tab['Dec']


    fig = plt.figure()
        #fig.suptitle('Season '+str(season+1),fontsize=12)

    ax = fig.add_subplot(211)
    m = Basemap(projection='moll',lon_0=180)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    m.drawmeridians(meridians,labels=[0,1,1,0],fontsize=10)
    for i in np.arange(len(meridians)):
        plt.annotate(np.str(int(meridians[i]))+'$^o$',xy=m(meridians[i],30),xycoords='data')

    x, y = m(lons,lats)
    
    norm=tab['Nvisits']
    m.scatter(x,y,c=norm,marker='s',s=10,cmap=plt.cm.jet)
 
    toprint='WFD fields - Season '+str(season)
    plt.annotate(toprint, xy=(0.30, 1.1), xycoords='axes fraction')
    plt.colorbar(fraction=0.02, pad=0.04)

    ax = fig.add_subplot(212)

        #print norm
    fontsize=12
    ax.hist(norm,bins=int(np.max(norm)),histtype='step',range=(0,int(np.max(norm))))
        #ax.set_xlim(xmin=0.)
    ax.set_ylabel('Number of Entries',{'fontsize': fontsize})
    ax.set_xlabel('Number of visits',{'fontsize': fontsize})
    #plt.gcf().savefig('Obs_Plots/Cadence_WFD_all_'+str(season)+'.png')
    
    plt.show()


"""
name='Observations_'+fieldname+'_'+str(fieldid)+'.txt'

for seas in range(1,2):
    myobs=Observations(fieldid=fieldid, filename=thedir+suffix+'/Season_'+str(seas)+'/'+fieldname+'/'+name)
    print len(myobs.seasons)
    for season in range(len(myobs.seasons)):
    #for season in range(2):
        obs_tot=myobs.seasons[season]
        idx = obs_tot['band']=='LSSTPG::i'
        obs=obs_tot[idx]
        diff=obs['mjd'][1:]-obs['mjd'][:-1]
        print np.min(diff),np.max(diff),np.max(obs['mjd'])-np.min(obs['mjd']),np.median(diff)
        plt.plot(obs_tot['mjd'],obs_tot['airmass'],'bo')
"""

def Ana_File(fi,j,output_q):
    r=[]
    #time_begin=time.time()
    summary=np.load(fi)
    idxa = summary['type']=='pixel'
    idxb = summary['type']=='point'
    for val in summary[idxa]:
        r.append((val['RA'],val['Dec'],len(summary[idxb])))
    #print('after proc',time.time()-time_begin,len(summary))
    if output_q is not None:
        output_q.put({j : r})

def Ana_File_new(fi,j=-1,output_q=None):
    r=[]
    time_begin=time.time()
    summary=np.load(fi)
    #print summary.dtype
    r.append((np.asscalar(summary['fieldRA']),np.asscalar(summary['fieldDec']),np.asscalar(summary['Nvisits'])))
    print('after proc',time.time()-time_begin,np.asscalar(summary['Nvisits']))
    if output_q is not None:
        output_q.put({j : r})

def Ana_File_Observations(fi,j=-1,output_q=None):
    r=[]
    time_begin=time.time()
    summary=np.load(fi)
    #print summary.dtype                                                                                                                                     
    #r.append((np.asscalar(summary['fieldRA']),np.asscalar(summary['fieldDec']),np.asscalar(summary['Nvisits'])))
    #print('after proc',time.time()-time_begin)
    if output_q is not None:
        output_q.put({j : summary})


def Multiproc(n_per_batch,sumfiles,target):
    n_per_batch=n_per_batch

    inter=range(0,len(sumfiles),n_per_batch)
    inter=np.append(inter,len(sumfiles))

    restot=[]
    time_begin=time.time()
#print range(len(sumfiles)-1)
    for jo in range(len(inter)-1):
        ida=inter[jo]
        idb=inter[jo+1]
        result_queue = multiprocessing.Queue()

        if (jo%10) == 0:
            print('Processing',jo,time.time()-time_begin)

        for j in range(ida,idb):
        #print val
            p=multiprocessing.Process(name='Subprocess-'+str(j),target=target,args=(sumfiles[j],j,result_queue))
        #restot=Get_coadd(val,key,season)
            p.start()
 
        resultdict = {}
    
        for j in range(ida,idb):
            resultdict.update(result_queue.get())

        for p in multiprocessing.active_children():
            p.join()

        for j in range(ida,idb):
            restot+=resultdict[j]

    return restot


def Multiproc_Observations(n_per_batch,sumfiles,target):
    n_per_batch=n_per_batch

    inter=range(0,len(sumfiles),n_per_batch)
    inter=np.append(inter,len(sumfiles))

    restot={}
    time_begin=time.time()
#print range(len(sumfiles)-1)                                                                                                                               
    for jo in range(len(inter)-1):
        ida=inter[jo]
        idb=inter[jo+1]
        result_queue = multiprocessing.Queue()

        if (jo%10) == 0:
            print('Processing obs',jo,time.time()-time_begin)

        for j in range(ida,idb):
        #print val                                                                                                                                          
            p=multiprocessing.Process(name='Subprocess-'+str(j),target=target,args=(sumfiles[j],j,result_queue))
        #restot=Get_coadd(val,key,season)                                                                                                                   
            p.start()

        resultdict = {}
        rescon=None
        for j in range(ida,idb):
            resultdict.update(result_queue.get())

        for p in multiprocessing.active_children():
            p.join()

        for j in range(ida,idb):
            if rescon is None:
                rescon=resultdict[j]
            else:
                rescon=np.concatenate((rescon,resultdict[j]))

        restot[jo]=rescon
        #if jo > 2000:
            #break
    
    return restot
            

parser = OptionParser()
parser.add_option("--simu_name", type="string", default='feature_baseline_10yrs', help="filter [%default]")
parser.add_option("--nside", type="int", default=64, help="filter [%default]")
parser.add_option("--year", type="int", default=2, help="filter [%default]")

opts, args = parser.parse_args()

thedir='OpSimFiles'
fieldname='DD'
fieldid=''
suffix='_'+opts.simu_name

nside=opts.nside
year=opts.year

#dir_files='/sps/lsst/data/dev/pgris/'+thedir+suffix+'/N_healpix_'+str(nside)
dir_files='/sps/lsst/data/dev/pgris/'+thedir+suffix+'/Season_'+str(year)
if nside <  0:

    sumfiles=glob.glob(dir_files+'_10.0'+'/Summary_*.npy')
    
    restot=Multiproc(10,sumfiles,target=Ana_File)
    res=np.rec.fromrecords(restot,names=['RA','Dec','Nvisits'])
    Plot_Cadence(res,'Nvisits',2)

else:


    what='Observations'
    sumfilesb=glob.glob(dir_files+'/N_healpix_'+str(nside)+'/'+what+'*.npy')

    #print('hello',sumfilesb,dir_files+'/N_healpix_'+str(nside)+'/'+what+'*.npy')
    
    if what == 'Summary':
        restotb=Multiproc(10,sumfilesb,target=Ana_File_new)
        resb=np.rec.fromrecords(restotb,names=['RA','Dec','Nvisits','MJD'])
        Plot_Cadence(resb,'Nvisits',2)
    else:
        resb=Multiproc_Observations(10,sumfilesb,target=Ana_File_Observations)
        for key,val in resb.items():
            plt.plot(val['observationStartMJD'],val['airmass'],'k.')
        

    """
    Ana_File_new('/sps/lsst/data/dev/pgris/OpSimFiles_feature_baseline_10yrs/Season_2/N_healpix_64/Observations_1200.npy')
    """
    

"""
print obs['fieldRA']
plt.plot(obs['fieldRA'],obs['fieldDec'],'ko')
"""


plt.show()
