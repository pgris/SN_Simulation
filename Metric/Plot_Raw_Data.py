import numpy as np
import glob
from Observations import * 
import pylab as pl
import multiprocessing

def Multiproc_Load(nbatch,files):

    ivals=range(0,len(files),nbatch)
    ivals=np.append(ivals,len(files))
    
    restot=None
    for i in range(len(ivals)-1):
    
        ia=ivals[i]
        ib=ivals[i+1]
        
        result_queue = multiprocessing.Queue()
        #print('processing',ia,ib)
        if ia > 0 and (ia%5000)==0:
            print('Dump in file',i)
            

        for j in range(ia,ib):
            #print('alors',files[j])
            p=multiprocessing.Process(name='Subprocess-'+str(j),target=Load_file,args=(files[j],j,result_queue))
            p.start()

        resultdict = {}
    
        for j in range(ia,ib):
            resultdict.update(result_queue.get())

        for p in multiprocessing.active_children():
            p.join()
    
        for j in range(ia,ib):
            if resultdict[j] is not None:
                if restot is None:
                    if resultdict[j] is not None:
                        restot=resultdict[j]
                else:
                    if resultdict[j] is not None:
                        restot=np.concatenate((restot,resultdict[j]))

    return restot

def Load_file(fi,j,out_q=None):

    #print('here',fi)
    fieldid=(fi.split('/')[-1]).split('_')[1]
    fieldid=int(fieldid.split('.')[0])
    obs=Observations(fieldid,filename=fi,season_length=365.,names=['observationStartMJD','fieldRA','fieldDec'])

    #print(len(obs.seasons))
    res=None
    if len(obs.seasons) > 0:
        res=obs.seasons[0]
    
    if out_q is not None:
        out_q.put({j : res})
    else:
        return res   


def Plot(x,y,legx='',legy=''):

    """
    fig = pl.figure()
    grid = pl.GridSpec(2, 2)
    print(grid)
    ax = fig.add_subplot(grid[0, :])
    pl.plot(tab[whatx],tab[whaty],'k.')

    ax = fig.add_subplot(grid[1, 0])
    ax.hist(tab[whatx],histtype='step')
    
    ax = fig.add_subplot(grid[1, 1])
    ax.hist(tab[whaty],histtype='step')
    """
    fig = pl.figure(figsize=(10, 10))
    grid = pl.GridSpec(4, 4, hspace=0.8, wspace=0.8)
    main_ax = fig.add_subplot(grid[:-1, 1:])
    y_hist = fig.add_subplot(grid[:-1, 0], xticklabels=[], sharey=main_ax)
    x_hist = fig.add_subplot(grid[-1, 1:], yticklabels=[], sharex=main_ax)

    # scatter points on the main axes
    #main_ax.plot(x, y, 'ok', markersize=3, alpha=0.2)
    main_ax.plot(x, y, 'ok', markersize=1.)

    # histogram on the attached axes
    x_hist.hist(x, 40, histtype='step',
                orientation='vertical', color='red')
    #x_hist.invert_yaxis()

    y_hist.hist(y, 40, histtype='step',
                orientation='horizontal', color='blue')
    y_hist.invert_xaxis()

    main_ax.set_xlabel(legx)
    main_ax.set_ylabel(legy)





season=3
fieldname='WFD'

thedir='../../Read_sqlite_python3/OpSimLogs_alt_sched_rolling_all/'
thedir+=fieldname+'/'
thedir+='Season_'+str(season)+'/'

files=glob.glob(thedir+fieldname+'_*.txt')

#print('hello',files[10])
tot_seas=Multiproc_Load(10,files)
"""
for fi in files:

    fieldid=(fi.split('/')[-1]).split('_')[1]
    fieldid=int(fieldid.split('.')[0])
        
    obs=Observations(fieldid,filename=fi,season_length=365.,names=['observationStartMJD','fieldRA','fieldDec'])

    print(len(obs.seasons))
    if len(obs.seasons) > 0:
        seas=obs.seasons[0]
        
        if tot_seas is None:
            tot_seas=seas
        else:
            tot_seas=np.concatenate((tot_seas,seas))
"""

#pl.plot(tot_seas['fiveSigmaDepth'],tot_seas['finSeeing'],'k.')
for band in 'gri':
    idx = tot_seas['filter']== band
    sel=tot_seas[idx]
    Plot(sel['fiveSigmaDepth'],sel['finSeeing'],'m$_{5}$ [mag]','seeing [\'\']')
    Plot(sel['fiveSigmaDepth'],sel['skyBrightness'],'m$_{5}$ [mag]','msky [mag]')
    Plot(sel['fiveSigmaDepth'],sel['cloud'],'m$_{5}$ [mag]','cloud')

pl.show()
    
