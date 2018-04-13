import numpy as np
import pickle as pkl
from mpl_toolkits.basemap import Basemap
import pylab as plt
import heapq
import multiprocessing
import time
import os
from optparse import OptionParser
import glob

parallels = np.arange(-90.,90,30.)
meridians = np.arange(0.,360.,60.)

def median(sel,var):
    
    if var != '':
        selb=np.sort(sel,order=var)
        selb=selb[var]
    else:
        selb=np.sort(sel)

    num=len(selb)

    #print('Median on',num,selb)
    if num >=5:
        n_lower=num/2-1.96*np.sqrt(float(num))/2.
        n_upper=1+num/2+1.96*np.sqrt(float(num))/2.

        if int(n_lower) < 0:
            n_lower = 0
        if int(n_upper) >=len(selb):
            n_upper=len(selb)-1
        #print 'hello',num,np.median(selb),int(n_lower),int(n_upper),selb[int(n_lower)],selb[int(n_upper)],selb
        res=np.median(selb),selb[int(n_lower)],selb[int(n_upper)]
    else:
        return -1.,-1.,-1.

    return res

def Plot_Ra_Dec():
    fieldnames=['WFD','NorthEclipticSpur-18c','SouthCelestialPole-18','GalacticPlane','DD']
    colors=['b','y','g','k','r']

    corresp=dict(zip(fieldnames,['WFD','NorthEclipticSpur','SouthCelestialPole','GalacticPlane','DD']))
    m = Basemap(projection='moll',lon_0=180)
    parallels = np.arange(-90.,90,30.)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
    meridians = np.arange(0.,360.,60.)
    m.drawmeridians(meridians,labels=[0,1,1,0],fontsize=10)
    for i in np.arange(len(meridians)):
        plt.annotate(np.str(int(meridians[i]))+'$^o$',xy=m(meridians[i],30),xycoords='data')

    for i,fieldname in enumerate(fieldnames):
        fich=fieldname+'_Ra_Dec.txt'
        res=np.loadtxt(fich,dtype={'names': ('fieldid','Ra','Dec'),'formats': ('i4','f8','f8')})
    #print res['Ra']
        x, y = m(np.rad2deg(res['Ra']),np.rad2deg(res['Dec']))
        m.scatter(x,y,marker='s',color=colors[i], label=corresp[fieldname])
    

    plt.legend(bbox_to_anchor=(0.05, -0.20), ncol=3,loc=2, borderaxespad=0.,fontsize=12.)   
    #plt.gcf().savefig('Map_Plots/Ra_Dec.png')

def Plot_Hist(what,tab_all,season):

    conv=dict(zip([band for band in 'ugrizy'],[1,2,3,4,5,6]))
    class_map = lambda x: conv[x]
    figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
    for band in 'ugrizy':
        idx = (tab_all['band']== band)&(tab_all['season']== season)&(tab_all['median_cadence']>0.)
        sel=tab_all[idx]
        print len(sel),season,band,sel['median_cadence'],sel['fieldid'],np.median(sel['median_cadence'])
        axa.hist(sel['median_cadence'],histtype='step',bins=35)
        #axa.scatter(sel['median_cadence'],map(class_map,sel['band']))

    axa.set_xlim([0,35])

def Plot_Map(what,tab_all,band,season,thresh,legend):

#print tab['ra']


    idx = (tab_all['band']==band)&(tab_all['season']==season)&(tab_all[what]>thresh)
    tab=tab_all[idx]
    
    #print tab['fieldid'],len(tab['fieldid'])

    lons = np.rad2deg(tab['ra'])
    lats = np.rad2deg(tab['dec'])

    m = Basemap(projection='moll',lon_0=180)

    x, y = m(lons,lats)
#x, y = m(*np.meshgrid(lons, lats))
#m.drawmapboundary(fill_color='#99ffff')
    parallels = np.arange(-90.,90,30.)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
    meridians = np.arange(0.,360.,60.)
    m.drawmeridians(meridians,labels=[0,1,1,0],fontsize=10)
    for i in np.arange(len(meridians)):
        plt.annotate(np.str(int(meridians[i]))+'$^o$',xy=m(meridians[i],10),xycoords='data')
#m.fillcontinents(color='#cc9966',lake_color='#99ffff')
    m.scatter(x,y,c=tab[what],marker='s',s=10,cmap=plt.cm.jet)
#print 'hello',len(x),tab['median_cadence']
#cs = m.contourf(x,y,tab['median_cadence'],cmap=plt.cm.jet)
#m.pcolormesh(x, y,tab['median_cadence'] , cmap='Paired')
    toprint=legend+' - '+band+' band - season '+str(season+1)
    plt.annotate(toprint, xy=(0.30, 1.1), xycoords='axes fraction')
    plt.colorbar(fraction=0.02, pad=0.04)
    plt.grid()

def Old_Plot():
    tab_all=None

    for fieldname in ['WFD']:
        if tab_all is None:
            tab_all=pkl.load(open(fieldname+'.pkl','rb'))
        else:
            tab_all=np.concatenate((tab_all,pkl.load(open(fieldname+'.pkl','rb'))))


    what='median_cadence'
    for season in range(10):
    #Plot_Map('m5',tab_all,'r',season,10.,'5-$\sigma$ depth')
        """
        for band in 'ugrizy':
        Plot_Map('median_cadence',tab_all,band,season,0.,'Cadence')
        #Plot_Hist('median_cadence',tab_all,season)
        plt.show()
        """
        for band in 'ugrizy':
            idx = (tab_all['band']==band)&(tab_all['season']==season)&(tab_all[what]>-0.5)
            tab=tab_all[idx]
            print season,band,np.percentile(np.array(tab['median_cadence']),50.),np.median(np.array(tab['median_cadence']))
         
    idb = (tab_all['fieldid'] == 311)&(tab_all['season'] == 0)
    tab=tab_all[idb]
    print len(tab),tab[['fieldid','band','median_cadence']]

#Plot_Ra_Dec()
def Plot_Cadence_old(tab,what):

    radec=np.unique((tab[['Ra','Dec']]))
    print radec[:]
    
    norm=[]
    normb=[]
    for val in radec:
        idd = (tab['Ra']==val[0])&(tab['Dec']==val[1])
        sel=tab[idd]
        #if what == 'cadence':
        normb.append(np.median(sel['mjd'][1:]-sel['mjd'][:-1]))
        #if what == 'visits':
        norm.append(np.median(sel['Nvisits']))

    lons = np.rad2deg([radec[i][0] for i in range(len(radec))])
    lats = np.rad2deg([radec[i][1] for i in range(len(radec))])
    
    m = Basemap(projection='moll',lon_0=180)
    
    x, y = m(lons,lats)
    
    m.scatter(x,y,c=norm,marker='s',s=10,cmap=plt.cm.jet)
    
    pp=[x for x in norm if str(x) != 'nan']
    v = np.linspace(np.min(pp), np.max(pp), 10, endpoint=True)
    #cbar = plt.colorbar(fraction=0.02, pad=0.04)
    print pp
    cbar = plt.colorbar(ticks=v,fraction=0.02, pad=0.04)
    if what == 'visits':
        cbar.ax.set_ylabel('Median number of Visits per day')
    if what == 'cadence':
        cbar.ax.set_ylabel('Median cadence ')

    mb = Basemap(projection='moll',lon_0=180)
    
    xb, yb = mb(lons,lats)
    
    mb.scatter(xb,yb,c=normb,marker='s',s=10,cmap=plt.cm.jet)
    

    plt.grid()


def Process_Vals(sel,season,j,out_q):

    r=[]
    sel.sort(order='mjd')
        #print 'hhh',sel
    diff=sel['mjd'][1:]-sel['mjd'][:-1]
    med_cad,low_cad,upp_cad=median(diff,'')
    medmjd={}
    for bband in 'ugrizy':
        #print(sel['mjd_'+bband])
        idx = sel['mjd_'+bband] > 0.
        pol=sel[idx]
        if len(pol) > 1:
            pol.sort(order='mjd_'+bband)
            medmjd[bband]=np.median(pol['mjd_'+bband][1:]-pol['mjd_'+bband][:-1])
            #print('hhh',bband,pol['mjd_'+bband][1:]-pol['mjd_'+bband][:-1])
        else:                           
            medmjd[bband]=-1.0

    #print('hello',med_cad,np.median(sel['mjd'][1:]-sel['mjd'][:-1]),sel['band'],sel['Nvisits'],sel['mjd'])
    duration=np.max(sel['mjd'])-np.min(sel['mjd'])
    m5=np.median(sel['m5sigmadepth'])
    medm5={}
    #print('hello',len(sel))
    visits={}
    visits_daily={}
    for bband in 'ugrizy':
        #if bband not in visits_tot.keys():
            #visits_tot[bband]=[]
        visits[bband]=np.sum(sel['Nvisits_'+bband])
        #visits_tot[bband].append(np.sum(sel['Nvisits_'+bband]))
        idd=sel['Nvisits_'+bband] > 0.
        selday=sel[idd]
        if len(selday) > 0:
            if len(selday) ==1:
                visits_daily[bband]=selday['Nvisits_'+bband]
            else:
                visits_daily[bband]=np.median(selday['Nvisits_'+bband])
        else:
            visits_daily[bband]=0
        #print('hello',bband,sel['Nvisits_'+bband])
        idc= sel['median_m5sigmadepth_'+bband] >=0.
        ssel=sel[idc]
        if len(ssel) > 1:
            medm5[bband]=np.median(ssel['median_m5sigmadepth_'+bband])
        else:
           medm5[bband]=-1. 
        

    airmass=np.median(sel['airmass'])
    nvisits=np.sum(sel['Nvisits'])
    """
    if len(sel)>= 2:
    largest=heapq.nlargest(2,diff)
    idxa = diff == largest[0]
    idxb = diff == largest[1]
    largest_1=largest[0]
    largest_2=largest[1]
    diff_time_1=np.asarray(sel['mjd'][idxa]-np.min(sel['mjd']))
    diff_time_2=np.asarray(sel['mjd'][idxb]-np.min(sel['mjd']))
    else:
    largest_1=-1
    largest_2=-1
    diff_time_1=-1
    diff_time_2=-1
    """
    largest_1=-1
    largest_2=-1
    diff_time_1=-1
    diff_time_2=-1
    fi=[season,sel['band'][0],np.mean(sel['Ra']),np.mean(sel['Dec']),np.median(sel['Nvisits']),med_cad,low_cad,upp_cad,int(sel['fieldid'][0]),duration,m5,airmass,largest_1,largest_2,diff_time_1,diff_time_2,nvisits]
    fi+=[medm5[b] for b in 'ugrizy']
    fi+=[medmjd[b] for b in 'ugrizy']
    fi+=[visits[b] for b in 'ugrizy'] 
    fi+=[visits_daily[b] for b in 'ugrizy']
    #fi+=[np.sum(visits_tot[b]) for b in 'ugrizy']
    #print('hhh',visits,medmjd)
    #r.append((season,sel['band'][0],np.mean(sel['Ra']),np.mean(sel['Dec']),np.median(sel['Nvisits']),med_cad,low_cad,upp_cad,int(sel['fieldid'][0]),duration,m5,airmass,largest_1,largest_2,diff_time_1,diff_time_2,nvisits,[medm5[b] for b in 'ugrizy']))
    r.append(tuple(fi))


    names=['season','band','Ra','Dec','median_visits','median_cadence','median_cadence_low','median_cadence_upp','fieldid','duration','median_m5','median_airmass','max_timeobs_1','max_timeobs_2','diff_time_1','diff_time_2','Nvisits','median_m5_u','median_m5_g','median_m5_r','median_m5_i','median_m5_z','median_m5_y','median_cadence_u','median_cadence_g','median_cadence_r','median_cadence_i','median_cadence_z','median_cadence_y','Nvisits_u','Nvisits_g','Nvisits_r','Nvisits_i','Nvisits_z','Nvisits_y','Nvisits_day_u','Nvisits_day_g','Nvisits_day_r','Nvisits_day_i','Nvisits_day_z','Nvisits_day_y']
    #print('names',len(names))
     
    dtypes=[]
    
    for val in names:
        if val != 'band':
            dtypes.append((val,'f8'))
        else:
            dtypes.append((val,'S15'))
    #return np.rec.fromrecords(r,names=['season','band','Ra','Dec','median_visits','median_cadence','median_cadence_low','median_cadence_upp','fieldid','duration','median_m5','median_airmass'])
    #print 'seaeon ',season,r
    res=np.rec.array(r,dtype=dtypes)
    if out_q is not None:
        out_q.put({j : res})
    else:
        return res

def Get_Vals(tab,season,pkl_file,band=''):

    fieldids=np.unique((tab['fieldid']))

    print('number of fields',len(fieldids),fieldids)
    #fieldids=fieldids[0:500]
    nbatch=15
    #nbatch=1
    #fieldids=[586]

    ivals=range(0,len(fieldids),nbatch)
    ivals=np.append(ivals,len(fieldids))

    restot=None
    time_begin=time.time()
    for i in range(len(ivals)-1):

        ia=ivals[i]
        ib=ivals[i+1]

        result_queue = multiprocessing.Queue()

        if (i%100)==0 and i > 1:
            print('Processing and dumping',i,time.time()-time_begin,len(restot))
            pkl.dump(restot, pkl_file)
            restot=None

        for j in range(ia,ib):
            idd = (tab['fieldid']==fieldids[j])
            sel=tab[idd]
            p=multiprocessing.Process(name='Subprocess-'+str(j),target=Process_Vals,args=(sel,season,j,result_queue))
            p.start()

        resultdict = {}
    
        for j in range(ia,ib):
            resultdict.update(result_queue.get())

        for p in multiprocessing.active_children():
            p.join()
    
        for j in range(ia,ib):
            if restot is None:
                restot=resultdict[j]
            else:
                restot=np.concatenate((restot,resultdict[j]))

        
    if restot is not None:
        print('Final dumping',len(restot))
        pkl.dump(restot, pkl_file)
    #return restot

def Get_Vals_old(tab,season,band=''):
    radec=np.unique((tab[['Ra','Dec']]))
    r=[]
    print('Nfields',len(radec))
    for val in radec:
       
        idd = (tab['Ra']==val[0])&(tab['Dec']==val[1])
        if band !='':
            idd&=(tab['band']=='LSSTPG::'+band)
        sel=tab[idd]
        #print 'here sel',len(sel)
        sel.sort(order='mjd')
        #print 'hhh',sel
        diff=sel['mjd'][1:]-sel['mjd'][:-1]
        med_cad,low_cad,upp_cad=median(diff,'')
        #print 'hello',med_cad,np.median(sel['mjd'][1:]-sel['mjd'][:-1]),sel['band'],sel['Nvisits']
        duration=np.max(sel['mjd'])-np.min(sel['mjd'])
        m5=np.median(sel['m5sigmadepth'])
        airmass=np.median(sel['airmass'])
        """
        if len(sel)>= 2:
            largest=heapq.nlargest(2,diff)
            idxa = diff == largest[0]
            idxb = diff == largest[1]
            largest_1=largest[0]
            largest_2=largest[1]
            diff_time_1=np.asarray(sel['mjd'][idxa]-np.min(sel['mjd']))
            diff_time_2=np.asarray(sel['mjd'][idxb]-np.min(sel['mjd']))
        else:
            largest_1=-1
            largest_2=-1
            diff_time_1=-1
            diff_time_2=-1
        """
        largest_1=-1
        largest_2=-1
        diff_time_1=-1
        diff_time_2=-1
        r.append((season,sel['band'][0],sel['Ra'][0],sel['Dec'][0],np.median(sel['Nvisits']),med_cad,low_cad,upp_cad,int(sel['fieldid'][0]),duration,m5,airmass,largest_1,largest_2,diff_time_1,diff_time_2))

    names=['season','band','Ra','Dec','median_visits','median_cadence','median_cadence_low','median_cadence_upp','fieldid','duration','median_m5','median_airmass','max_timeobs_1','max_timeobs_2','diff_time_1','diff_time_2']

    dtypes=[]
    
    for val in names:
        if val != 'band':
            dtypes.append((val,'f8'))
        else:
            dtypes.append((val,'S15'))
    #return np.rec.fromrecords(r,names=['season','band','Ra','Dec','median_visits','median_cadence','median_cadence_low','median_cadence_upp','fieldid','duration','median_m5','median_airmass'])
    #print 'seaeon ',season,r
    return np.rec.array(r,dtype=dtypes)

def Plot_Cadence(tab,season,band):

    radec=np.unique((tab[['Ra','Dec']]))
    #print radec[:]
    
    norm=[]
    normb=[]
   
    for val in radec:
        idd = (tab['Ra']==val[0])&(tab['Dec']==val[1])
        sel=tab[idd]
        #if what == 'cadence':
        normb.append(np.median(sel['mjd'][1:]-sel['mjd'][:-1]))
        #if what == 'visits':
        norm.append(np.median(sel['Nvisits']))
       
    lons = np.rad2deg([radec[i][0] for i in range(len(radec))])
    lats = np.rad2deg([radec[i][1] for i in range(len(radec))])
    
    lons=[radec[i][0] for i in range(len(radec))]
    lats=[radec[i][1] for i in range(len(radec))]

    fig = plt.figure()
    fig.suptitle('Season '+str(season+1)+' - '+band+' band', fontsize=12)

    ax = fig.add_subplot(211)
    m = Basemap(projection='moll',lon_0=180)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    m.drawmeridians(meridians,labels=[0,1,1,0],fontsize=10)
    for i in np.arange(len(meridians)):
        plt.annotate(np.str(int(meridians[i]))+'$^o$',xy=m(meridians[i],30),xycoords='data')

    x, y = m(lons,lats)
    
    m.scatter(x,y,c=norm,marker='s',s=10,cmap=plt.cm.jet)
    
    pp=[x for x in norm if str(x) != 'nan']
    v = np.linspace(np.min(pp), np.max(pp), 10, endpoint=True)
    #cbar = plt.colorbar(fraction=0.02, pad=0.04)
    v=[float(int(x*10.))/10. for x in v]
    cbar = plt.colorbar(ticks=v,fraction=0.02, pad=0.04)
    #cbar.ax.set_ylabel('Median number of Visits per day')
    ax.set_title('Median number of Visits per day', fontsize=10)
    ax = fig.add_subplot(212)

    mb = Basemap(projection='moll',lon_0=180)
    mb.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    mb.drawmeridians(meridians,labels=[0,1,1,0],fontsize=10)
    for i in np.arange(len(meridians)):
        plt.annotate(np.str(int(meridians[i]))+'$^o$',xy=m(meridians[i],30),xycoords='data')
    xb, yb = mb(lons,lats)
    
    mb.scatter(xb,yb,c=normb,marker='s',s=10,cmap=plt.cm.jet)
    ax.set_title('Median cadence [day-1]', fontsize=10)
    pp=[x for x in normb if str(x) != 'nan']
    v = np.linspace(np.min(pp), np.max(pp), 10, endpoint=True)
    #cbar = plt.colorbar(fraction=0.02, pad=0.04)
    #print pp
    cbar = plt.colorbar(ticks=v,fraction=0.02, pad=0.04)

    #plt.grid()    

def nonan(val,cutval):
    return [x for x in val if str(x) != 'nan' if x >= cutval]

def Plot_Cadence_DD(tab_sel,fieldid):
    figa, axa = plt.subplots(ncols=1, nrows=3, figsize=(10,9))
    figa.suptitle('DD - fieldid '+str(fieldid))
    tot_label=[]
    fontsize=12.
    
    myls=['-','--']
    colors=dict(zip([i for i in range(10)],['k','k','r','r','b','b','g','g','m','m']))
    corresp_inverted=dict(zip([i for i in range(len(bands))],bands))
    print corresp_inverted
    ylim=[]
    for season in range(0,10):
    #for band in 'ugrizy':
        idx = (tab_sel['season']== season)&(tab_sel['median_cadence']>0.)
        sel=tab_sel[idx]
        print len(sel),season,sel['median_cadence'],sel['fieldid'],np.median(sel['median_cadence']),sel['band'],map(class_map,sel['band'])
        axa[0].plot(map(class_map,sel['band']),sel['median_cadence'],ls=myls[season%2],color=colors[season])

        ll='Y'+str(season+1)
        tot_label.append(axa[1].errorbar(map(class_map,sel['band']),sel['median_cadence_low'],ls=myls[season%2],color=colors[season],label=ll))
        
        axa[2].errorbar(map(class_map,sel['band']),sel['median_cadence_upp'],ls=myls[season%2],color=colors[season])

    axa[0].set_ylabel('Median cadence \n [day$^{-1}$]',{'fontsize': fontsize})
    axa[1].set_ylabel('95% C.L. Lower cadence \n [day$^{-1}$]',{'fontsize': fontsize})  
    axa[2].set_ylabel('95% C.L. Upper cadence \n [day$^{-1}$]',{'fontsize': fontsize})

    labs = [l.get_label() for l in tot_label]
    axa[1].legend(tot_label, labs, ncol=5,loc='lower right',prop={'size':fontsize},frameon=False)
    for io in range(3):
        #axa[io].set_xlim([-0.1,5.1])
        axa[io].set_xticklabels([corresp_inverted[i] for i in range(len(bands))])
        #axa[io].set_xlim([-0.1,6.1])
        ylim = axa[io].get_ylim()
        axa[io].set_ylim([ylim[0]-0.5,ylim[1]+0.5])
        axa[io].set_xlabel('band',{'fontsize': fontsize})

    plt.gcf().savefig('Obs_Plots/Cadence_'+str(fieldid)+'.png')

def Plot_Vars_DD(tab_sel,fieldid):
    figa, axa = plt.subplots(ncols=1, nrows=3, figsize=(10,9))
    figa.suptitle('DD - fieldid '+str(fieldid))
    tot_label=[]
    fontsize=12.
    
    myls=['-','--']
    colors=dict(zip([i for i in range(10)],['k','k','r','r','b','b','g','g','m','m']))
    corresp_inverted=dict(zip([i for i in range(len(bands))],bands))
    print corresp_inverted
    ylim=[]
    for season in range(10):
    #for band in 'ugrizy':
        idx = tab_sel['season']== season
        sel=tab_sel[idx]
        
        axa[0].plot(map(class_map,sel['band']),sel['duration'],ls=myls[season%2],color=colors[season])

        ll='Y'+str(season+1)
        tot_label.append(axa[1].errorbar(map(class_map,sel['band']),sel['median_m5'],ls=myls[season%2],color=colors[season],label=ll))
        
        axa[2].errorbar(map(class_map,sel['band']),sel['median_airmass'],ls=myls[season%2],color=colors[season])
        #axa[2].errorbar(map(class_map,sel['band']),sel['max_diff_timeobs'],ls=myls[season%2],color=colors[season])
    axa[0].set_ylabel('Season length \n [day$^{-1}$]',{'fontsize': fontsize})
    axa[1].set_ylabel('Median 5$\sigma$-depth [mag]',{'fontsize': fontsize})  
    axa[2].set_ylabel('Median airmass',{'fontsize': fontsize})
    #axa[2].set_ylabel('Max $\Delta T$ (two obs)\n [day]',{'fontsize': fontsize})
    labs = [l.get_label() for l in tot_label]
    axa[1].legend(tot_label, labs, ncol=5,loc='lower right',prop={'size':fontsize},frameon=False)
    for io in range(3):
        #axa[io].set_xlim([-0.1,5.1])
        axa[io].set_xticklabels([corresp_inverted[i] for i in range(len(bands))])
        #axa[io].set_xlim([-0.1,6.1])
        ylim = axa[io].get_ylim()
        axa[io].set_ylim([ylim[0]-0.5,ylim[1]+0.5])
        axa[io].set_xlabel('band',{'fontsize': fontsize})

    plt.gcf().savefig('Obs_Plots/Vars_'+str(fieldid)+'.png')

def Plot_Varsb_DD(tab_sel,fieldid):
    figa, axa = plt.subplots(ncols=2, nrows=2, figsize=(10,9))
    figa.suptitle('DD - fieldid '+str(fieldid))
    tot_label=[]
    fontsize=12.
    
    myls=['-','--']
    colors=dict(zip([i for i in range(10)],['k','k','r','r','b','b','g','g','m','m']))
    corresp_inverted=dict(zip([i for i in range(len(bands))],bands))
    print corresp_inverted
    ylim=[]
    for season in range(10):
    #for band in 'ugrizy':
        idx = tab_sel['season']== season
        sel=tab_sel[idx]
        
        

        ll='Y'+str(season+1)
        tot_label.append(axa[0][0].errorbar(map(class_map,sel['band']),sel['max_timeobs_1'],ls=myls[season%2],color=colors[season],label=ll))
        axa[0][1].errorbar(map(class_map,sel['band']),sel['diff_time_1'],ls=myls[season%2],color=colors[season],label=ll)
        
        axa[1][0].errorbar(map(class_map,sel['band']),sel['max_timeobs_2'],ls=myls[season%2],color=colors[season])
        axa[1][1].errorbar(map(class_map,sel['band']),sel['diff_time_2'],ls=myls[season%2],color=colors[season])
        #axa[2].errorbar(map(class_map,sel['band']),sel['max_diff_timeobs'],ls=myls[season%2],color=colors[season])
    axa[0][0].set_ylabel('1st max $\Delta T^{obs}$  [day]',{'fontsize': fontsize})
    axa[0][1].set_ylabel('$\Delta T (1st max,min^{obs})$ [day]',{'fontsize': fontsize})  
    axa[1][0].set_ylabel('2nd max $\Delta T^{obs}$  [day]',{'fontsize': fontsize})
    axa[1][1].set_ylabel('$\Delta T (2nd max,min^{obs})$ [day]',{'fontsize': fontsize})
    #axa[2].set_ylabel('Max $\Delta T$ (two obs)\n [day]',{'fontsize': fontsize})
    labs = [l.get_label() for l in tot_label]
    axa[0][0].legend(tot_label, labs, ncol=3,loc='upper center',prop={'size':fontsize},frameon=False)
    for io in range(2):
        for jo in range(2):
        #axa[io].set_xlim([-0.1,5.1])
            axa[io][jo].set_xticklabels([corresp_inverted[i] for i in range(len(bands))])
        #axa[io].set_xlim([-0.1,6.1])
            ylim = axa[io][jo].get_ylim()
            axa[io][jo].set_ylim([ylim[0]-0.5,ylim[1]+0.5])
            axa[io][jo].set_xlabel('band',{'fontsize': fontsize})

    plt.gcf().savefig('Obs_Plots/Varsb_'+str(fieldid)+'.png')

def Plot_Cadence_all_WFD(tab,season,what,leg,legunit,extent):
    
    radec=np.unique((tab[['Ra','Dec']]))
    fieldids=np.unique((tab['fieldid']))
    norm=[]
    #print('radec',radec)
    ra_filt=[]
    dec_filt=[]
    for fieldid in fieldids:
        #idd = (tab['Ra']==val[0])&(tab['Dec']==val[1])
        idd=tab['fieldid']==fieldid
        sel=tab[idd]
        #norm.append(np.median(sel['median_cadence']))
        #print('hallo',sel['median_cadence'],val[0],val[1],len(sel))
        #print('here',what,sel[what],sel[['Ra','Dec']],sel['fieldid'])
        ro=np.asscalar(sel[what])
        if ro > 0:
            ra_filt.append(sel['Ra'])
            dec_filt.append(sel['Dec'])
            norm.append(ro)
    #print(norm)
    lons = np.rad2deg([radec[i][0] for i in range(len(radec))])
    lats = np.rad2deg([radec[i][1] for i in range(len(radec))])
    lons=np.rad2deg(ra_filt)
    lats=np.rad2deg(dec_filt)

    #lons = [radec[i][0] for i in range(len(radec))]
    #lats = [radec[i][1] for i in range(len(radec))]

    #lons=ra_filt
    #lats=dec_filt


    fig = plt.figure(figsize=(10,9))
        #fig.suptitle('Season '+str(season+1),fontsize=12)

    ax = fig.add_subplot(211)
    m = Basemap(projection='moll',lon_0=180)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    m.drawmeridians(meridians,labels=[0,1,1,0],fontsize=10)
    for i in np.arange(len(meridians)):
        plt.annotate(np.str(int(meridians[i]))+'$^o$',xy=m(meridians[i],30),xycoords='data')

    x, y = m(lons,lats)
    
    m.scatter(x,y,c=norm,marker='s',s=20,edgecolor='none',cmap=plt.cm.jet)
 
    toprint='WFD fields - Season '+str(season+2)+' \n '+ extent+'\n '+leg
    plt.annotate(toprint, xy=(0.30, 1.1), xycoords='axes fraction')
    plt.colorbar(fraction=0.02, pad=0.04)

    ax = fig.add_subplot(212)

        #print norm
    fontsize=12
    #print(norm)
    if what == 'median_cadence':
        ax.hist(norm,bins=int(np.max(norm)),histtype='step',range=(0,int(np.max(norm))))
    else:
        ax.hist(norm,bins=20,histtype='step')
    #ax.hist(norm,bins=5)
    #print(norm)
        #ax.set_xlim(xmin=0.)
    ax.set_ylabel('Number of Entries',{'fontsize': fontsize})
    ax.set_xlabel(leg+' '+legunit,{'fontsize': fontsize})
    name_plot='Cadence_WFD_'+what+'_'+str(season)+'.png'
    outdir='Obs_Plots/'+extent
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    plt.gcf().savefig(outdir+'/'+name_plot)
    
    #plt.show()

def Plot_Cadence_Hist_WFD(tab,what):

    print('rr',len(tab),tab['band'])

    for w in what:
        for band in 'ugrizy':
            idx = tab['band']=='LSSTPG::'+band
            Plot_Hist(tab[idx],w+'_'+band)

def Plot_m5_Summary(tab):

    fieldids=np.unique(tab['fieldid'])

    print('test',tab.dtype)
    r=[]
    bands='ugrizy'
    for fieldid in fieldids:
        idx = tab['fieldid']==fieldid
        sel=tab[idx]
        for band in bands:
            what='median_m5_'+band
            med,do,up=median(sel,what)
            r.append((np.mean(sel['Ra']),np.mean(sel['Dec']),band,med,med-do,up-med))
            if do < 10:
                print('pb here ?',med,do,up)
    rt=np.rec.fromrecords(r,names=['Ra','Dec','band','med','med_low','med_upp'])
    print(rt)

    for band in bands:
        idx = (rt['band']==band)&(rt['med_low']<20.)
        sel=rt[idx]
        Plot_Moll(sel['Ra'],sel['Dec'],sel['med'],'',0,band+' band','Median 5$\sigma$ depth [mag]','median_m5_'+band)
        #Plot_Moll(sel['Ra'],sel['Dec'],sel['med_low'],'',0,band+' band','Median 5$\sigma$ depth [mag]')
        #Plot_Moll(sel['Ra'],sel['Dec'],sel['med_upp'],'',0,band+' band','Median 5$\sigma$ depth [mag]')

def Plot_Summary(tab,what,legx,savename):

    fieldids=np.unique(tab['fieldid'])

    print('test',tab.dtype)
    r=[]
    bands='ugrizy'
    #fieldids=[586]
    for fieldid in fieldids:
        idx = tab['fieldid']==fieldid
        sel=tab[idx]
        for band in bands:
            whatb=what+'_'+band
            med,do,up=median(sel,whatb)
            r.append((np.mean(sel['Ra']),np.mean(sel['Dec']),band,med,med-do,up-med))
            print(what,med,do,up)
    rt=np.rec.fromrecords(r,names=['Ra','Dec','band','med','med_low','med_upp'])
    print(rt)

    for band in bands:
        idx = (rt['band']==band)&(rt['med_low']<20.)
        sel=rt[idx]
        Plot_Moll(sel['Ra'],sel['Dec'],sel['med'],'',-1,band+' band',legx,savename+'_'+band)
        #Plot_Moll(sel['Ra'],sel['Dec'],sel['med_low'],'',0,band+' band','Lower 95% C.L. '+legx)
        #Plot_Moll(sel['Ra'],sel['Dec'],sel['med_upp'],'',0,band+' band','Upper 95% C.L. '+legx)

    

        

def Plot_Cadence_Summary(tab):

    r=[]
    
    for season in range(10):
        idb = tab['season']==season
        sel=tab[idb]
        for band in 'ugrizy':
            med,do,up=median(sel,'median_cadence_'+band)
            r.append((band,med,up,do,season))
    rt=np.rec.fromrecords(r,names=['band','med','med_upp','med_low','season'])

    fontsize=12
    for band in 'ugrizy':
        figa, axa = plt.subplots()
        ik = rt['band']==band
        sel=rt[ik]

        axa.plot(sel['season']+1,sel['med'],color='k',label='median')
        axa.plot(sel['season']+1,sel['med_upp'],color='k',ls='--',label='95% C.L. upper')
        axa.plot(sel['season']+1,sel['med_low'],color='k',ls='-.',label='95% C.L. lower')

        axa.set_xlabel('Season',{'fontsize': fontsize})
        axa.set_ylabel('Median cadence [day-1]',{'fontsize': fontsize})
        axa.legend(loc='best',prop={'size':fontsize},frameon=False)
        axa.set_title(band+' band')
        plt.gcf().savefig('Obs_Plots_Summary/Median_Cadence_'+band+'.png')


def Plot_Hist(tabo,what):

    figa, axa = plt.subplots()
    
    max_norm=np.max(tabo[what])
    min_norm=np.min(tabo[what])
    print('go',what,len(what),max_norm)
    for season in range(10):
        idb = tabo['season']==season
        sel=tabo[idb]
        if what.count('cadence'):
            axa.hist(sel[what],bins=int(max_norm),histtype='step',range=(0,int(max_norm)))
        else:
            axa.hist(sel[what],bins=10,histtype='step',normed=True,stacked=True,range=(min_norm,max_norm))



def Plot_Cadence_all_WFD_multiple(tab,season,what,leg,legunit,extent):
    
    radec=np.unique((tab[['Ra','Dec']]))
    fieldids=np.unique((tab['fieldid']))
    norm={}
    ra_filt={}
    dec_filt={}

    for w in what:
        norm[w]=[]
        ra_filt[w]=[]
        dec_filt[w]=[]
    #print('radec',radec)
    
    for fieldid in fieldids:
        #idd = (tab['Ra']==val[0])&(tab['Dec']==val[1])
        idd=tab['fieldid']==fieldid
        sel=tab[idd]
        #norm.append(np.median(sel['median_cadence']))
        #print('hallo',sel['median_cadence'],val[0],val[1],len(sel))
        #print('here',what,sel[what],sel[['Ra','Dec']],sel['fieldid'])
        for w in what:
            ro=np.asscalar(sel[w])
            if ro > 0:
                ra_filt[w].append(sel['Ra'])
                dec_filt[w].append(sel['Dec'])
                norm[w].append(ro)
   
    for i,w in enumerate(what):
        Plot_Moll(ra_filt[w],dec_filt[w],norm[w],extent,season,leg[i],legunit[i])


    #plt.show()

def Plot_Moll(ra,dec,norm,extent,season,leg,legunit,savename):

    #fig = plt.figure(figsize=(10,9))
    fig=plt.figure()
        #fig.suptitle('Season '+str(season+1),fontsize=12)
    lons=np.rad2deg(ra)
    lats=np.rad2deg(dec)
    #ax = fig.add_subplot(211)
    m = Basemap(projection='moll',lon_0=180)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    m.drawmeridians(meridians,labels=[0,1,1,0],fontsize=10)
    for i in np.arange(len(meridians)):
        plt.annotate(np.str(int(meridians[i]))+'$^o$',xy=m(meridians[i],30),xycoords='data')

    x, y = m(lons,lats)
    
    m.scatter(x,y,c=norm,marker='s',s=20,edgecolor='none',cmap=plt.cm.jet)
 
    if season >=0:
        toprint='WFD fields - Season '+str(season+1)+' \n '+ extent+'\n '+leg
    else:
       toprint='WFD fields - All seasons \n '+ extent+'\n '+leg 
    plt.annotate(toprint, xy=(0.30, 1.1), xycoords='axes fraction')
    cbar=m.colorbar(location='bottom',pad="5%")
    cbar.set_label(legunit)
    plt.gcf().savefig('Obs_Plots_Summary/'+savename+'.png')

def Plot_Bi(tab,season,whatx,legx,legunitx,whaty,legy,legunity,extent):

    figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(7,5))
    toprint='WFD fields - Season '+str(season+2)+' \n '+ extent
    figa.suptitle(toprint)
    idx=(tab[whatx]>0)&(tab[whaty]>0)
    sel=tab[idx]
    axa.plot(sel[whatx],sel[whaty],'k.')
    
    axa.set_xlabel(legx+' '+legunitx)
    axa.set_ylabel(legy+' '+legunity)
    name_plot='Cadence_WFD_'+whatx+'_'+whaty+'_'+str(season)+'.png'
    outdir='Obs_Plots/'+extent
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    plt.gcf().savefig(outdir+'/'+name_plot)


parser = OptionParser()
parser.add_option("--simu_name", type="string", default='feature_baseline_10yrs', help="filter [%default]")
parser.add_option("--year", type="int", default=2, help="filter [%default]")
parser.add_option("--fieldname", type="string", default='WFD', help="filter [%default]")

opts, args = parser.parse_args()


fieldname=opts.fieldname

#tab_all=pkl.load(open('Cadence_'+fieldname+'_all.pkl','rb'))

add='_all'
#add=''
extent=opts.simu_name
"""
extent='feature_baseline_10yrs'
extent='feature_rolling_half_mask_10yrs'
extent='feature_rolling_twoThird_10yrs'
extent='alt_sched'
extent='alt_sched_rolling'
"""
#num=dict(zip(['feature_baseline_10yrs','feature_rolling_half_mask_10yrs','feature_rolling_twoThird_10yrs','alt_sched','alt_sched_rolling',''],[8,5,5,6,4,2]))
tab_all=None

finame=''
if extent != '':
    finame='Cadence_'+extent+'_'+fieldname+add
else:
    finame='Cadence_'+fieldname+add

files=glob.glob(finame+'*.pkl')

print('there pal',files)
for fi in files: 
    f=open(fi,'rb')
    objs = []
    while 1:
        try:
        #print 'trying'
            objs.append(pkl.load(f))
        except EOFError:
            break
        
    for obj in objs:
        if tab_all is None:
            tab_all=obj
        else:
            tab_all=np.concatenate((tab_all,obj))

print(np.unique(tab_all['season']))

bands='ugrizy'
#bands='g'

bbands=['rgizy']

if add == '':
    tot_tab=None
    conv=dict(zip(['LSSTPG::'+b for b in bands],range(len(bands))))
    class_map = lambda x: conv[x]
    for season in range(1):
        print 'season',season
    
        for band in bands:
            idx = (tab_all['season']==season)&(tab_all['band']=='LSSTPG::'+band)
            tab=tab_all[idx]
            
            Plot_Cadence(tab,season,band)
            res=Get_Vals(tab,season,band)
            print band,res
            fig = plt.figure()
            fig.suptitle(band+' band - Season '+str(season))
            ax = fig.add_subplot(211)
            ax.hist(nonan(res['median_visits'],0.),bins=50,histtype='step')
            ax.set_xlabel('N$_{visits}$ (median)')
            ax.set_ylabel('Number of Entries')
            ax = fig.add_subplot(212)
            ax.hist(nonan(res['median_cadence'],0.),bins=50,histtype='step') 
            ax.set_xlabel('Median cadence [days-1]')
            ax.set_ylabel('Number of Entries')
            plt.show()
            
            #tab=tab_all[idx]
            
            """
            res=Get_Vals(tab,season,band)
            if tot_tab is None:
                tot_tab=np.array(res,dtype=res.dtype)
            else:
                tot_tab=np.vstack([tot_tab,np.array(res,dtype=res.dtype)])
            """
    """
    #for fieldid in [290,744,1427,2412,2786]:
    #for fieldid in [2786]:
    for fieldid in [100,101,102,103]:
        ida=tot_tab['fieldid']==fieldid
        tab_sel=tot_tab[ida]
        print('selection',len(tab_sel))
        Plot_Cadence_DD(tab_sel,fieldid)
        Plot_Vars_DD(tab_sel,fieldid)
        Plot_Varsb_DD(tab_sel,fieldid)
    """
"""
tot_tab=None
for season in range(0,10):
    print 'season',season
    for band in bands:
        idx = (tab_all['season']==season)&(tab_all['band']=='LSSTPG::'+band)
        
        tab=tab_all[idx]
        
        res=Get_Vals(tab,season,band)
        if tot_tab is None:
            tot_tab=np.array(res,dtype=res.dtype)
        else:
           tot_tab=np.vstack([tot_tab,np.array(res,dtype=res.dtype)])

print tot_tab
conv=dict(zip(bands,range(len(bands))))
class_map = lambda x: conv[x]



for fieldid in [290,744,1427,2412,2786]:
    ida=tot_tab['fieldid']==fieldid
    tab_sel=tot_tab[ida]
    
    Plot_Cadence_DD(tab_sel,fieldid)    
    #Plot_Vars_DD(tab_sel,fieldid)
"""
tot_tab=None
"""
idbb= tab_all['fieldid']==321
tab_all=tab_all[idbb]
"""

if add.count('all'):
    tot_tab=None
    Process=False
    for season in range(10):
        print 'season',season
        if Process:
        #tot_tab=None
            pkl_file = open('Summary_'+extent+'_'+str(season)+'.pkl','wb')

            idx = (tab_all['season']==season)
            sel=tab_all[idx]
        #print np.unique(sel['band']),'LSST::'+bbands[0],type(sel),sel.dtype
            tot_season=None
            """
            for band in np.unique(sel['band']):
                print 'hello',band
            #if (fieldname=='DD' and band.split(':')[2].count('rg') and band.split(':')[2].count('zy')) or fieldname=='WFD':
                passed=True
                if fieldname=='DD' :
                    for b in 'rgizy':
                        passed&=band.split(':')[2].count(b)
            
           
                if passed:
                    idxb = sel['band']==band
                    selb=sel[idxb]
                
                    if tot_season is None:
                        tot_season=np.array(selb,dtype=sel.dtype)
                    else:
                        tot_season=np.concatenate((tot_season,np.array(selb,dtype=sel.dtype)))
             
       
                res=Get_Vals(tot_season,season)
              """
            print('number of data',len(sel))
            Get_Vals(sel,season,pkl_file)
           
            pkl_file.close()
    #print 'res',res['median_cadence']
            """
            if tot_tab is None:
            tot_tab=res
            else:
            tot_tab=np.concatenate((tot_tab,res))
            """
        #print(tot_tab,tot_tab.dtype)
        res=None
        filename='Summary_'+extent+'_'+str(season)+'.pkl'
        print('opening',filename)
        f=open(filename,'rb') 
        objs = []
        while 1:
            try:
        #print 'trying'
                objs.append(pkl.load(f))
            except EOFError:
                break
        for obj in objs:
            if res is None:
                res=obj
            else:
                res=np.concatenate((res,obj))
        
        print('Nfields',len(res))
        #Plot_Cadence_all_WFD(res,season,'median_cadence','Median Cadence','[Day$^{-1}$]',extent)
        #Plot_Cadence_all_WFD(res,season,'median_m5','Median m5','[mag]',extent)
        #Plot_Cadence_all_WFD(res,season,'median_airmass','Median airmass','',extent)
        #Plot_Cadence_all_WFD(res,season,'Nvisits','Number of visits','',extent) 
        #for bband in 'ugrizy':
        """
        for bband in 'g':
            what=['median_m5_'+bband,'median_cadence_'+bband,'Nvisits_'+bband]
            #leg=['Median m5-'+bband,'median_cadence_'+bband,'Median cadence -'+bband,'Nvisits_'+bband,'Number of Visits -'+bband]
            leg=[bband+' band']*len(what)
            legunit=['Median m5 [mag]','Median cadence [Day$^{-1}$]','Nvisits']
            Plot_Cadence_all_WFD_multiple(res,season,what,leg,legunit,extent)
            #Plot_Cadence_all_WFD(res,season,'median_m5_'+bband,'Median m5-'+bband,'[mag]',extent)
            #Plot_Cadence_all_WFD(res,season,'median_cadence_'+bband,'Median cadence -'+bband,'[Day$^{-1}$]',extent)
            #Plot_Cadence_all_WFD(res,season,'Nvisits_'+bband,'Number of Visits -'+bband,'',extent)
            #Plot_Bi(res,season,'median_m5_'+bband,'Median m5-'+bband,'[mag]','median_cadence_'+bband,'Median cadence -'+bband,'[Day$^{-1}$]',extent)
        #plt.show()
        """
    files=glob.glob('Summary_'+extent+'_*.pkl')
    
    #files=['Summary__2.pkl','Summary__2.pkl']
    print('files',files)
    rest=None
    objs = []
    for fi in files:
        f=open(fi,'rb') 
        
        while 1:
            try:
                objs.append(pkl.load(f))
            except EOFError:
                break
        
    for obj in objs:
        if rest is None:
            rest=obj
        else:
            rest=np.concatenate((rest,obj))

    #Plot_Cadence_Hist_WFD(rest,['median_cadence','median_m5'],)
    Plot_Cadence_Summary(rest)
    Plot_m5_Summary(rest)
    Plot_Summary(rest,'Nvisits','Median number of visits','Nvisits')
    #Plot_Summary(rest,'Nvisits_tot','Total number of visits')
    Plot_Summary(rest,'Nvisits_day','Median number of visits per obs day','Nvisits_day')
    plt.show()

    
#for fieldid in [290,744,1427,2412,2786]:
    if fieldname=='DD':
        figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        figa.suptitle('10*g+20*r+20*i+26*z+20*y')
        tot_label=[]
        fontsize=12
        for fieldid in [290,744,1427,2412,2786]:
        #for fieldid in [1427]:
#for fieldid in [1427]:
            ida=tot_tab['fieldid']==fieldid
            sel=tot_tab[ida]
            ll='Field '+str(fieldid)
            tot_label.append(axa.errorbar(sel['season']+1,sel['median_cadence'],ls='-',label=ll))

        labs = [l.get_label() for l in tot_label]
        axa.legend(tot_label, labs, ncol=3,loc='lower right',prop={'size':fontsize},frameon=False)
        axa.set_xlabel('Year',{'fontsize': fontsize})
        axa.set_ylabel('Median Cadence [Day$^{-1}$]',{'fontsize': fontsize})


#plt.show()
