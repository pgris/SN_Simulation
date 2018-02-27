from Generate_LC import Generate_LC
from Telescope import *
import multiprocessing
from Parameters import parameters
import pylab as pl
import glob
from Observations import *
import multiprocessing
import os 

def f5_cadence_lims_sncosmo(SNR=dict(zip([b for b in "griz"], 
                                 [30., 40., 30., 20.])),
                    zs=[0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1],
                    X1=-2., Color=0.2,
                            expTime=dict(zip([b for b in "griz"], 
                                               [30., 30., 30., 30.])),telescope=None):

    bands=SNR.keys()
    restframe_phase_range = (-20., 40.)	
    pmin, pmax = restframe_phase_range
    r=[]
    lims={}
    for z in zs:
        r.append((z,-2.0,0.2,0.))

    #params=np.rec.fromrecords(r, names = ['z', 'X1', 'Color', 'DayMax'])
    params=np.rec.array(r, dtype=[('z', 'f8'),('X1', 'f8'), ('Color', 'f8'),('DayMax','f8')])
    #print 'before',params['z']
    cadence=1.
    
    #m5_lim={'u':23.61,'g':24.83,'r':24.35,'i':23.88,'z':23.30,'y':22.43}

    #this m5_lim values correspond to an exposure time of 30s...needs to be corrected for DDF

    obs_param=parameters()
    m5_lim=obs_param.m5
    for key, vals in expTime.items():
        m5_lim[key]=m5_lim[key]+1.25*np.log10(vals/30.)

    for param in params:
        rb=[]
        z=param['z']
        lims[z]={}
        mjd_min = np.floor(pmin * (1.+z) + param['DayMax'])
        mjd_max = np.ceil(pmax * (1.+z) + param['DayMax'])
        mjd = np.arange(mjd_min, mjd_max, cadence)
        lc_sncosmo=Generate_LC(param,telescope=telescope)
        #lc_sncosmo=Generate_LC(param,mjd,expTime).lc
        #airmass_val=np.repeat(airmass, len(mjd), 0)
       
        for mjd_val in mjd:
            for band in bands:
                rb.append((mjd_val,band,obs_param.seeing[band],m5_lim[band],expTime[band],obs_param.msky[band],))
        table_obs=np.rec.fromrecords(rb,names=['mjd','band','seeing','m5sigmadepth','exptime','sky'])
        lc=lc_sncosmo(table_obs)
        #print(table_obs['seeing'])
        """
        result_queue = multiprocessing.Queue()
        process=[]
        for b in bands:
            expTime_val=np.repeat(expTime[b], len(mjd), 0)
            m5_lim_val=np.repeat(m5_lim[b], len(mjd), 0)
        #ret[b]=lc_sncosmo(mjd,b,expTime[b],result_queue)
            p=multiprocessing.Process(name='Subprocess-'+b,target=lc_sncosmo,args=(mjd,airmass_val,m5_lim_val,b,expTime_val,result_queue))
            process.append(p)
            p.start()

        resultdict = {}
        for b in bands:
            resultdict.update(result_queue.get())

        for p in multiprocessing.active_children():
            p.join()
        """
        print('redshift',z)
        for filtre in bands:
           
            #res=resultdict[filtre][0]
            idx=lc['band']==filtre
            res=lc[idx]
            #print('bou',filtre,res['flux_e']/res['flux_e_err'])
            idx = res['flux_e']/res['flux_e_err'] > 5.
            res=res[idx]
            Li2=np.sqrt(np.sum(np.square(res['flux_e'])))
            lim = 5. * Li2 / SNR[filtre]
            #print(filtre,lim,Li2,SNR[filtre])
            lims[z][filtre]=lim

    return lims 

def f5_cadence_plot(lsstpg, band, lims=None, median_log_summary=None, 
                    median_log_summaryb=None,dict_marker=None,
                    lims_sncosmo=None,
                    mag_range=(23., 26.5), 
                    dt_range=(0.5, 15.), 
                    target={},targetb={},alpha=1.,simu_name=''):
    #                    SNR=dict(zip(['LSSTPG::' + b for b in 'grizy'], 
    #                                 [20., 20., 20., 20., 10.]))):
    
    dt = np.linspace(dt_range[0], dt_range[1], 50)
    m5 = np.linspace(mag_range[0], mag_range[1], 50)
    b = [band] * len(m5)
    f5 = lsstpg.mag_to_flux(m5, b)
    
    F5,DT = np.meshgrid(f5, dt)
    M5,DT = np.meshgrid(m5, dt)
    metric = np.sqrt(DT) * F5
    #print('metric',metric)
    sorted_keys=np.sort([k  for  k in  lims.keys()])[::-1]

    # draw limits 
    pl.figure()
    pl.imshow(metric, extent=(mag_range[0], mag_range[1], 
                              dt_range[0], dt_range[1]), 
              aspect='auto', alpha=0.25)
    
    if median_log_summary is not None:
        idx = median_log_summary['band'] == band
        m = median_log_summary[idx]
        pl.plot(m['m5'], m['cadence'], 'r+',alpha=alpha)
        m5_exp=np.median(m['m5'])
        #print 'h1',band,np.median(m['m5'])
        
    if median_log_summaryb is not None:
        """
        idx = median_log_summaryb['band'] == band
        m = median_log_summaryb[idx]
        pl.plot(m['m5'], m['cadence'], 'b+',alpha=alpha)
        m5_exp=np.median(m['m5'])
        """
        #print 'h2',band,np.median(m['m5'])
        for key, val in median_log_summaryb.items():
            idx = val['band'] == band
            m = val[idx]
            pl.plot(m['m5'], m['cadence'],dict_marker[key],alpha=alpha)
            m5_exp=np.median(m['m5'])

    """
    if median_log_summaryc is not None:
        idx = median_log_summaryc['band'] == band
        m = median_log_summaryc[idx]
        pl.plot(m['m5'], m['cadence'], 'rs',alpha=0.7)
        #print 'h3',band,np.median(m['m5']),m['m5']

    if median_log_summaryd is not None:
        idx = median_log_summaryd['band'] == band
        m = median_log_summaryd[idx]
        pl.plot(m['m5'], m['cadence'], 'bs',alpha=0.7)
        #print 'h4',band,np.median(m['m5'])
     """
    cadence_ref=3
    
    if lims is not None:
        fmt = {}
        ll = [lims[zz][band] for zz in sorted_keys]
        print('limits',ll)
        cs = pl.contour(M5, DT, metric, ll, colors='k')
        dict_target_snsim=Get_target(cs,sorted_keys,cadence_ref,m5_exp)

        strs = ['$z=%3.1f$' % zz for zz in sorted_keys]
        for l,s in zip(cs.levels, strs):
            fmt[l] = s
        pl.clabel(cs, inline=True, fmt=fmt, fontsize=16, use_clabeltext=True)
    
    
    if lims_sncosmo is not None:

        llc = [lims_sncosmo[zz][band] for zz in sorted_keys]
        b = [band]*len(llc)
        print 'baba',band,llc,[zz for zz in sorted_keys],lsstpg.flux_to_mag(llc,b)
        csc = pl.contour(M5, DT, metric, llc, colors='b') 
        dict_target_sncosmo=Get_target(csc,sorted_keys,cadence_ref,m5_exp)

        strs = ['$z=%3.1f$' % zz for zz in sorted_keys]
        for l,s in zip(csc.levels, strs):
            fmt[l] = s
        pl.clabel(csc, inline=1, fmt=fmt, fontsize=16)

    print 'Band',band
    #ref_z=dict(zip([b for b in 'rizy'],[0.8,1.0,0.7,0.6]))
    #ref_z=dict(zip([b[-1:] for b in target.keys()],[0.3,0.5,0.4,0.4]))
    """
    ref_z=dict(zip([b[-1:] for b in targetb.keys()],[val for val in targetb.values()]))
    print 'allo',targetb.keys(),ref_z
    key=ref_z[band[-1:]]
    snsim_target=dict_target_snsim[key]
    if lims_sncosmo is not None:
        sncosmo_target=dict_target_sncosmo[key]
        print key,snsim_target[1]-snsim_target[2],sncosmo_target[1]-sncosmo_target[2]
    else:
        print key,snsim_target[1]-snsim_target[2]

    """
    t = target.get(band, None)
    print target, t

    if t is not None:
        pl.plot(t[0], t[1], 
                color='r', marker='*', 
                markersize=15)
    
    """
    pl.plot(snsim_target[1], cadence_ref, 
                color='r', marker='*', 
                markersize=15)
    """
    pl.xlabel('$m_{5\sigma}$', fontsize=18)
    pl.ylabel(r'Observer frame cadence $^{-1}$ [days]', fontsize=18)
    pl.title('$%s$ - %s' % (band.split(':')[-1],simu_name), fontsize=18)
    pl.xlim(mag_range)
    pl.ylim(dt_range)
    pl.grid(1)
    pl.gcf().savefig('Plots_'+simu_name+'/cadence_'+band+'.png', bbox_inches='tight')

    if median_log_summary is not None:
        
        figa, axa = pl.subplots(ncols=2, nrows=1,figsize=(15,9))
        figa.suptitle('$%s$ - %s' % (band.split(':')[-1],simu_name), fontsize=18)
        idx = (median_log_summary['band'] == band)&(median_log_summary['cadence']<=30)
        m = median_log_summary[idx]
        #print(len(m),m['cadence'])
        axa[0].hist(m['m5'],bins=20,histtype='step')
        axa[1].hist(m['cadence'],bins=range(int(np.min(m['cadence'])),int(np.max(m['cadence'])),1),histtype='step')
        axa[0].set_xlabel('$m_{5\sigma}$', fontsize=18)
        axa[0].set_ylabel('Number of Events', fontsize=18)
        axa[1].set_xlabel('Observer frame cadence $^{-1}$ [days]', fontsize=16)
        axa[1].set_ylabel('Number of Events', fontsize=18)

        #m5_exp=np.median(m['m5'])
        
        #print 'h1',band,np.median(m['m5'])

        pl.gcf().savefig('Plots_'+simu_name+'/histo_'+band+'.png', bbox_inches='tight')


def Get_target(cs,sorted_keys,cadence_ref,m5_exp):
    dict_target={}
    ptot = cs.allsegs
    for i,p in enumerate(ptot):
        m5_ref=-1
        diffref=999.
        for val in p:
            for vv in val:
                #print 'oh',vv
                diff=np.abs(cadence_ref-vv[1])
                if diff < diffref:
                    diffref=diff
                    m5_ref=vv[0]
        #print sorted_keys[i],'Target for snsim',cadence_ref,'days-1:',band,diffref,m5_ref,'exp',m5_exp,'diff',m5_ref-m5_exp
        dict_target[sorted_keys[i]]=(diffref,m5_ref,m5_exp)

    return dict_target

def Median_log_summary(fi,bands,j,out_q=None):
    
    r=[]
   
    fieldid=(fi.split('/')[-1]).split('_')[1]
    fieldid=int(fieldid.split('.')[0])
        
    obs=Observations(fieldid,filename=fi,season_length=365.,names=['observationStartMJD','fieldRA','fieldDec'])
        #print(fi,len(obs.seasons))
    for i in range(len(obs.seasons)):
        season=obs.seasons[i]
            #print(i,season)
        for band in bands:
            idx = season['filter']==band
            sel=season[idx]
            if len(sel) >= 2:
                m5=np.median(sel['fiveSigmaDepth'])
                airmass=np.median(sel['airmass'])
                sky=np.median(sel['skyBrightness'])
                seeing=np.median(sel['finSeeing'])
                moon=np.median(sel['moonPhase'])
                cadence=np.median(sel['observationStartMJD'][1:]-sel['observationStartMJD'][:-1])
                ra_sn=np.median(sel['RA_SN'])
                dec_sn=np.median(sel['Dec_SN'])
                r.append((fieldid,ra_sn,dec_sn,band,cadence,m5,airmass,sky,seeing,moon))
                
    res=None
    if len(r) > 0:
        res=np.rec.fromrecords(r,names=['fieldid','RA_SN','Dec_SN','band','cadence','m5','airmass','skyBrightness','finSeeing','moonPhase'])
    if out_q is not None:
        out_q.put({j : res})
    else:
        return res   
    
def Median_log_summary_old(files,bands,j,out_q=None):
    
    r=[]
    #print('files',files)
    #for fi in files[:1000]:
    for fi in files:
        print('hello',fi)
        fieldid=(fi.split('/')[-1]).split('_')[1]
        fieldid=int(fieldid.split('.')[0])
        
        obs=Observations(fieldid,filename=fi,season_length=365.,names=['observationStartMJD','fieldRA','fieldDec'])
        #print(fi,len(obs.seasons))
        for i in range(len(obs.seasons)):
            season=obs.seasons[i]
            #print(i,season)
            for band in bands:
                idx = season['filter']==band
                sel=season[idx]
                if len(sel) >= 2:
                    m5=np.median(sel['fiveSigmaDepth'])
                    airmass=np.median(sel['airmass'])
                    sky=np.median(sel['skyBrightness'])
                    seeing=np.median(sel['finSeeing'])
                    moon=np.median(sel['moonPhase'])
                    cadence=np.median(sel['observationStartMJD'][1:]-sel['observationStartMJD'][:-1])
                    ra_sn=np.median(sel['RA_SN'])
                    dec_sn=np.median(sel['Dec_SN'])
                    r.append((fieldid,ra_sn,dec_sn,band,cadence,m5,airmass,sky,seeing,moon))
                
    
    res=np.rec.fromrecords(r,names=['fieldid','RA_SN','Dec_SN','band','cadence','m5','airmass','skyBrightness','finSeeing','moonPhase'])
    if out_q is not None:
        out_q.put({j : res})
    else:
        return res      

def MultiProc_Median_log_summary(files,bands,simu_name,nbatch=20):

    restot=None
    ivals=range(0,len(files),nbatch)
    ivals=np.append(ivals,len(files))
    dirout='Summary_'+simu_name
    if not os.path.isdir(dirout) :
        os.makedirs(dirout)

    ifich=-1

    for i in range(len(ivals)-1):
    #for i in range(174,175,1):
        ia=ivals[i]
        ib=ivals[i+1]

        result_queue = multiprocessing.Queue()
        #print('processing',ia,ib)
        if ia > 0 and (ia%5000)==0:
            print('Dump in file',i)
            ifich+=1
            np.save(dirout+'/Sum_'+str(ifich)+'.npy',restot)
            restot=None

        for j in range(ia,ib):
            #print('alors',files[j])
            p=multiprocessing.Process(name='Subprocess-'+str(j),target=Median_log_summary,args=(files[j],bands,j,result_queue))
            p.start()

        resultdict = {}
    
        for j in range(ia,ib):
            resultdict.update(result_queue.get())

        for p in multiprocessing.active_children():
            p.join()
    
        for j in range(ia,ib):
            if resultdict[j] is not None:
                if restot is None:
                    restot=resultdict[j]
                else:
                    restot=np.concatenate((restot,resultdict[j]))


    if restot is not None:
        print('Dump in file final')
        ifich+=1
        np.save(dirout+'Sum_'+str(ifich)+'.npy',restot)
    
    #return restot
        
 
def Plot(med_log,bands,whata,whatb):

    for band in bands:
        idxb=med_log['band']==band
        sel=med_log[idxb]
        figa, axa = pl.subplots(ncols=1, nrows=1)
        figa.suptitle(band+' band')
        axa.plot(sel[whata],sel[whatb],'b.')
        axa.set_xlabel(whata)
        axa.set_ylabel(whatb)

def Plot_Stat(med_log,bands,whata,simu_name):

    r=[]
    corresp_band=dict(zip([b for b in 'ugrizy'],[i for i in range(0,7)]))
    print('corr',corresp_band)
    for band in bands:
        idxbb=med_log['band']==band
        sel=med_log[idxbb]
        for val in range(2,8,1):
            idxc=sel['cadence']<=val+1
            if val < 7:
                idxb=(sel['cadence']>=val)&(sel['cadence']<val+1)
            else:
                idxb=(sel['cadence']>=val)
            ssel=sel[idxb]
            sselc=sel[idxc]
            r.append((band,corresp_band[band],val,float(len(ssel))/float(len(sel)),float(len(sselc))/float(len(sel))))
    
    print(r)
    res=np.rec.fromrecords(r, names = ['band','iband','cadence','frac','frac_upper'])
    dict_plot={}
    dict_plot['legend_colz']='Fraction of Events'
    dict_plot['title']='Cadence vs band \n '+simu_name
    dict_plot['xt']=range(1,5)
    dict_plot['xticks']=['g','r','i','z']
    dict_plot['yt']=range(2,8)
    dict_plot['yticks']=['2 days','3 days','4 days','5 days',' 6 days',' >= 7 days' ]
    dict_plot['typesel']='int'
    dict_plot['ylims']=[1.5,7.5]
    dict_plot['simu_name']=simu_name
    dict_plot['save_name']='Cadence_Distrib'

    dict_plotn=dict_plot.copy()
    dict_plotn['yt']=range(2,8)
    dict_plotn['yticks']=['$\leq$ 3 days','$\leq$ 4 days','$\leq$ 5 days','$\leq$  6 days',' $\leq$ 7 days' , ' $\leq$ 8 days']
    """
    dict_plot['typesel']='int'
    dict_plot['ylims']=[1.5,7.5]
    dict_plot['simu_name']=simu_name
    """
    dict_plotn['save_name']='Cadence_Distrib_upper'
    
    #Plot_Map(res,'iband','cadence','frac',width=0.5,legend_colz='Fraction of Events',title='Cadence vs band',xt=range(1,5),xticks=['g','r','i','z'],yt=range(2,7),yticks=['2 days','3 days','4 days','5 days',' >= 6 days' ],typesel='int',ylims=[1.5,6.5])
    Plot_Map(res,'iband','cadence','frac',width_x=0.5,width_y=0.5,dict_plot=dict_plot)
    Plot_Map(res,'iband','cadence','frac_upper',width_x=0.5,width_y=0.5,dict_plot=dict_plotn)
    
    ra_step=10.
    dec_step=10.
    for band in bands:
        ro=[]
        idxb=med_log['band']==band
        sel=med_log[idxb]
       
        #figa, axa = pl.subplots(ncols=1, nrows=1)
        #axa.plot(sel['RA_SN'],sel['Dec_SN'],'ro')
        for ra in np.arange(np.min(sel['RA_SN']),np.max(sel['RA_SN'])+ra_step,ra_step):
            for dec in np.arange(np.min(sel['Dec_SN']),np.max(sel['Dec_SN'])+dec_step,dec_step):
                idxc=(sel['RA_SN']>=ra)&(sel['RA_SN']<ra+ra_step)
                idxc&=(sel['Dec_SN']>=dec)&(sel['Dec_SN']<dec+dec_step)
                ro.append((band,ra,dec,len(sel[idxc])))

        resb=np.rec.fromrecords(ro, names = ['band','ra','dec','N_SN'])
        #print resb
        dict_plotb={}
        dict_plotb['legend_colz']='Number of SN'
        dict_plotb['title']=band+' band - '+simu_name
        dict_plotb['typesel']='float'
        dict_plotb['xlims']=[np.min(sel['RA_SN']),np.max(sel['RA_SN'])]
        dict_plotb['ylims']=[np.min(sel['Dec_SN']),np.max(sel['Dec_SN'])]
        dict_plotb['xlab']='RA [deg]'
        dict_plotb['ylab']='Dec [deg]'
        dict_plotb['simu_name']=simu_name
        dict_plotb['save_name']='Control_'+band
        #print('band',band)
        Plot_Map_new(sel,'RA_SN','Dec_SN',width_x=ra_step,width_y=dec_step,dict_plot=dict_plotb)
        #break
    """
    x=[]
    y=[]

    width=0.5
    
    for val in res:
        x.append(val['iband']-width)
        y.append(val['cadence']+width)
            
    x.append(np.max(x)+2.*width)
    y.append(np.min(y)-2.*width)
    
    
    x=np.unique(x)
    y=np.unique(y)
 
    XX, YY = np.meshgrid(x, y)
        
    fig = plt.figure()
    #ax = fig.add_subplot(211)
    ax = fig.add_subplot()

    z = np.zeros((len(x)-1, len(y)-1))
        
    for i in range(len(x)-1):
        for j in range(len(y)-1):
            idx=Select(x[i]+width,y[j]+width,res,'iband','cadence')
            if len(res[idx]) > 0.:
                z[i,j]=res[idx]['frac']
            else:
                z[i,j]=0.
            print('here pal',i,j,x[i],y[j],z[i,j])
 
    plt.pcolormesh(XX,YY,z.T,cmap=plt.cm.jet)

    cbar=plt.colorbar(fraction=0.02, pad=0.04)
    cbar.ax.set_ylabel('Fraction of events', rotation=270,labelpad=20)
    plt.ylim([1.5,6.5])
    xt=range(1,5)
    plt.xticks(xt, ['g','r','i','z'])
    yt=range(2,7)
    plt.yticks(yt, ['2 days','3 days','4 days','5 days',' >= 6 days' ]) 
    plt.title('Cadence vs band')
    """
    """
        figa, axa = pl.subplots(ncols=1, nrows=1)
        figa.suptitle(band+' band')
        axa.plot(sel[whata],sel[whatb],'b.')
        axa.set_xlabel(whata)
        axa.set_ylabel(whatb)
    """
   
def Plot_Map(res,xstr,ystr,norm,width_x,width_y,dict_plot={}):

    x=[]
    y=[]

    
    for val in res:
        x.append(val[xstr]-width_x)
        y.append(val[ystr]+width_y)
            
    x.append(np.max(x)+2.*width_x)
    y.append(np.min(y)-2.*width_y)
    
    
    x=np.unique(x)
    y=np.unique(y)
 
    XX, YY = np.meshgrid(x, y)
        
    fig = plt.figure()
    #ax = fig.add_subplot(211)
    ax = fig.add_subplot()

    z = np.zeros((len(x)-1, len(y)-1))
        
    for i in range(len(x)-1):
        for j in range(len(y)-1):
            if dict_plot['typesel'] == 'int':
                idx=Select(x[i]+width_x,y[j]+width_y,res,xstr,ystr)
            else:
               idx=Selectf(x[i]+width_x,y[j]+width_y,res,xstr,ystr) 
            if len(res[idx]) ==1 :
                z[i,j]=res[idx][norm]
            else:
                if len(res[idx]) >0. :
                    print('hhhh',x[i],y[j],res[idx][norm])
                z[i,j]=0.
            #print('here pal',i,j,x[i],y[j],z[i,j])
 
    plt.pcolormesh(XX,YY,z.T,cmap=plt.cm.jet)
    v = np.linspace(np.min(z), np.max(z), 10, endpoint=True)
    cbar=plt.colorbar(fraction=0.02, pad=0.04,ticks=v)
    
    if dict_plot.has_key('legend_colz'):
        cbar.ax.set_ylabel(dict_plot['legend_colz'], rotation=270,labelpad=20)
    if dict_plot.has_key('xlims'):
        plt.xlim(dict_plot['xlims'])

    if dict_plot.has_key('ylims'):
        plt.ylim(dict_plot['ylims'])
    if dict_plot.has_key('xt'):
        plt.xticks(dict_plot['xt'],dict_plot['xticks'])
    if dict_plot.has_key('yt'):
        plt.yticks(dict_plot['yt'],dict_plot['yticks']) 
    if dict_plot.has_key('title'):
        plt.title(dict_plot['title'])
    if dict_plot.has_key('xlab'):
        plt.xlabel(dict_plot['xlab'])
    if dict_plot.has_key('ylab'):
        plt.ylabel(dict_plot['ylab'])
    if dict_plot.has_key('save_name'):    
        pl.gcf().savefig('Plots_'+dict_plot['simu_name']+'/'+dict_plot['save_name']+'.png', bbox_inches='tight')

def Plot_Map_new(data,xstr,ystr,width_x,width_y,dict_plot={}):

    x=[]
    y=[]

    xmin=round(np.min(data[xstr]))
    xmax=round(np.max(data[xstr]))
    ymin=round(np.min(data[ystr]))
    ymax=round(np.max(data[ystr]))
        
    x=np.arange(xmin,xmax,width_x)
    y=np.arange(ymin,ymax,width_y)
    
    #x.append(xmax+width_x)
    #y.append(ymax+width_y)
    
    #print('here',x,y)
    x=np.unique(x)
    y=np.unique(y)
 
    XX, YY = np.meshgrid(x, y)
        
    #print(XX,YY)
    fig = plt.figure()
    #ax = fig.add_subplot(211)
    ax = fig.add_subplot()

    z = np.zeros((len(x)-1, len(y)-1))
        
    for i in range(len(x)-1):
        for j in range(len(y)-1):
            idx=Select_new(x[i],x[i+1],y[j],y[j+1],data,xstr,ystr) 
            z[i,j]=len(data[idx])
           
            #z[i,j]=1.
            #print('here pal',i,j,x[i],x[i+1],y[j],y[j+1],z[i,j])
 
    plt.pcolormesh(XX,YY,z.T,cmap=plt.cm.jet)
    v = np.linspace(int(np.min(z)), int(np.max(z)), 10, endpoint=True,dtype='i8')
    cbar=plt.colorbar(fraction=0.02, pad=0.04,ticks=v)
    if 'legend_colz' in dict_plot.keys():
        cbar.ax.set_ylabel(dict_plot['legend_colz'], rotation=270,labelpad=20)
    if 'xlims' in dict_plot.keys():
        plt.xlim(dict_plot['xlims'])

    if 'ylims' in dict_plot.keys():
        plt.ylim(dict_plot['ylims'])
    if 'xt' in dict_plot.keys():
        plt.xticks(dict_plot['xt'],dict_plot['xticks'])
    if 'yt' in dict_plot.keys():
        plt.yticks(dict_plot['yt'],dict_plot['yticks']) 
    if 'title' in dict_plot.keys():
        plt.title(dict_plot['title'])
    if 'xlab' in dict_plot.keys():
        plt.xlabel(dict_plot['xlab'])
    if 'ylab' in dict_plot.keys():
        plt.ylabel(dict_plot['ylab'])
    if dict_plot.has_key('save_name'):    
        pl.gcf().savefig('Plots_'+dict_plot['simu_name']+'/'+dict_plot['save_name']+'.png', bbox_inches='tight')

def Plot_Map_old(res,xstr,ystr,norm,width,legend_colz,title='',xt=None,xticks=None,yt=None,yticks=None,typesel='int',xlims=None,ylims=None):

    x=[]
    y=[]

    
    for val in res:
        x.append(val[xstr]-width)
        y.append(val[ystr]+width)
            
    x.append(np.max(x)+2.*width)
    y.append(np.min(y)-2.*width)
    
    
    x=np.unique(x)
    y=np.unique(y)
 
    XX, YY = np.meshgrid(x, y)
        
    fig = plt.figure()
    #ax = fig.add_subplot(211)
    ax = fig.add_subplot()

    z = np.zeros((len(x)-1, len(y)-1))
        
    for i in range(len(x)-1):
        for j in range(len(y)-1):
            if typesel == 'int':
                idx=Select(x[i]+width,y[j]+width,res,xstr,ystr)
            else:
               idx=Selectf(x[i]+width,y[j]+width,res,xstr,ystr) 
            if len(res[idx]) > 0.:
                z[i,j]=res[idx][norm]
            else:
                z[i,j]=0.
            print('here pal',i,j,x[i],y[j],z[i,j])
 
    plt.pcolormesh(XX,YY,z.T,cmap=plt.cm.jet)

    cbar=plt.colorbar(fraction=0.02, pad=0.04)
    cbar.ax.set_ylabel(legend_colz, rotation=270,labelpad=20)
    if xlims is not None:
        plt.xlim(xlims)

    if ylims is not None:
        plt.ylim(ylims)
    if xt is not None:
        plt.xticks(xt,xticks)
    if yt is not None:
        plt.yticks(yt,yticks) 
    plt.title(title)
    

def Select(raval,decval,tab,xstr,ystr):

    min_ra=0.99999
    max_ra=1.00001
    min_dec=0.99999
    max_dec=1.00001
    if raval < 0:
        min_ra=1.00001
        max_ra=0.99999
    if decval < 0:
        min_dec=1.00001
        max_dec=0.99999

    idx=(tab[xstr] ==raval)
    idx&=(tab[ystr] >= min_dec*decval)&(tab[ystr]< max_dec*decval)

    return idx

def Selectf(raval,decval,tab,xstr,ystr):

    min_ra=0.999999
    max_ra=1.000001
    min_dec=0.999999
    max_dec=1.000001
    if raval < 0:
        min_ra=1.000001
        max_ra=0.999999
    if decval < 0:
        min_dec=1.000001
        max_dec=0.999999

    idx=(tab[xstr] >= min_ra*raval)&(tab[xstr]< max_ra*raval)
    idx&=(tab[ystr] >= min_dec*decval)&(tab[ystr]< max_dec*decval)

    return idx

def Select_new(ramin,ramax,decmin,decmax,tab,xstr,ystr):

    #print('hhhhe',ramin,ramax,decmin,decmax,xstr,ystr)
    idx=(tab[xstr] >= ramin)&(tab[xstr]< ramax)
    idx&=(tab[ystr] >= decmin)&(tab[ystr]< decmax)

    return idx

def Plotb(med_log,bands):
    r=[]

    for val in np.unique(med_log['fieldid']):
        ro=[]
        names=[]
        idx=med_log['fieldid']==val
        sel=med_log[idx]
        for band in bands:
            idxb=sel['band']==band
            selb=sel[idxb]
            for var in ['cadence','m5']:
                names.append(var+'_'+band)
                if len(selb) >= 1:
                    ro.append(np.asscalar(selb[var]))
                else:
                    ro.append(-1.)
        r.append(ro)

    print(r,names)
    res=np.rec.fromrecords(r, names = names)

    for ba in bands:
        for bb in bands:
            if bb != ba:
               figa, axa = pl.subplots(ncols=1, nrows=1)
               idx = (res['cadence_'+ba]<=30.)&(res['cadence_'+bb]<=30.)&(res['cadence_'+ba]>0.)&(res['cadence_'+bb]>0.)
               sel=res[idx]
               axa.plot(sel['cadence_'+ba],sel['cadence_'+bb],'k.')
               axa.set_xlabel('cadence_'+ba)
               axa.set_ylabel('cadence_'+bb)



atmos=False
instr=Telescope(atmos=atmos,aerosol=False,airmass=-1.)
instrb=Telescope(atmos=True,aerosol=False,airmass=1.2)

zs=[0.1, 0.2, 0.3, 0.4, 0.5]
SNR_sncosmo=dict(zip([b for b in "griz"],
                     [30., 40., 30., 20.]))

bands=SNR_sncosmo.keys()

target={'g': (24.83, 3.), # 23.86
        'r': (24.35, 3.), # 23.82
        'i': (23.88, 3.), # 23.51
        'z': (23.30, 3.)}

lims_sncosmo = f5_cadence_lims_sncosmo(zs=zs, SNR=SNR_sncosmo,telescope=instr)
lims_sncosmob = f5_cadence_lims_sncosmo(zs=zs, SNR=SNR_sncosmo,telescope=instrb)
lims_sncosmob = None
#                           bands=[instr_name + '::' + b for b in bands])

simu_name='alt_sched_rolling'
simu_name='alt_sched'
opsim_dir='/pbs/throng/lsst/users/gris/Read_sqlite_python3/OpSimLogs_'+simu_name+'/WFD'
opsim_files=glob.glob(opsim_dir+'/WFD_*.txt')

#median_log_summary=Median_log_summary(opsim_files,bands)

dump=False
if dump:
    MultiProc_Median_log_summary(opsim_files,bands,simu_name)

sum_files=glob.glob('Summary_'+simu_name+'/*.npy')

median_log_summary=None
for fi in sum_files:
    fich=np.load(fi)
    if median_log_summary is None:
        median_log_summary=fich
    else:
        median_log_summary=np.concatenate((median_log_summary,fich))
        

for bn in bands:
    f5_cadence_plot(instr, bn, # instr_name + '::' + bn,
                    lims=lims_sncosmo,lims_sncosmo=lims_sncosmob, 
                    mag_range=(21., 25.5),
                    dt_range=(0.5, 30.),
                    median_log_summary=median_log_summary,median_log_summaryb=None,
                    target=target,alpha=0.1,simu_name=simu_name)

#Plot(median_log_summary,bands,'m5','moonPhase')
Plot_Stat(median_log_summary,bands,'cadence',simu_name)
pl.show()

