import numpy as np
from Observations import *
import pylab as plt
from matplotlib import cm
import sncosmo
from Telescope import *
from astropy import (cosmology, units as u, constants as const)
import pickle as pkl

def median(sel,var):
    
    if var != '':
        selb=np.sort(sel,order=var)
        selb=selb[var]
    else:
        selb=np.sort(sel)

    num=len(selb)

    
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

def Plot_Obs_per_Field(resu):
    fontsize=12.
    
    myls=['-','--']
    colors=dict(zip([i for i in range(10)],['k','k','r','r','b','b','g','g','m','m']))
    
    for fieldid in fieldids:
        figa, axa = plt.subplots(ncols=1, nrows=3, figsize=(10,9))
        idx = resu['fieldid']==fieldid
        figa.suptitle(fieldname+' field - '+str(fieldid))
        tot_label=[]
        sela=resu[idx]
        m5_var={}
        for season in range(10):
            idxb = sela['season']==season
            selb=sela[idxb]
            x=[]
            y=[]
            z=[]
            w=[]
            
            for b in bands:
                if not m5_var.has_key(b):
                    m5_var[b]=[]
                idxc = selb['band']==b
                selc=selb[idxc]
                x.append(selc['ib'][0])
                y.append(selc['airmass'][0])
                z.append(selc['m5'][0])
                w.append(selc['seeing'][0])
                m5_var[b].append(selc['m5'][0])
            print fieldid,season+1,z
            axa[0].plot(x,y,ls=myls[season%2],color=colors[season])
            
            ll='Y'+str(season+1)
            tot_label.append(axa[1].errorbar(x,z,ls=myls[season%2],color=colors[season],label=ll))
                
            axa[2].plot(x,w,ls=myls[season%2],color=colors[season])
        

        axa[0].set_ylabel('Median airmass',{'fontsize': fontsize})
        axa[1].set_ylabel('Median m5 [mag]',{'fontsize': fontsize})
        
        labs = [l.get_label() for l in tot_label]
        axa[1].legend(tot_label, labs, ncol=5,loc='best',prop={'size':12},frameon=False)
        
        axa[2].set_ylabel('Median seeing [\'\']',{'fontsize': fontsize}) 

        print fieldid,'m5 variations'
        for b in bands:
            print b,np.max(m5_var[b])-np.min(m5_var[b])
        
        
        for j in range(3):
            axa[j].set_xlabel('band',{'fontsize': fontsize})
            axa[j].set_xlim([-0.1,4.1])
            axa[j].set_xticks([i for i in range(len(bands))])
            axa[j].set_xticklabels([corresp_inverted[i] for i in range(len(bands))])
        plt.gcf().savefig('Obs_Plots/'+fieldname+'_'+str(fieldid)+'.png')

def Plot_Cadence_per_Field(resu):
    fontsize=12.
    
    myls=['-','--']
    colors=dict(zip([i for i in range(10)],['k','k','r','r','b','b','g','g','m','m']))
    
    for fieldid in fieldids:
        figa, axa = plt.subplots(ncols=1, nrows=3, figsize=(10,9))
        idx = resu['fieldid']==fieldid
        figa.suptitle(fieldname+' field - '+str(fieldid))
        tot_label=[]
        sela=resu[idx]
        for season in range(10):
            idxb = sela['season']==season
            selb=sela[idxb]
            x=[]
            y=[]
            z=[]
            w=[]
            for b in bands:
                idxc = selb['band']==b
                selc=selb[idxc]
                x.append(selc['ib'][0])
                y.append(selc['mean_cadence'][0])
                z.append(selc['rms_cadence'][0])
                w.append(selc['duration'][0])
                print fieldid,season,x,y
            axa[0].plot(x,y,ls=myls[season%2],color=colors[season])
            
            ll='Y'+str(season+1)
            tot_label.append(axa[1].errorbar(x,z,ls=myls[season%2],color=colors[season],label=ll))
                
            axa[2].plot(x,w,ls=myls[season%2],color=colors[season])
        

        axa[0].set_ylabel('Mean cadence [day]',{'fontsize': fontsize})
        axa[1].set_ylabel('RMS cadence [day]',{'fontsize': fontsize})
        
        labs = [l.get_label() for l in tot_label]
        axa[1].legend(tot_label, labs, ncol=5,loc='best',prop={'size':12},frameon=False)
        
        #axa[2].set_ylabel('Duration [days]',{'fontsize': fontsize}) 
        axa[2].set_ylabel('Observation Period [day]',{'fontsize': fontsize}) 
        
        for j in range(3):
            axa[j].set_xlabel('band',{'fontsize': fontsize})
            axa[j].set_xlim([-0.1,4.1])
            axa[j].set_xticks([i for i in range(len(bands))])
            axa[j].set_xticklabels([corresp_inverted[i] for i in range(len(bands))])
        plt.gcf().savefig('Cadence_Plots/'+fieldname+'_'+str(fieldid)+'.png')

def Plot_per_Field(resu,what=('mean_cadence','rms_cadence'),myleg=('Mean cadence [day]','RMS cadence [day]')):
    fontsize=12.
    
    myls=['-','--']
    colors=dict(zip([i for i in range(10)],['k','k','r','r','b','b','g','g','m','m']))
    
    for fieldid in fieldids:
        figa, axa = plt.subplots(ncols=1, nrows=2, figsize=(10,9))
        idx = resu['fieldid']==fieldid
        figa.suptitle(fieldname+' field - '+str(fieldid))
        tot_label=[]
        sela=resu[idx]
        for season in range(10):
            idxb = sela['season']==season
            selb=sela[idxb]
            
            x=[]
            y=[]
            z=[]
            w=[]
            for b in bands:
                idxc = selb['band']==b
                selc=selb[idxc]
                x.append(selc['ib'][0])
                y.append(selc[what[0]][0])
                z.append(selc[what[1]][0])
                
            
            print fieldid,season,what[0],y,what[1],z,np.max(selb[what[0]])-np.min(selb[what[0]])
            axa[0].plot(x,y,ls=myls[season%2],color=colors[season])
            
            ll='Y'+str(season+1)
            tot_label.append(axa[1].errorbar(x,z,ls=myls[season%2],color=colors[season],label=ll))
                
            
        idxc = sela['band']=='a'
        print fieldid,season,sela['duration'][idxc]
        
        ax2 = axa[0].twiny()
        ax2.plot(sela['season'][idxc]+1,sela['duration'][idxc],ls='-',color='k',marker='s',label='grizy')
        
        axa[0].set_ylabel(myleg[0],{'fontsize': fontsize})
        axa[1].set_ylabel(myleg[1],{'fontsize': fontsize})
        ax2.set_xlabel('Year',{'fontsize': fontsize})

        labs = [l.get_label() for l in tot_label]
        axa[1].legend(tot_label, labs, ncol=5,loc='best',prop={'size':12},frameon=False)
        ax2.legend(loc='best',prop={'size':12})
        for j in range(2):
            axa[j].set_xlabel('band',{'fontsize': fontsize})
            axa[j].set_xlim([-0.1,4.1])
            axa[j].set_xticks([i for i in range(len(bands))])
            axa[j].set_xticklabels([corresp_inverted[i] for i in range(len(bands))])

        
        #plt.gcf().savefig('Cadence_Plots/'+fieldname+'_'+str(fieldid)+'.png')

def Plot_Cadence_per_Year(resu,nyears=10):

    #myls=['-','--']
    colors=dict(zip(fieldids,['k','r','b','g','m']))
    fontsize=12

    for season in range(nyears):
        idx = resu['season']==season
        sela=resu[idx]
        figa, axa = plt.subplots(ncols=1, nrows=3, figsize=(10,9))
        figa.suptitle(fieldname+' - Year '+str(season+1))
        tot_label=[]
        for fieldid in fieldids:
            idxb = sela['fieldid'] == fieldid
            selb=sela[idxb]
            x=[]
            y=[]
            z=[]
            w=[]
            for b in bands:
                idxc = selb['band']==b
                selc=selb[idxc]
                x.append(selc['ib'][0])
                y.append(selc['mean_cadence'][0])
                z.append(selc['rms_cadence'][0])
                w.append(selc['duration'][0])
                #print fieldid,season,x,y
            axa[0].plot(x,y,ls='-',color=colors[fieldid])
            axa[0].set_xlabel('band',{'fontsize': fontsize})
            axa[0].set_ylabel('Mean cadence [day]',{'fontsize': fontsize})
            axa[0].set_xlim([-0.1,4.1])
            axa[0].set_xticks([i for i in range(len(bands))])
            axa[0].set_xticklabels([corresp_inverted[i] for i in range(len(bands))])
            ll=str(fieldid)
            tot_label.append(axa[1].errorbar(x,z,ls='-',color=colors[fieldid],label=ll))
                
            axa[2].plot(x,w,ls='-',color=colors[fieldid])


        axa[0].set_ylabel('Mean cadence [day]',{'fontsize': fontsize})
        axa[1].set_ylabel('RMS cadence [day]',{'fontsize': fontsize})
        
        labs = [l.get_label() for l in tot_label]
        axa[1].legend(tot_label, labs, ncol=5,loc='best',prop={'size':12},frameon=False)
        
        #axa[2].set_ylabel('Duration [days]',{'fontsize': fontsize}) 
        axa[2].set_ylabel('Season length [day]',{'fontsize': fontsize})
        
        for j in range(3):
            axa[j].set_xlabel('band',{'fontsize': fontsize})
            axa[j].set_xlim([-0.1,4.1])
            axa[j].set_xticks([i for i in range(len(bands))])
            axa[j].set_xticklabels([corresp_inverted[i] for i in range(len(bands))])
        
        plt.gcf().savefig('Cadence_Plots/'+fieldname+'_Y'+str(season+1)+'.png')

def Plot_Cadence_vs_Year(resu,fieldname,fieldid,nyears=10):

    #myls=['-','--']
    colors=dict(zip(bands,['b','c','g','y','r','m']))
    fontsize=12

    res=[]

    for season in range(nyears):
        idx = (resu['season']==season)&(resu['fieldid']==fieldid)
        sela=resu[idx]
        for b in bands+'a':
            idxc = sela['band']==b
            selc=sela[idxc]
            #print 'aouh',season, b,fieldid,len(selc),selc['median_cadence']
            if len(selc) >0:
                res.append((season+1,b,np.asscalar(selc['median_cadence']),np.asscalar(selc['upper_cadence']),np.asscalar(selc['lower_cadence'])))
    
       
    resu=np.rec.fromrecords(res,names=['season','band','median_cadence','upper_cadence','lower_cadence'])

    f = open('Cadence_tex/Cadence_'+fieldname+'_'+str(fieldid)+'.tex','w')
    
    f.write('\\begin{table}[ht]'+'\n')
    f.write('\\begin{center}'+'\n')
    f.write('\\caption{'+fieldname+' '+str(fieldid)+'. Median cadence ($\pm$ 95\% confidence intervals) (in day$^{-1}$) for each band and each season. }\label{tab:cad'+str(fieldid)+'}'+'\n')
    
    f.write('\\begin{tabular}{ccccccc}'+'\n')
    f.write('\\hline'+'\n')

    chap='season'
    chap+=' & \multicolumn{6}{c}{band} \\\\'
    f.write(chap+'\n')
    chap=' '
    for b in bands:
        chap+=' & '+b
    chap+=' & all'
    chap+=' \\\\'
    f.write(chap+'\n')
    print chap
    f.write('\\hline'+'\n')
    for seas in range(1,11):
        idx=resu['season']==seas
        resb=resu[idx]
        toprint=str(seas)
        for b in bands+'a':
            idd=resb['band']==b
            resc=resb[idd]
            if len(resc) > 0:
                med=str(round(np.asscalar(resc['median_cadence']),1))
                upp=str(round(np.asscalar(resc['upper_cadence']-resc['median_cadence']),1))
                low=str(round(np.asscalar(resc['median_cadence']-resc['lower_cadence']),1))
                if med == '0.0' and upp == '0.0' and low == '0.0':
                    toprint+=' & - '
                else:
                    toprint+=' & '+med+'$^{+'+str(upp)+'}_{-'+str(low)+'}$'
            else:
                med='-'
                toprint+=' & '+med
            #toprint+=' & '+med+'$^{+'+str(upp)+'}_{-'+str(low)+'}$'
        toprint+=' \\\\'
        print toprint
        f.write(toprint+' \n')
    
    f.write('\\hline'+'\n')
    f.write('\\end{tabular}'+'\n')
    f.write('\\end{center}'+'\n')
    f.write('\\end{table}'+'\n')


    figa, axa = plt.subplots(ncols=1, nrows=3, figsize=(10,9))
    figa.suptitle(fieldname+' - '+str(fieldid))

    tot_label=[]
    for b in bands:
        ixc = resu['band']==b
        resb=resu[ixc]
        axa[0].plot(resb['season'],resb['median_cadence'],ls='-',color=colors[b])
        tot_label.append(axa[1].errorbar(resb['season'],resb['upper_cadence']-resb['median_cadence'],ls='-',color=colors[b],label=b+' band'))
        axa[2].plot(resb['season'],resb['median_cadence']-resb['lower_cadence'],ls='-',color=colors[b])
        
        #print b,resb['median_cadence']
        #print b,resb['upper_cadence']
        #print b,resb['lower_cadence']

   

    for i in range(3):
        axa[i].set_xlabel('Year',{'fontsize': fontsize})

    axa[0].set_ylabel('Median cadence [day$^{-1}$]')
    axa[1].set_ylabel('Upper-median cadence [day$^{-1}$]')
    axa[2].set_ylabel('Median-lower cadence [day$^{-1}$]')

    labs = [l.get_label() for l in tot_label]
    axa[1].legend(tot_label, labs, ncol=6,loc='upper center',prop={'size':10},frameon=False)


def Plot_Nobs(Nobs):

    for season in range(1):
        Plot_Nobs_Indiv(Nobs[Nobs['season']==season],season)

def Plot_Nobs_Indiv(Nobs,season):
 
   for band in ['g','r','i','z','y','all']:    
        figb, axb = plt.subplots(ncols=1, nrows=2, figsize=(10,9))
            #figb.suptitle(band+' band')
        figb.suptitle(fieldname+' - '+str(fieldid)+' - Year '+str(season+1)+' '+band+' band') 
        idx = Nobs['band']==band
            #sel=Nobs[np.where(np.logical_and(Nobs['T0']>=62.,Nobs['T0']<63.))]
        sel=Nobs[idx]
        ll=''
        if band == 'all':
            ll='grizy'
            
        axb[0].plot(sel['T0'],sel['Nbef'],label=ll)
        axb[1].plot(sel['T0'],sel['Naft'],label=ll)

        if band == 'all':
            idx = Nobs['band']==band+'_no_g'
            #sel=Nobs[np.where(np.logical_and(Nobs['T0']>=62.,Nobs['T0']<63.))]
            sel=Nobs[idx]
            ll='rizy'
            axb[0].plot(sel['T0'],sel['Nbef'],label=ll,color='r')
            axb[1].plot(sel['T0'],sel['Naft'],label=ll,color='r')
            axb[0].plot(sel['T0'],[4.]*len(sel['T0']),color='k')
            axb[1].plot(sel['T0'],[10.]*len(sel['T0']),color='k')

        axb[0].set_xlabel('T0 [day]',{'fontsize': fontsize})
        axb[0].set_ylabel('Nobs in [T0-20, T0]',{'fontsize': fontsize})
        axb[1].set_xlabel('T0 [day]',{'fontsize': fontsize})
        axb[1].set_ylabel('Nobs in [T0, T0+40]',{'fontsize': fontsize})

        if band =='all':
            axb[0].legend(loc='upper left',prop={'size':fontsize})
            axb[1].legend(loc='upper left',prop={'size':fontsize})

def Plot_Diffs(diffs,fieldname,fieldids,colors):
    
        fontsize=12
        r=[]
        for fieldid in fieldids:
            for key,diff in diffs[fieldid].items():
            #print diff.keys()
                
                figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
                figa.suptitle(fieldname+' - '+str(fieldid)+' - Year '+str(key+1))
                
                #for band in 'grizy':
                for band in 'z':
                    ll=band+' band'
                    axa.errorbar([val for val in diff[band][0]],[val for val in diff[band][1]],color=filtercolors[band],marker='o',label=ll)
                    axa.errorbar([val for val in diff[band][0]],[np.median([val for val in diff[band][1]])]*len(diff[band][0]),color=filtercolors[band],ls='--')
                    totag=sorted(diff[band][1])
                    
                    median=np.median(totag)
                    nvals=float(len(totag))
                    
                    rank_cl_upper=1.+nvals/2.+1.96*np.sqrt(nvals)/2.
                    rank_cl_lower=nvals/2.-1.96*np.sqrt(nvals)/2.
                    #print 'alors',key,band,nvals,rank_cl_upper,rank_cl_lower
                    median_cl_upper=totag[int(rank_cl_upper)]
                    median_cl_lower=totag[int(rank_cl_lower)]
                    
                    r.append((fieldid,key+1,band,median,median_cl_lower,median_cl_upper))
                axa.legend(loc='best',prop={'size':fontsize})
                axa.set_ylabel('$\Delta$T = T$_{obs}$-T$_{obs-1}$ [day]',{'fontsize': fontsize})
                axa.set_xlabel('MJD [day]',{'fontsize': fontsize})
                
                #plt.gcf().savefig('Obs_Plots/DeltaT_'+fieldname+'_'+str(fieldid)+'_Y'+str(key+1)+'.png')
                #plt.close(figa)
        med_diffs=np.rec.fromrecords(r,names=['fieldid','season','band','median_diff','median_diff_lower','median_diff_upper'])

        myband='z'
        figb, axb = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        figb.suptitle(myband+' band')
        if len(fieldids) ==1:
            figb.suptitle('Field '+str(fieldids[0])+' '+myband+' band')
        for fieldid in fieldids:
            idx=(med_diffs['fieldid']==fieldid)&(med_diffs['band']==myband)
            sela=med_diffs[idx]
            print fieldid
            
            if len(fieldids) == 1:
                ll='median'
                axb.plot(sela['season'],sela['median_diff'],color=colors[fieldid],label=ll)
                ll='95% CL (lower)'
                axb.plot(sela['season'],sela['median_diff_lower'],color=colors[fieldid],ls='--',label=ll)
                ll='95% CL (upper)'
                axb.plot(sela['season'],sela['median_diff_upper'],color=colors[fieldid],ls=':',label=ll)
            else:
                ll='Field '+str(fieldid)
                axb.plot(sela['season'],sela['median_diff'],color=colors[fieldid],label=ll)
            for season in range(1,11):
                selb=sela[np.where(sela['season']==season)]
                print season,selb['median_diff'][0],selb['median_diff_lower'][0],selb['median_diff_upper'][0]
        axb.set_xlabel('Year',{'fontsize': fontsize})
        axb.set_ylabel('Median $\Delta$T [day]',{'fontsize': fontsize})
        axb.legend(loc='best',prop={'size':fontsize})
        ylim = axb.get_ylim()
        axb.set_ylim([ylim[0]-0.1,ylim[1]])
        xlim = axb.get_xlim()
        axb.set_xlim([xlim[0]-0.1,xlim[1]+0.1])
        major_ticks = np.arange(ylim[0],ylim[1], 1)
        axb.set_yticks(major_ticks)
        major_ticks = np.arange(xlim[0],xlim[1]+1,1)
        axb.set_xticks(major_ticks)
        axb.grid(which='both')

        figc, axc = plt.subplots(ncols=1, nrows=2, figsize=(10,9))
        figc.suptitle(myband+' band')
        tot_label=[]
        
        for fieldid in fieldids:
            idx=(med_diffs['fieldid']==fieldid)&(med_diffs['band']==myband)
            sela=med_diffs[idx]
           
            ll='Field '+str(fieldid)
            tot_label.append(axc[0].errorbar(sela['season'],sela['median_diff']-sela['median_diff_lower'],color=colors[fieldid],label=ll))
            axc[1].plot(sela['season'],sela['median_diff_upper']-sela['median_diff'],color=colors[fieldid],label=ll)

        axc[0].set_xlabel('Year',{'fontsize': fontsize})
        axc[0].set_ylabel('Median $\Delta$T - Lower (95% C.L.) $\Delta$T [day]',{'fontsize': fontsize})
        axc[0].legend(loc='best',prop={'size':fontsize})
        ylim = axc[0].get_ylim()
        axc[0].set_ylim([ylim[0]-0.5,ylim[1]])
        xlim = axc[0].get_xlim()
        axc[0].set_xlim([xlim[0]-0.1,xlim[1]+0.1])
        major_ticks = np.arange(ylim[0],ylim[1], 1)
        axc[0].set_yticks(major_ticks)
        major_ticks = np.arange(xlim[0],xlim[1]+1,1)
        axc[0].set_xticks(major_ticks)
        axc[0].grid(which='both') 
        labs = [l.get_label() for l in tot_label]
        axc[0].legend(tot_label, labs, ncol=5,loc='best',prop={'size':10},frameon=False)

        axc[1].set_xlabel('Year',{'fontsize': fontsize})
        axc[1].set_ylabel(' Upper (95% C.L.) $\Delta$T - Median $\Delta$T [day]',{'fontsize': fontsize})
        axc[1].legend(loc='best',prop={'size':fontsize})
        ylim = axc[1].get_ylim()
        axc[1].set_ylim([ylim[0]-0.5,ylim[1]])
        xlim = axc[1].get_xlim()
        axc[1].set_xlim([xlim[0]-0.1,xlim[1]+0.1])
        major_ticks = np.arange(ylim[0],ylim[1], 1)
        axc[1].set_yticks(major_ticks)
        major_ticks = np.arange(xlim[0],xlim[1]+1,1)
        axc[1].set_xticks(major_ticks)
        axc[1].grid(which='both')
        axc[1].legend(tot_label, labs, ncol=5,loc='best',prop={'size':10},frameon=False)

def Plot_median_m5(res,fieldname,fieldids,filtercolors,what='m5',leg='m$_5$'):


    lleg=['Median - Lower (95% C.L.) '+leg+' [mag]','Upper (95% C.L.) -median '+leg+'[mag]']
    fontsize=10

    for fieldid in fieldids:
        idx = res['fieldid']==fieldid
        sela=res[idx]
        ras=[]
        figa, axa = plt.subplots(ncols=1, nrows=2, figsize=(10,9))
        figa.suptitle(fieldname+' - '+str(fieldid))
        figb, axb = plt.subplots(ncols=1, nrows=2, figsize=(10,9))
        figb.suptitle(fieldname+' - '+str(fieldid))
        tot_label=[]
        tot_labelb=[]
        tot_labelc=[]
        for band in 'ugrizy':
            idxb = sela['band']==band
            selb=sela[idxb]
            median_m5=np.median(selb[what])
            tot_label.append(axa[0].errorbar(selb['season']+1,selb[what],color=filtercolors[band],label=band+' band'))
            tot_labelc.append(axa[1].errorbar(selb['season']+1,selb[what]-median_m5,color=filtercolors[band],label=band+' band'))
            tot_labelb.append(axb[0].errorbar(selb['season']+1,selb[what]-selb[what+'_lower'],color=filtercolors[band],label=band+' band'))
            axb[1].errorbar(selb['season']+1,selb[what+'_upper']-selb[what],color=filtercolors[band],label=band+' band')

        axa[0].set_xlabel('Year',{'fontsize': fontsize})
        axa[0].set_ylabel(leg+'$^{median}$ [mag]',{'fontsize': fontsize})
        labs = [l.get_label() for l in tot_label]
        axa[0].legend(tot_label, labs, ncol=5,loc='best',prop={'size':12},frameon=False)

        axa[1].set_xlim([0.8,10.2])
        axa[1].set_xlabel('Year',{'fontsize': fontsize})
        axa[1].set_ylabel(leg+'$^{median}$(Year)-m$_{5}$$^{median}$ [mag]',{'fontsize': fontsize})
        labs = [l.get_label() for l in tot_labelc]
        axa[1].legend(tot_labelc, labs, ncol=5,loc='best',prop={'size':12},frameon=False)
        axa[1].set_xticks([i for i in range(1,11,1)])
        """
        plt.gcf().savefig('Obs_Plots/'+what+'_med_year_'+fieldname+'_'+str(fieldid)+'.png')
        plt.close(figa)
        """
        
        for i in range(2):
            axb[i].set_xlim([0.8,10.2])
            axb[i].set_xlabel('Year',{'fontsize': fontsize})
            axb[i].set_ylabel(lleg[i],{'fontsize': fontsize})
            labs = [l.get_label() for l in tot_labelb]
            axb[i].legend(tot_labelb, labs, ncol=5,loc='best',prop={'size':12},frameon=False)
            axb[i].set_xticks([i for i in range(1,11,1)])

        #plt.gcf().savefig('Obs_Plots/'+what+'_upper_lower_year_'+fieldname+'_'+str(fieldid)+'.png')
        #plt.close(figa)
        

def Print_begin(f,fieldname,fieldid):
    
    f.write('\\begin{table}[ht]'+'\n')
    f.write('\\begin{center}'+'\n')
    f.write('\\caption{'+fieldname+' '+str(fieldid)+'. Total number of visits per year.}\label{tab:visit'+str(fieldid)+'}\n')
    
    f.write('\\begin{tabular}{cccccccc}'+'\n')
    f.write('\\hline'+'\n')
    
    chap='Year'
    for band in bands:
        chap+=' & '+band
    chap+= ' & all \\\\'
    f.write(chap+'\n')
    f.write('\\hline'+'\n')
        
def Print_begin_multi(f,fieldname):
    
    f.write('\\begin{table}[ht]'+'\n')
    f.write('\\begin{center}'+'\n')
    f.write('\\caption{'+fieldname+' fields. Median 5-$\sigma$ depth and sky magnitudes.}\label{tab:'+fieldname+'_median_mag}\n')
    
    f.write('\\begin{tabular}{ccccccccccccc}'+'\n')
    f.write('\\hline'+'\n')
    f.write(' & \\multicolumn{6}{c}{5-$\sigma$ depth [mag]} & \\multicolumn{6}{c}{m$_{sky}$ [mag]} \\\\'+'\n')
    chap='Year'
    for band in bands:
        chap+=' & '+band
    for band in bands:
        chap+=' & '+band
    chap+= ' & \\\\'
    f.write(chap+'\n')
    f.write('\\hline'+'\n')

def Print_end(f):
    
    f.write('\\hline'+'\n')
    f.write('\\end{tabular}'+'\n')
    f.write('\\end{center}'+'\n')
    f.write('\\end{table}'+'\n')   

def Plot_median(res,fieldids,fieldcolors,nyears=10):

    figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
    figb, axb = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
    fontsize=12

    print 'there pal',fieldids
    for fieldid in fieldids:
        idx = res['fieldid']==fieldid
        sela=res[idx]
        ras=[]
        
        
        for band in 'ugrizy':
            idxb = sela['band']==band
            selb=sela[idxb]
            print 'hello',fieldid,selb,selb.dtype
            median_len, lower_len,upper_len= median(selb,'duration')
            median_obs, lower_obs,upper_obs= median(selb,'obstime')
            ras.append((np.mean(selb['ib']),median_len, lower_len,upper_len,median_obs/3600., lower_obs/3600.,upper_obs/3600.))
        
        rasrec=np.rec.fromrecords(ras,names=['iband','duration','duration_lower','duration_upper','obstime','obstime_lower','obstime_upper'])
        ll='Fieldid '+str(fieldid)
        #axa[0].plot([vv[0] for vv in ras],[vv[1] for vv in ras],color=fieldcolors[fieldid],label=ll)
        axa.plot(rasrec['iband'],rasrec['duration'],color=fieldcolors[fieldid],label=ll)
        axa.plot(rasrec['iband'],rasrec['duration_upper'],color=fieldcolors[fieldid],ls='--')
        axa.plot(rasrec['iband'],rasrec['duration_lower'],color=fieldcolors[fieldid],ls='--')

        #axa[1].plot([vv[0] for vv in ras],[vv[2] for vv in ras],color=fieldcolors[fieldid],label=ll)
        axb.plot(rasrec['iband'],rasrec['obstime'],color=fieldcolors[fieldid],label=ll)
        axb.plot(rasrec['iband'],rasrec['obstime_upper'],color=fieldcolors[fieldid],ls='--')
        axb.plot(rasrec['iband'],rasrec['obstime_lower'],color=fieldcolors[fieldid],ls='--')

    
    #axa[0].set_ylabel('Median duration [day] ',{'fontsize': fontsize})
    axa.set_ylabel('Median season length [day] ',{'fontsize': fontsize})
    axa.legend(loc='best',prop={'size':fontsize},frameon=False)
    axa.set_xlabel('band',{'fontsize': fontsize})
    axa.set_xlim([-0.1,4.1])
    axa.set_xticks([i for i in range(len(bands))])
    axa.set_xticklabels([corresp_inverted[i] for i in range(len(bands))])

    axb.set_ylabel('Median expTime [h]',{'fontsize': fontsize})
    axb.legend(loc='best',prop={'size':fontsize},frameon=False)
    
    axb.set_xlabel('band',{'fontsize': fontsize})
    axb.set_xlim([-0.1,4.1])
    axb.set_xticks([i for i in range(len(bands))])
    axb.set_xticklabels([corresp_inverted[i] for i in range(len(bands))])
    figa.savefig('Obs_Plots/season_length_'+fieldname+'_'+str(fieldid)+'.png')
    figb.savefig('Obs_Plots/exptime_'+fieldname+'_'+str(fieldid)+'.png')

    for band in 'grizy':
        figc, axc = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        figc.suptitle(band+' band')
        for fieldid in fieldids:
            ll='Fieldid '+str(fieldid)
            idx = res['fieldid']==fieldid
            sela=res[idx]
            idxb = sela['band']==band
            sela=sela[idxb]
            axc.plot(sela['season']+1,sela['duration'],color=fieldcolors[fieldid],label=ll)

        axc.set_ylabel('Season length [day]',{'fontsize': fontsize})
        axc.set_xlabel('Year',{'fontsize': fontsize})
        axc.legend(loc='best',prop={'size':fontsize},frameon=False)

    
    figc, axc = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
    r=[]
    for fieldid in fieldids:
        ll='DD '+str(fieldid)
        idx = res['fieldid']==fieldid
        sela=res[idx]
        idxb = (sela['band'] != 'u')&(sela['band'] != 'a')
        sela=sela[idxb]
        axc.plot(sela['season']+1,sela['duration'],color=fieldcolors[fieldid],label=ll)

        axc.set_ylabel('Season length [day]',{'fontsize': fontsize})
        axc.set_xlabel('Year',{'fontsize': fontsize})
        axc.legend(loc='best',prop={'size':fontsize},frameon=False)
        med,low,upp=median(sela,'duration')
        r.append((fieldid,med,low,upp))
    plt.gcf().savefig('DD_season_length.png')
        

    med_len=np.rec.fromrecords(r,names=['fieldid','median_len','lower_len','upper_len'])  

    
    #put this is a tex file
    f = open('Cadence_tex/DD_Season_Length.tex','w')
    
    f.write('\\begin{table}[ht]'+'\n')
    f.write('\\begin{center}'+'\n')
    f.write('\\caption{Season length ($\pm$ 95\% confidence intervals) (in day) for the DD fields. All bands but u were considered.}\label{tab:seasonlength'+str(fieldid)+'}\n')
    
    f.write('\\begin{tabular}{cc}'+'\n')
    f.write('\\hline'+'\n')

    chap='Field'
    chap+=' & Season length (day)\\\\'
    f.write(chap+'\n')
    f.write('\\hline'+'\n')
    for field in fieldids:
        toprint=''
        idx= med_len['fieldid']==field
        med_lenf=med_len[idx]
        med=str(round(np.asscalar(med_lenf['median_len']),1))
        upp=str(round(np.asscalar(med_lenf['upper_len']-med_lenf['median_len']),1))
        low=str(round(np.asscalar(med_lenf['median_len']-med_lenf['lower_len']),1))
        toprint+=str(field)+' & '+med+'$^{+'+str(upp)+'}_{-'+str(low)+'}$'
        toprint+=' \\\\'
        print toprint
        f.write(toprint+' \n')
    
    f.write('\\hline'+'\n')
    f.write('\\end{tabular}'+'\n')
    f.write('\\end{center}'+'\n')
    f.write('\\end{table}'+'\n')    


    for fieldid in fieldids:
        idx = res['fieldid'] == fieldid
        sela=res[idx]
        f = open('Cadence_tex/DD_'+str(fieldid)+'_Nvisits.tex','w')
        Print_begin(f,'DD',fieldid)
        for season in range(nyears):
            idxb = sela['season']==season
            selb=sela[idxb]
            toprint=str(season+1)
            for b in bands+'a':
                idxc = selb['band']==b
                selc= selb[idxc]
                #print fieldid,season,b,selc['nexp_sum']
                toprint+=' & '+str(int(np.asscalar(selc['nexp_sum'])))
            toprint+=' \\\\'
            f.write(toprint+' \n')
        Print_end(f)
        f.close()
        #break
    
   

   



def Plot_Observations(obs,fieldname,fieldid):

    fontsize=12

    figc, axc = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
    figc.suptitle(fieldname+' field '+str(fieldid))
    axc.scatter(obs.all_seasons['mjd'],obs.all_seasons['m5sigmadepth'],s=80,facecolors='none', edgecolors='r')
    axc.set_xlim([np.min(obs.all_seasons['mjd'])-50.,np.max(obs.all_seasons['mjd'])+50.])
    axc.set_xlabel('MJD [day]',{'fontsize': fontsize})
    axc.set_ylabel('$m_{5\sigma}$ [mag]',{'fontsize': fontsize})
    plt.gcf().savefig('Obs_Plots/m5_vs_mjd_'+fieldname+'_'+str(fieldid)+'.png')

    for i in range(len(obs.seasons)):
        figc, axc = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        figc.suptitle(fieldname+' field '+str(fieldid)+'- Season '+str(i+1))
        axc.scatter(obs.all_seasons['mjd'],obs.all_seasons['m5sigmadepth'],s=80,facecolors='none', edgecolors='r')
        axc.plot(obs.seasons[i]['mjd'],obs.seasons[i]['m5sigmadepth'],'*r')
        axc.set_xlim([np.min(obs.all_seasons['mjd'])-50.,np.max(obs.all_seasons['mjd'])+50.])
        axc.set_xlabel('MJD [day]',{'fontsize': fontsize})
        axc.set_ylabel('$m_{5\sigma}$ [mag]',{'fontsize': fontsize})
        plt.gcf().savefig('Obs_Plots/m5_vs_mjd_'+fieldname+'_'+str(fieldid)+'_season_'+str(i+1)+'.png')

    
    Plot_Obs_per_filter(data=obs.seasons[0],fieldname=fieldname,fieldid=fieldid,season=1)

def Plot_Obs_per_filter(data=None, fieldname='',fieldid=0,season=0, flux_name='m5sigmadepth',xfigsize=None, yfigsize=None, figtext=None,figtextsize=1.,ncol=2,color=None,cmap=None, cmap_lims=(3000., 10000.)):
       
    telescope=Telescope(airmass=np.median(1.2))

    transmission=telescope.throughputs
    
    for filtre in 'ugrizy':
        band=sncosmo.Bandpass(transmission.atmosphere[filtre].wavelen,transmission.atmosphere[filtre].sb, name='LSSTPG::'+filtre,wave_unit=u.nm)
        sncosmo.registry.register(band)

    toff=0.
    bands=set(data['band'])
    
    #tmin, tmax = [], []
    if data is not None:
        tmin=np.min(data['mjd']) - 10.
        tmax=np.max(data['mjd']) + 10.
        #tmin.append(np.min(data['mjd']) - 10.)
        #tmax.append(np.max(data['mjd']) + 10.)
    print 'hello',tmin,tmax

        # Calculate layout of figure (columns, rows, figure size). We have to
        # calculate these explicitly because plt.tight_layout() doesn't space the
        # subplots as we'd like them when only some of them have xlabels/xticks.
    wspace = 0.6  # All in inches.
    hspace = 0.3
    lspace = 1.0
    bspace = 0.7
    trspace = 0.2
    nrow = (len(bands) - 1) // ncol + 1
    if xfigsize is None and yfigsize is None:
        hpanel = 2.25
        wpanel = 3.
    elif xfigsize is None:
        hpanel = (yfigsize - figtextsize - bspace - trspace -
                  hspace * (nrow - 1)) / nrow
        wpanel = hpanel * 3. / 2.25
    elif yfigsize is None:
        wpanel = (xfigsize - lspace - trspace - wspace * (ncol - 1)) / ncol
        hpanel = wpanel * 2.25 / 3.
    else:
        raise ValueError('cannot specify both xfigsize and yfigsize')

    figsize = (lspace + wpanel * ncol + wspace * (ncol - 1) + trspace,
               bspace + hpanel * nrow + hspace * (nrow - 1) + trspace +
               figtextsize)  

        # Create the figure and axes.
    fig, axes = plt.subplots(nrow, ncol, figsize=figsize, squeeze=False)
    fig.suptitle(fieldname+' Field '+str(fieldid)+' - Season '+str(season))

    fig.subplots_adjust(left=lspace / figsize[0],
                        bottom=bspace / figsize[1],
                        right=1. - trspace / figsize[0],
                        top=1. - (figtextsize + trspace) / figsize[1],
                        wspace=wspace / wpanel,
                        hspace=hspace / hpanel)
        # Color options.
    if color is None:
        if cmap is None:
            cmap = cm.get_cmap('jet_r')
        # Loop over bands
    bands = list(bands)
    waves = [sncosmo.get_bandpass(b).wave_eff for b in bands]
    waves_and_bands = sorted(zip(waves, bands))

    for axnum in range(ncol * nrow):
        row = axnum // ncol
        col = axnum % ncol
        ax = axes[row, col]
        
        if axnum >= len(waves_and_bands):
            ax.set_visible(False)
            ax.set_frame_on(False)
            continue

        wave, band = waves_and_bands[axnum]

        bandname_coords = (0.92, 0.92)
        bandname_ha = 'right'
        if color is None:
            bandcolor = cmap((cmap_lims[1] - wave) /
                                 (cmap_lims[1] - cmap_lims[0]))
        else:
            bandcolor = color

        if data is not None:
            mask = data['band'] == band
            time = data['mjd'][mask]
            flux = data[flux_name][mask]
            ax.errorbar(time, flux, ls='None',
                        color=bandcolor, marker='.', markersize=10.)    
            # Band name in corner
            ax.text(bandname_coords[0], bandname_coords[1], band,
                    color='k', ha=bandname_ha, va='top', transform=ax.transAxes)

            ax.axhline(y=0., ls='--', c='k')  # horizontal line at flux = 0.
            ax.set_xlim((tmin, tmax))
            ax.set_ylim((np.min(flux)-1.,np.max(flux)+1.))
            if (len(bands) - axnum - 1) < ncol:
                ax.set_xlabel('time')
               
            else:
                for l in ax.get_xticklabels():
                    l.set_visible(False)
            if col == 0:
                if flux_name == 'flux':
                    ax.set_ylabel('flux ($ZP_{{{0}}} = {1}$)'
                                  .format(sncosmo.get_magsystem(zpsys).name.upper(), zp))
                if flux_name == 'flux_e_sec':
                   ax.set_ylabel('flux (e/sec)')
                
                if flux_name == 'm5sigmadepth':
                    ax.set_ylabel('$m_{5\sigma}$ [mag]')

    plt.gcf().savefig('Obs_Plots/m5_vs_filter_'+fieldname+'_'+str(fieldid)+'_season_'+str(season)+'.png')

def Fill_Cadence(full_season,season):
    ra=[]
    min_season=np.min(full_season['mjd'])
    max_season=np.max(full_season['mjd'])
        
    for val in np.arange(min_season,max_season,1.):
        idxa = np.logical_and(full_season['mjd']> val -20. ,full_season['mjd'] < val)
        idxb = np.logical_and(full_season['mjd']> val,full_season['mjd'] < val+40.)
        sela=full_season[idxa]
        selb=full_season[idxb]
        ra.append((season+1,'all',val,len(sela),len(selb),len(sela)+len(selb)))
        ppa=sela[sela['band'] != 'LSSTPG::g']
        ppb=selb[selb['band'] != 'LSSTPG::g']
        ra.append((season+1,'all_no_g',val,len(ppa),len(ppb),len(ppa)+len(ppb)))
        #print val,len(sela),len(selb)
        for b in 'grizy':
            selbf=sela[np.where(sela['band']=='LSSTPG::'+b)]
            selaft=selb[np.where(selb['band']=='LSSTPG::'+b)]
            #print 'hello',b,selbf
            ra.append((season+1,b,val,len(selbf),len(selaft),len(selbf)+len(selaft)))
    return ra


def Fill_Medians(myseason,season):

    r=[]

    alldiffs={}
    for b in bands:
        idx = myseason['band']=='LSSTPG::'+b
        sel = myseason[idx]
            
        sel.sort(order='mjd')
        """
        if b == 'g':
            iseason+=1
        """
        diff=[io-jo for jo,io in zip(sel['mjd'][:-1], sel['mjd'][1:])]
        alldiffs[b]=(sel['mjd'][1:],diff)
        #print 'hello',b,diff,len(sel),sel

        if len(sel)>1:
            #print 'season',season
            m5_med,m5_lower,m5_upper=median(sel,'m5sigmadepth')
            msky_med,msky_lower,msky_upper=median(sel,'sky')
            airmass_med,airmass_lower,airmass_upper=median(sel,'airmass')
            seeing_med,seeing_lower,seeing_upper=median(sel,'seeing')
            median_cad,lower_cad,upper_cad=median(diff,'')
            #nexp_med,nexp_low,nexp_upper=median(sel,'Nexp')
            r.append((key,season,b,np.mean(diff),np.std(diff),np.max(sel['mjd'])-np.min(sel['mjd']),corresp[b],m5_med,m5_lower,m5_upper,np.sum(sel['exptime']),median_cad,upper_cad,lower_cad,np.sum(sel['Nexp']),msky_med,msky_lower,msky_upper,airmass_med,airmass_lower,airmass_upper,seeing_med,seeing_lower,seeing_upper,np.mean(sel['Ra']),np.mean(sel['Dec'])))

    idxc = myseason['band']!='LSSTPG::nn'
    selcb=myseason[idxc]
    
    selcb.sort(order='mjd')
    diff=[io-jo for jo,io in zip(selcb['mjd'][:-1], selcb['mjd'][1:])]
    if len(selcb) > 1:
        m5_med,m5_lower,m5_upper=median(selcb,'m5sigmadepth')
        msky_med,msky_lower,msky_upper=median(selcb,'sky')
        airmass_med,airmass_lower,airmass_upper=median(selcb,'airmass')
        median_cad,lower_cad,upper_cad=median(diff,'')
        seeing_med,seeing_lower,seeing_upper=median(sel,'seeing')
        #nexp_med,nexp_low,nexp_upper=median(sel,'Nexp')
        r.append((key,season,'a',np.mean(diff),np.std(diff),np.max(selcb['mjd'])-np.min(selcb['mjd']),corresp['a'],m5_med,m5_lower,m5_upper,np.sum(selcb['exptime']),median_cad,upper_cad,lower_cad,np.sum(selcb['Nexp']),msky_med,msky_lower,msky_upper,airmass_med,airmass_lower,airmass_upper,seeing_med,seeing_lower,seeing_upper,np.mean(selcb['Ra']),np.mean(selcb['Dec'])))

    return r,alldiffs


dirmeas='Mean_Obs_newrefs'
dirmeas='WFD'
fieldname='WFD'

thedir='OpSimLogs/'+dirmeas

myobs={}

#fieldids=[120,121,122,123]
#fieldids=[124,125,126,127]
#fieldids=[128,129,130,131]
#fieldids=[field+4 for field in fieldids]
#print 'alors',fieldids
fieldids=[310]
#fieldids=[1427]
#fieldids=[290,744,1427,2412,2786]
#fieldids=[744,1427]
#fieldids=[2786]

#fieldids=np.loadtxt('WFD_Ra_Dec.txt',dtype={'names': ('fieldid', 'ra', 'dec'),'formats': ('i4','f8','f8')})

#fieldids=fieldids['fieldid']

filelist='fieldIDs_minion_1016_'+fieldname+'.txt'
fields=np.loadtxt(filelist,dtype={'names': ('name','fieldid'),'formats': ('S8','i4')})

fieldids=fields['fieldid'][:]

for fieldid in fieldids:
    #print 'hello',fieldid
    name='Observations_'+fieldname+'_'+str(fieldid)+'.txt'
    myobs[fieldid]=Observations(fieldid=fieldid, filename=thedir+'/'+name)
    
#Plot_Observations(myobs[744],'DD',744)

#plt.show()

r=[]
bands='ugrizy'
corresp=dict(zip(bands+'a',[i for i in range(len(bands)+1)]))
corresp_inverted=dict(zip([i for i in range(len(bands)+1)],bands+'a'))

filtercolors = {'u':'c', 'g':'b', 'r':'g', 'i':'y', 'z':'r', 'y':'m'}
fieldcolors=dict(zip([309,310,1427,2412,2786],'bgyrm'))
#fieldcolors=dict(zip([fieldids[0],744,1427,2412,2786],'bgyrm'))
ra=[]
all_diff={}

for key, vals in myobs.items():
    print key,len(vals.seasons)
    all_diff[key]={}
    print 'nombre de saison',key,len(vals.seasons)
    #break
    nyears=len(vals.seasons)
    for season in range(nyears):
        
        myseason=vals.seasons[season]
        all_diff[key][season]={}
        full_season=myseason.copy()
        full_season.sort(order='mjd')
        idx = full_season['band'] != 'LSSTPG::u'
        full_season=full_season[idx]
        #print full_season[full_season['band'] == 'LSSTPG::i']
        #print full_season['band']
        
        ra+=Fill_Cadence(full_season,season)
        ro,alldiffs =Fill_Medians(myseason,season)
        r+=ro
        all_diff[key][season]=alldiffs
        
resu=np.rec.fromrecords(r,names=['fieldid','season','band','mean_cadence','rms_cadence','duration','ib','m5','m5_lower','m5_upper','obstime','median_cadence','upper_cadence','lower_cadence','nexp_sum','msky','msky_lower','msky_upper','airmass','airmass_lower','airmass_upper','seeing','seeing_lower','seeing_upper','ra','dec'])
Nobs=np.rec.fromrecords(ra,names=['season','band','T0','Nbef','Naft','Nmeas'])

pkl_file = open(fieldname+'.pkl','wb')
pkl.dump(resu, pkl_file)
pkl_file.close()
"""
     
r=[]
for season in range(10):
    fields_all=None
    for key in myobs.keys():
        if len(myobs[key].seasons) > season :
            myseason=myobs[key].seasons[season]
            if fields_all is None:
                fields_all=myseason
            else:
                fields_all=np.concatenate((fields_all,myseason))
    #print fields_all['Ra']
    for band in 'ugrizy':
        idx = fields_all['band']=='LSSTPG::'+band
        sel= fields_all[idx]
        m5_med,m5_lower,m5_upper=median(sel,'m5sigmadepth') 
        msky_med,msky_lower,msky_upper=median(sel,'sky')
        print band,m5_med,m5_lower,m5_upper,msky_med,msky_lower,msky_upper
        r.append((season,band,m5_med,m5_lower,m5_upper,msky_med,msky_lower,msky_upper))

resu=np.rec.fromrecords(r,names=['season','band','m5','m5_l','m5_u','msky','msky_l','msky_u'])


f = open('Cadence_tex/'+fieldname+'_m5msky.tex','w')
Print_begin_multi(f,fieldname)
        
for season in range(10):
    iid = resu['season']==season
    toprint=str(season+1)
    sel_seas=resu[iid]
    for what in ['m5','msky']:
        for band in 'ugrizy':
            ik=sel_seas['band']==band
            ssel=sel_seas[ik]
            med=np.asscalar(ssel[what])
            med_p=np.asscalar(ssel[what+'_u']-ssel[what])
            med_m=np.asscalar(ssel[what]-ssel[what+'_l'])
            #toprint+=' & '+str(round(med,2))+'$^{+'+str(round(med_p,2))+'}_{-'+str(round(med_m,2))+'}$'
            toprint+=' & '+str(round(med,1))
    print toprint
    toprint+=' \\\\'
    f.write(toprint+ '\n')
Print_end(f)

"""

#print 'alors?',len(resu)



"""
print resu
iseason=0
idx= resu['season'] == iseason
sela=resu[idx]
for fieldid in fieldids:
    idxb=sela['fieldid']==fieldid['fieldid']
    selb=sela[idxb]
    for band in 'grizy':
        idxc = selb['band']==band
        selc = selb[idxc]
        print fieldid, band, selc['mean_cadence'],selc['rms_cadence'],selc['duration']
"""

#Plot_per_Field(resu)
#Plot_per_Field(resu,what=('duration','obstime'),myleg=('Duration [day]','Observing Time [s]'))
#Plot_Diffs(all_diff,fieldname,fieldids,fieldcolors)
#Plot_Obs_per_Field(resu)
#Plot_Cadence_per_Year(resu)

"""
for field in fieldids:
    Plot_Cadence_vs_Year(resu,'WFD',field)
"""


#Plot_Nobs(Nobs)
#Plot_median(resu,fieldids,fieldcolors)
#Plot_median_m5(resu,fieldname, fieldids,filtercolors)
#Plot_median_m5(resu,fieldname, fieldids,filtercolors,what='airmass',leg='airmass')
#Plot_median_m5(resu,fieldname, fieldids,filtercolors,what='msky',leg='msky')
#Plot_median_m5(resu,fieldname, fieldids,filtercolors,what='seeing',leg='seeing')
plt.show()
