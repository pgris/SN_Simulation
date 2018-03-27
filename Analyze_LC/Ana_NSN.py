import cPickle as pkl
import numpy as np
import glob
import matplotlib.pyplot as plt
"""
from matplotlib import rc

# activate latex text rendering
rc('text', usetex=True)
"""
def Read_File(name_added=''):
   tot_nsn=None
   for fieldid in fieldids:
       name='N_SN_'+fieldname+'_'+str(fieldid)
       if name_added != '':
           name='N_SN_'+fieldname+'_'+str(fieldid)+'_'+name_added

       pkl_file = open(name+'.pkl','rb')
       if tot_nsn is None:
           tot_nsn=pkl.load(pkl_file)
       else:
           tot_nsn=np.vstack((tot_nsn,pkl.load(pkl_file))) 


   return tot_nsn

def Plot_NSN(axa,tot_nsn,ls='-',draw_legend=True,draw_numbers=True,add_info=None):

    tot_label=[]
    fontsize=12
    minval=220.
    interv=15
    res=[]
    for fieldid in fieldids:
        idx = tot_nsn['fieldid']==fieldid
        sel=tot_nsn[idx]
        nsn=[]
        err=[]
        season=[]
        for i in range(10):
            iid = sel['season']==i
            if sel['n_sn_detected'][iid] >0:
               season.append(i+1)
               nsn.append(sel['n_sn_detected'][iid])
               err.append(sel['err_detected'][iid])
        
        ll='Field '+str(fieldid)
        if draw_legend:
            tot_label.append(axa.errorbar(season,nsn,yerr=err,label=ll,color=colors[fieldid]))
        else:
            axa.errorbar(season,nsn,yerr=err,label=ll,color=colors[fieldid],ls=ls)
        
    print 'boo',tot_nsn.dtype,tot_nsn['season']
    for i in range(10):
        idf = tot_nsn['season']==i
        selb= tot_nsn[idf]
        print i, np.sum(selb['n_sn_detected']),np.sqrt(np.sum(np.power(selb['err_detected'],2.)))
        res.append((i, np.sum(selb['n_sn_detected']),np.sqrt(np.sum(np.power(selb['err_detected'],2.)))))

    nsn_per=np.rec.fromrecords(res,names=['period','nsn','err_nsn'])   
    
    if draw_legend:
        axa.set_xlabel('Year',{'fontsize': fontsize})
        axa.set_ylabel('Number of Type Ia Supernovae',{'fontsize': fontsize})
        axa.set_xlim([0.5,10.5])

        labs = [l.get_label() for l in tot_label]
        axa.legend(tot_label, labs, ncol=5,loc='best',prop={'size':fontsize},frameon=False)
        
    if draw_numbers:
        for i in range(10):
            idx = nsn_per['period']==i
            selp=nsn_per[idx]
            n_sn=int(selp['nsn'])
            err_sn=int(selp['err_nsn'])
            print 'hello',i,n_sn,err_sn
            thetext='Year '+str(i+1)+' : $\mathrm{N_{SN\/ Ia}}$ = '+str(n_sn)+'$\pm$'+str(err_sn)
            if add_info is not None:
                idxb = add_info['period']==i
                selpb=add_info[idxb]
                n_snb=int(selpb['nsn'])
                err_snb=int(selpb['err_nsn'])
                thetext+=' / '+str(n_snb)+'$\pm$'+str(err_snb)

            if i < 5:
                axa.text(1,minval-(i%5)*interv,thetext)
            else:
                axa.text(5,minval-(i%5)*interv,thetext) 

        ntot_nsn=np.sum(nsn_per['nsn'])
        errtot_nsn=np.sqrt(np.sum(np.power(selp['err_nsn'],2.)))
        myt='Perrett rate - 10 years - $\mathrm{N_{SN\/Ia}}$ = '+str(int(ntot_nsn))+' $\pm$ '+str(int(errtot_nsn))
        if add_info is not None:
            ntot_nsnb=np.sum(add_info['nsn'])
            errtot_nsnb=np.sqrt(np.sum(np.power(add_info['err_nsn'],2.)))
            myt+=' / '+str(int(ntot_nsnb))+' $\pm$ '+str(int(errtot_nsnb))
        figa.suptitle(myt) 

    return nsn_per

fieldname='DD'
fieldids=[290,744,1427,2412,2786]
#fieldids=[100,101,102,103]

tot_nsn=Read_File()
#tot_nsn_cut=Read_File('colorcut')

print tot_nsn

figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(12,9))
colors = dict(zip(fieldids,['b','k','g','r','m']))

nsn_per=Plot_NSN(axa,tot_nsn,draw_legend=True,draw_numbers=True)
#Plot_NSN(axa,tot_nsn_cut,draw_legend=False,ls='--',add_info=nsn_per)


figa.savefig('Plots_NSN/Summary_'+fieldname+'.png')
plt.show()
