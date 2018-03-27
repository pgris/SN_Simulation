import cPickle as pkl
import numpy as np
import glob
import matplotlib.pyplot as plt
"""
from matplotlib import rc

# activate latex text rendering
rc('text', usetex=True)
"""
def Read_File(fieldname,fieldids,name_added='',thedir=''):
   tot_nsn=None
   for fieldid in fieldids:
       name='N_SN_'+fieldname+'_'+str(fieldid)
       if name_added != '':
           name='N_SN_'+fieldname+'_'+str(fieldid)+'_'+name_added

       files=glob.glob(thedir+name+'_Season_*')
       #print('Loading',thedir,files)
       for fi in files:
          pkl_file = open(fi,'rb')
          if tot_nsn is None:
             tot_nsn=pkl.load(pkl_file)
          else:
             tot_nsn=np.vstack((tot_nsn,pkl.load(pkl_file))) 

   return tot_nsn

def Get_Numbers(sel):
   
   nsn=[]
   err=[]
   season=[]
        #print 'there',sel
   for i in range(10):
      iid = sel['season']==i
      if sel['n_sn_detected'][iid] >0:
         season.append(i+1)
         nsn.append(np.asscalar(sel['n_sn_detected'][iid]))
         err.append(np.asscalar(sel['err_detected'][iid]))

   return season,nsn,err

def Plot_NSN_new(axa,tot_nsn,fieldid,ls='-',tot_label=[],labname='',draw_legend=False,res=[]):

   #tot_label=[]
   fontsize=12
    #minval=220.
   minval=140.
    #interv=15
   interv=5
   #res=[]
      
   season,nsn,err=Get_Numbers(tot_nsn)
   nsn_tot=np.sum(nsn)
   err_tot=np.sqrt(np.sum([val*val for val in err]))
   ll='Field '+str(fieldid)
   #axa.errorbar(season,nsn,yerr=err,label=ll,color=colors[fieldid],ls=ls)
   tot_label.append(axa.errorbar(season,nsn,yerr=err,label=labname,color='k',ls=ls))
   for i in range(len(season)):
      res.append((fieldid,season[i],labname,nsn[i],err[i]))
  

def Plot_NSN(axa,tot_nsn,fieldids,ls='-',draw_title=True,draw_legend=True,draw_numbers=True,add_info=None,colors=None):

    tot_label=[]
    fontsize=12
    #minval=220.
    minval=140.
    #interv=15
    interv=5
    res=[]
    for fieldid in fieldids:
        
        season,nsn,err=Get_Numbers(fieldid,tot_nsn)
        nsn_tot=np.sum(nsn)
        err_tot=np.sqrt(np.sum(err*err))
        ll='Field '+str(fieldid)
        if draw_legend:
            tot_label.append(axa.errorbar(season,nsn,yerr=err,label=ll,color=colors[fieldid],ls=ls))
        else:
            #print 'hello',nsn
            axa.errorbar(season,nsn,yerr=err,label=ll,color=colors[fieldid],ls=ls)
        
    print 'boo',tot_nsn.dtype,tot_nsn['season']
    for i in range(10):
        idf = tot_nsn['season']==i
        selb= tot_nsn[idf]
        print i, np.sum(selb['n_sn_detected']),np.sqrt(np.sum(np.power(selb['err_detected'],2.)))
        res.append((i, np.sum(selb['n_sn_detected']),np.sqrt(np.sum(np.power(selb['err_detected'],2.)))))

    nsn_per=np.rec.fromrecords(res,names=['period','nsn','err_nsn'])   
    
    
    axa.set_xlabel('Year',{'fontsize': fontsize})
    axa.set_ylabel('Number of Type Ia Supernovae',{'fontsize': fontsize})
    axa.set_xlim([0.5,10.5])
    
    if draw_legend:
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
                axa.text(5,minval-(i%5)*interv,thetext)
            else:
                axa.text(8,minval-(i%5)*interv,thetext) 

    if draw_title:
        ntot_nsn=np.sum(nsn_per['nsn'])
        errtot_nsn=np.sqrt(np.sum(np.power(nsn_per['err_nsn'],2.)))
        myt='Perrett rate - 10 years - $\mathrm{N_{SN\/Ia}}$ = '+str(int(ntot_nsn))+' $\pm$ '+str(int(errtot_nsn))
        if add_info is not None:
            ntot_nsnb=np.sum(add_info['nsn'])
            errtot_nsnb=np.sqrt(np.sum(np.power(add_info['err_nsn'],2.)))
            myt+=' / '+str(int(ntot_nsnb))+' $\pm$ '+str(int(errtot_nsnb))
        figa.suptitle(myt) 

    return nsn_per


def Plot_Field(fieldid,tot_nsn,rtot=[]):

   figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(12,9))
   figa.suptitle('Perrett rate - Field '+str(fieldid))
   tot_label=[]
   io=-1
   lsty=['-','--','-.',':',(0, (5, 10)),(0, (5, 10))]
   for key,vals in tot_nsn.items():
      io+=1
      idb = vals['fieldid']==fieldid
      Plot_NSN_new(axa,vals[idb],fieldid,ls=lsty[io],tot_label=tot_label,labname=key,res=rtot)
      
   axa.set_xlabel('Year',{'fontsize': fontsize})
   axa.set_ylabel('Number of Type Ia Supernovae',{'fontsize': fontsize})
   axa.set_xlim([0.5,10.5])
   labs = [l.get_label() for l in tot_label]
   axa.legend(tot_label, labs, ncol=2,loc='best',prop={'size':fontsize},frameon=False)
   

fieldname='DD'
#fieldids=[744,1427,2412,2786]
fieldids=[100,101,102,103]
#fieldids=[1427]
simu_name='sncosmo'
#selection=dict(zip(['N_SN_sncosmo_Science_Book_no_colorcut','N_SN_sncosmo_Science_Book_with_colorcut','N_SN_sncosmo_Sela_no_color_cut','N_SN_sncosmo_Sela_with_color_cut'],['Science_Book','Science_Book_with_colorcut','Sela','Sela_with_color_cut']))
selection_names=['Sela','Sela_cc','Selb_10_2','Selb_10_2_cc']
tot_nsn={}
for selname in selection_names:
   tot_nsn[selname]=Read_File(fieldname,fieldids,thedir='N_SN_'+simu_name+'_'+selname+'/')

#tot_nsn=Read_File(fieldname,fieldids,thedir='Cad_duration_180day/')
#tot_nsn_sncosmo=Read_File(fieldname,fieldids,thedir='NSN_sncosmo_nocolorcut/')
#tot_nsn_cut=Read_File('colorcut')

#print tot_nsn_sncosmo

colors_minion = dict(zip(fieldids,['b','k','g','r','m']))

colors = dict(zip(fieldids,['b','k','g','r','m']))
fontsize=12

rtot=[]
for fieldid in fieldids:
   Plot_Field(fieldid,tot_nsn,rtot=rtot)
  
print rtot
nsn_stat=np.rec.fromrecords(rtot,names=['fieldname','season','selection','nsn','err_nsn'])


r=[]
for sel in np.unique(nsn_stat['selection']):
   ida = nsn_stat['selection']==sel
   sela=nsn_stat[ida]
   for season in range(1,11):
      idx=sela['season']==season
      selb=sela[idx]
      print(sel,season,np.sum(selb['nsn']),np.sqrt(np.sum(selb['err_nsn']*selb['err_nsn'])))
      r.append((sel,season,np.sum(selb['nsn']),np.sqrt(np.sum(selb['err_nsn']*selb['err_nsn']))))

nsn_all=np.rec.fromrecords(r,names=['selection','season','nsn','err_nsn'])


#now write this in latex
#selections=['Sela', 'Sela_with_color_cut','Science_Book' ,'Science_Book_with_colorcut']

"""
idx={}
for sel in selection:
   idx[sel] = nsn_all['selection'] == sel
"""
print('Season & Sela & Sela_cc &  Selb_10_2 &  Selb_10_2_cc')
totsum={}
for sel in selection_names:
   totsum[sel]={}
   totsum[sel]['nsn']=0
   totsum[sel]['err_nsn']=0
for season in range(1,11):
   str_tot=str(season)+' & '
   for j,sel in enumerate(selection_names):
      #print(nsn_all['selection'],nsn_all['season'],sel,season)
      idx= (nsn_all['selection'] == sel)&(nsn_all['season']==season)
      #print(nsn_all[idx]['nsn'],nsn_all[idx]['err_nsn'])
      str_tot+=str(int(nsn_all[idx]['nsn']))+' $\pm$ '+str(int(nsn_all[idx]['err_nsn']))
      totsum[sel]['nsn']+=nsn_all[idx]['nsn']
      totsum[sel]['err_nsn']+=nsn_all[idx]['err_nsn']*nsn_all[idx]['err_nsn']
      if j != len(selection_names)-1:
         str_tot+=' & '
      else:
         str_tot+=' \\\\'
   print(str_tot)
str_fi = 'all & '
for j,sel in enumerate(selection_names):
   str_fi+=str(int(totsum[sel]['nsn']))+' $\pm$ '+str(int(np.sqrt(totsum[sel]['err_nsn'])))
   if j != len(selection_names)-1:
      str_fi+=' & '
   else:
      str_fi+=' \\\\'
print(str_fi)



#figa.savefig('Plots_NSN/Summary_'+fieldname+'.png')

plt.show()
