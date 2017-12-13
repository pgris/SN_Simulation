from SN_Rate import *
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np

def Plot_Rate(zz,ratea,erra,zzb=None,rateb=None,errb=None):

    figbb, axbb = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
    axbb.errorbar(zz,1.e4*ratea,yerr=1.e4*erra,marker='.', mfc='black', mec='black', ms=8, linestyle='-',color='blue',label=rate_name+' rate')
    if rateb is not None:
        axbb.errorbar(zzb,1.e4*rateb,yerr=1.e4*errb,marker='.', mfc='black', mec='black', ms=8, linestyle='-',color='red',label=rate_nameb+' rate')

    axbb.set_yscale('log')
    axbb.set_xlabel('z')
    axbb.set_ylabel('10$^{-4}$ SNe Ia yr$^{-1}$ Mpc$^{-3}$ h$^{3}$$_{70}$')
    axbb.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    axbb.set_ylim(top=1.1)
    plt.legend(loc='best')
    
def Plot_NSN_per_z_bin(survey_area,zz,nsn,err_nsn,zzb=None,nsnb=None,err_nsnb=None,cumul=False):

    tot_sn=np.sum(nsn)
    err_tot_sn=np.power(np.sum(np.power(err_nsn,2.)),0.5)

    duration=0

    figc, axc = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
    ll='N$_{SN Ia}$ = '+str(int(tot_sn))+'$\pm$'+str(int(err_tot_sn))+' - '+str(int(duration))+' days'
    figc.suptitle('Survey area: '+str(survey_area)+' deg$^{2}$ - Duration: '+str(duration)+' year')
    
    if not cumul:
        axc.errorbar(zz,nsn,yerr=err_nsn,marker='.', mfc='black', mec='black', ms=8, linestyle='-',color='blue',label=rate_name+' rate - '+ll)
    else:
        axc.errorbar(zz,np.cumsum(nsn),yerr=err_nsn,marker='.', mfc='black', mec='black', ms=8, linestyle='-',color='blue',label=rate_name+' rate - '+ll)
    
    if zzb is not None:
        tot_snb=np.sum(nsnb)
        err_tot_snb=np.power(np.sum(np.power(err_nsnb,2.)),0.5)
        ll='N$_{SN Ia}$ = '+str(int(tot_snb))+'$\pm$'+str(int(err_tot_snb))+' - '+str(int(durationb))+' days'
        if not cumul:
            axc.errorbar(zzb,nsnb,yerr=err_nsnb,marker='.', mfc='black', mec='black', ms=8, linestyle='-',color='red',label=rate_nameb+' rate - '+ll)
        else:
            axc.errorbar(zzb,np.cumsum(nsnb),yerr=err_nsnb,marker='.', mfc='black', mec='black', ms=8, linestyle='-',color='red',label=rate_nameb+' rate - '+ll)  
    
    axc.set_xlabel('z')
    axc.set_ylabel('N(SNe Ia) / 0.01 z-bin')
    axc.set_xlim([0.,1.])
    axc.legend(loc='best',prop={'size':12})

def Plot_NSN_vs_z(axc,survey_area,data,cumul=False):

    for val in data:
         tot_sn=np.sum(val['nsn'])
         err_tot_sn=np.sqrt(np.sum(np.cumsum(np.power(val['err_nsn'],2.)),0.5))
         yerr=[]
         yerr.append(val['err_nsn'][0]**2)
         for io in range(1,len(val['err_nsn'])):
             suma=np.sum([val['err_nsn'][i]*val['err_nsn'][i] for i in range(1,io+1)])
             
             sumb=0.
             for i in range(1,io+1):
                 for j in range(1,io+1):
                     if j != i:
                         sumb+=val['err_nsn'][i]*val['err_nsn'][j]
                
             yerr.append(suma+sumb)
        
         err_tot_sn=np.sqrt(yerr[-1])

         ll='N$_{SN Ia}$ = '+str(int(tot_sn))+'$\pm$'+str(int(err_tot_sn))+' - '+str(int(val['duration']))+' days'
         if not cumul:
             
             #axc.errorbar(val['zz'],val['nsn'],marker='.', mfc=val['color'], mec=val['color'], linestyle='-',color=val['color'],label=val['rate_name']+' rate - '+ll)
             axc.errorbar(val['zz'],val['nsn'],marker=None, linestyle='-',color=val['color'],label=val['rate_name']+' rate - '+ll)
             axc.errorbar(val['zz'],val['nsn']+val['err_nsn'],marker=None, linestyle='--',color=val['color'])

             axc.errorbar(val['zz'],val['nsn']-val['err_nsn'],marker=None, linestyle='--',color=val['color'])

             #axc.errorbar(val['zz'],val['nsn'],marker='.', mfc='k', mec='k', ms=1, linestyle='-',color='k',label=val['rate_name']+' rate - '+ll)
             #axc.fill_between(val['zz'], val['nsn']-val['err_nsn'], val['nsn']+val['err_nsn'],color=val['color'])
         else:
             #axc.errorbar(val['zz'],np.cumsum(val['nsn']),yerr=np.sqrt(yerr),marker='.', mfc=val['color'], mec=val['color'], ms=8, linestyle='-',color=val['color'],label=val['rate_name']+' rate - '+ll)
             axc.errorbar(val['zz'],np.cumsum(val['nsn']),marker=None, linestyle='-',color=val['color'],label=val['rate_name']+' rate - '+ll)
             axc.errorbar(val['zz'],np.cumsum(val['nsn'])+np.sqrt(yerr),marker=None, linestyle='--',color=val['color'])
             axc.errorbar(val['zz'],np.cumsum(val['nsn'])-np.sqrt(yerr),marker=None,linestyle='--',color=val['color'])
             
             
             

    axc.grid()
    axc.set_xlabel('z')
    if cumul is False:
        axc.set_ylabel('N(SNe Ia) / 0.01 z-bin')
    else:
        axc.set_ylabel('N(SNe Ia) (z < )')
    axc.set_xlim([0.,1.])
    axc.legend(loc='best',prop={'size':12})

def Plot_NSN_Cumul_vs_z(axc,survey_area,data):
   
    #figc, axc = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
    #figc.suptitle('Survey area: '+str(survey_area)+' deg$^{2}$')

    for val in data:
        print 'aaa',len(val['err_nsn'])
        tot_sn=val['nsn'][-1]
        err_tot_sn=val['err_nsn'][-1]

        ll='N$_{SN Ia}$ = '+str(int(tot_sn))+'$\pm$'+str(int(err_tot_sn))+' - '+str(int(val['duration']))+' days'
        axc.errorbar(val['zz'],val['nsn'],yerr=val['err_nsn'],marker='.', mfc=val['color'], mec=val['color'], ms=8, linestyle='-',color=val['color'],label=val['rate_name']+' rate - '+ll)
        
        
    axc.set_xlabel('z')
    
    axc.set_ylabel('N(SNe Ia) (z < )')
    axc.set_xlim([0.,1.])
    axc.legend(loc='best',prop={'size':12})

nfields=2293
nfields=1.
survey_area= 9.6*nfields
duration = 140.
zmin=0.01
zmax=1.
zstep=0.01

rate_name='Perrett'
rate=SN_Rate(rate=rate_name,survey_area=survey_area,duration=duration/365.)

zz,rate_one,err_rate,nsn,err_nsn=rate(zmin,zmax,zstep)

print 'rate',rate_one,err_rate

rate_nameb='Perrett'
durationb=200.
rateb=SN_Rate(rate=rate_nameb,survey_area=survey_area,duration=durationb/365.)
zzb,rate_two,err_rateb,nsnb,err_nsnb=rateb(zmin,zmax,zstep)



#print zz,nsn,nsnb

#Plot_Rate(zz,rate_one,err_rate)


case=[]

case.append(dict(zip(['zz','nsn','err_nsn','duration','color','rate_name'],[zz,nsn,err_nsn,rate.duration*365.,'b',rate_name])))
case.append(dict(zip(['zz','nsn','err_nsn','duration','color','rate_name'],[zzb,nsnb,err_nsnb,rateb.duration*365.,'r',rate_nameb])))

print 'error',nsn,err_nsn,nsn,err_nsn/nsn
caseb=[]
zz_a=[]
nsn_a=[]
err_nsn_a=[]
nsn_b=[]
err_nsn_b=[]

#z=np.linspace(0.01,1.,99)
z=np.linspace(0.,1.,11)
nsn,err_nsn=rate.N_SN(z)
nsnb,err_nsnb=rateb.N_SN(z)

caseb.append(dict(zip(['zz','nsn','err_nsn','duration','color','rate_name'],[z,nsn,err_nsn,rate.duration*365.,'b',rate_name])))
caseb.append(dict(zip(['zz','nsn','err_nsn','duration','color','rate_name'],[z,nsnb,err_nsnb,rateb.duration*365.,'r',rate_nameb])))   

#print 'errors',len(nsn),len(err_nsn),z,nsn,err_nsn
figc, axc = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
figc.suptitle('Survey area: '+str(survey_area)+' deg$^{2}$')
Plot_NSN_vs_z(axc,survey_area,case,cumul=False)

#Plot_NSN_Cumul_vs_z(axc,survey_area,caseb)

plt.show()


"""
rate_namec='Perrett'
r=[]
f = open('nsn.tex','w')
f.write('\\begin{table}[ht]'+'\n')
#f.write('\\begin{center}'+'\n')
f.write('\\caption{Number of type Ia supernovae with ten years of observation and a field-of-view of 9.6 deg$^{2}}\label{tab:nsn}'+'\n')
f.write('\\begin{tabular}{llllllllll}'+'\n')
f.write('\\hline'+'\n')

chap='duration'
for z in np.linspace(0.1,1.0,10):
    chap+=' & z $\leq$ '+str(z)
chap+=' \\\\'
f.write(chap+'\n')
for duration in np.arange(130,210,10):
    thestr=str(duration)
    duration*=10 # for ten years of observations
    
    for z in np.linspace(0.1,1.1,11):
        ratec=SN_Rate(rate=rate_namec,survey_area=survey_area,duration=duration/365.)
        nsn,err_nsn=ratec.N_SN(z)
        print duration,z,nsn,err_nsn
        thestr+=' & '+str(int(nsn))+'$\pm$'+str(int(err_nsn))
        r.append((duration,z,nsn,err_nsn))
    thestr+=' \\\\'
    f.write(thestr+'\n')
    print thestr
f.write('\\hline'+'\n')
f.write('\\end{tabular}'+'\n')
#f.write('\\end{center}'+'\n')
f.write('\\end{table}'+'\n')

f.close()

res=np.rec.fromrecords(r,names=['duration','z','sns','err_nsn'])
"""
    
