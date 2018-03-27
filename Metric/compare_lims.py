import numpy as np
import pylab as plt
from saunerie import psf
from matplotlib.legend import Legend

bands='griz'
lims={}
instr_snsim=psf.find('LSSTPG')

for simu in ['snsim','sncosmo']:
    for air in ['noairmass','airmass_12']:
        name=simu+'_'+air
        lims[name]=np.load('f5_cadence_lims_'+name+'.npy').item()

print(lims.keys(),type(lims['snsim_noairmass']))

r=[]
for key,vals in lims['snsim_noairmass'].items():
    val_snsim_airmass=lims['snsim_airmass_12'][key]
    val_sncosmo_airmass=lims['sncosmo_airmass_12'][key]
    val_sncosmo_noairmass=lims['sncosmo_noairmass'][key]
    for band in bands:
        print(vals['LSSTPG::'+band],val_snsim_airmass['LSSTPG::'+band],val_sncosmo_noairmass[band],val_sncosmo_airmass[band]) 
        b= ['LSSTPG::'+band]
        maga=np.asscalar(instr_snsim.flux_to_mag(vals['LSSTPG::'+band],b))
        magb=np.asscalar(instr_snsim.flux_to_mag(val_snsim_airmass['LSSTPG::'+band],b))
        magc=np.asscalar(instr_snsim.flux_to_mag(val_sncosmo_noairmass[band],b))
        magd=np.asscalar(instr_snsim.flux_to_mag(val_sncosmo_airmass[band],b))

        r.append((band,key,maga,magb,magc,magd))
        print(maga,magb,magc,magd)
res=np.rec.fromrecords(r,names=['band','z','snsim_noairmass','snsim_airmass_12','sncosmo_noairmass','sncosmo_airmass_12'])

figa, axa = plt.subplots(ncols=1, nrows=1)
filtercolors = {'u':'b', 'g':'g', 'r':'r', 'i':'c', 'z':'m', 'y':'m'}
tot_label=[]
tot_labelb=[]
for band in bands:
    #figa, axa = plt.subplots(ncols=1, nrows=1)
    idx=res['band']==band
    sel=res[idx]
    sel.sort(order='z')
    if band == 'r':
        tot_label.append(axa.errorbar(sel['z'],sel['snsim_noairmass']-sel['sncosmo_noairmass'],color=filtercolors[band],ls='-',label='snsim_noairmass-sncosmo_noairmass'))
        tot_label.append(axa.errorbar(sel['z'],sel['snsim_airmass_12']-sel['sncosmo_airmass_12'],color=filtercolors[band],ls='--',label='snsim_airmass_12-sncosmo_airmass_12'))
        tot_label.append(axa.errorbar(sel['z'],sel['snsim_airmass_12']-sel['snsim_noairmass'],color=filtercolors[band],ls='-.',label='snsim_airmass_12-snsim_noairmass'))
        tot_label.append(axa.errorbar(sel['z'],sel['sncosmo_airmass_12']-sel['sncosmo_noairmass'],color=filtercolors[band],ls=':',label='sncosmo_airmass_12-sncosmo_noairmass'))
        tot_labelb.append(axa.errorbar(sel['z'],sel['snsim_noairmass']-sel['sncosmo_noairmass'],color=filtercolors[band],ls='-',label=band+' band'))
    else:
        tot_labelb.append(axa.errorbar(sel['z'],sel['snsim_noairmass']-sel['sncosmo_noairmass'],color=filtercolors[band],ls='-',label=band+' band'))
        axa.errorbar(sel['z'],sel['snsim_airmass_12']-sel['sncosmo_airmass_12'],color=filtercolors[band],ls='--')
        axa.errorbar(sel['z'],sel['snsim_airmass_12']-sel['snsim_noairmass'],color=filtercolors[band],ls='-.')
        axa.errorbar(sel['z'],sel['sncosmo_airmass_12']-sel['sncosmo_noairmass'],color=filtercolors[band],ls=':')  
        
labs = [l.get_label() for l in tot_label]
axa.legend(tot_label, labs, ncol=2,loc='lower left',prop={'size':12},frameon=False)

labsb = [l.get_label() for l in tot_labelb]

leg=Legend(axa,tot_labelb, labsb,ncol=4,loc='upper right',prop={'size':12},frameon=False)
axa.add_artist(leg)

axa.set_xlabel('z')
axa.set_ylabel('$\Delta m_{5}$ [mag]')
axa.set_xlim([0.095,0.505])

plt.show()
