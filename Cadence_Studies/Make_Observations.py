import numpy as np


bands='ugrizy'
m5=dict(zip(bands,[23.61,24.83,24.35,23.88,23.30,22.43]))
msky = dict(zip(bands,[22.95,22.24,21.20,20.47,19.60,18.63]))
kAtm = dict(zip(bands,[0.50,0.21,0.13,0.10,0.07,0.18]))
seeing=dict(zip(bands,[0.92,0.87,0.83,0.80,0.78,0.76]))
Tvisit=30. #in sec
Nvisits=dict(zip(bands,[10,10,20,20,20,26,20]))

cadence = 4. #in days

mjd_min=100.
mjd_max=mjd_min+140.


m5_coadd={}

for key, val in m5.items():
    m5_coadd[key]=val+1.25*np.log10(float(Nvisits[key])*Tvisit/30.)

fieldid=100
fi=open('Observations/Obs_'+str(fieldid)+'_'+str(int(cadence))+'.txt',"w")

entete=['band','mjd','exptime','rawSeeing','seeing','moon_frac','sky','kAtm','airmass','m5sigmadepth','Nexp','Ra','Dec']
towrite=''
for val in entete:
    towrite+='# '+val+' :\n'
towrite+='# end'
fi.write(towrite+'\n')
moon_frac=0.
airmass=1.2
Ra=0.0
Dec=0.

shifta=10./(60.*24.)
shiftb=cadence/2.
shift_days=dict(zip('grizy',[0.,shifta,2.*shifta,3.*shifta,4.*shifta]))
#shift_days=dict(zip('grizy',[0.,shifta,shiftb,shiftb+shifta,shiftb+2.*shifta]))

print(entete)
for mjd in np.arange(mjd_min,mjd_max+cadence,cadence):
    for i,band in enumerate('grizy'):
        mjd_val=mjd+shift_days[band]
        if mjd_val <= mjd_max+cadence:
            obs=''
            obs+='LSST::'+band+' '+str(mjd_val)+' '+str(float(Nvisits[band])*Tvisit)+' -1.00 '+str(seeing[band])+' '+str(moon_frac)+' '+str(msky[band])
            obs+=' '+str(kAtm[band])+' '+str(airmass)+' '+str(m5_coadd[band])+' '+str(Nvisits[band])+' '+str(Ra)+' '+str(Dec)
            print(obs)
            fi.write(obs+'\n')

fi.close()
