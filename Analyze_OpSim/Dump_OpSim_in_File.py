import pickle as pkl
#import cPickle as pkl
import numpy as np
from optparse import OptionParser
import os
from Parameters import parameters
import math

expMJD='observationStartMJD'
#skyb='filtSkyBrightness'
skyb='skyBrightness'
#visitTime='visitExpTime'
visitTime='visitExposureTime'
rawSeeing='seeingFwhmGeom'
#rawSeeing='rawSeeing'


def Get_mean_finSeeing(observations):
    
    res=[]
    params=parameters()
    for obs in observations:
        filtre=obs['filter'][0]
        seeing=obs[rawSeeing]
        airmass=obs['airmass']
        Filter_Wavelength_Correction = np.power(500.0 / params.filterWave[filtre], 0.3)
        Airmass_Correction = math.pow(obs['airmass'],0.6)
        FWHM_Sys = params.FWHM_Sys_Zenith * Airmass_Correction
        FWHM_Atm = seeing * Filter_Wavelength_Correction * Airmass_Correction
        finSeeing = params.scaleToNeff * math.sqrt(np.power(FWHM_Sys,2) + params.atmNeffFactor * np.power(FWHM_Atm,2))
        res.append(finSeeing)
    return np.mean(res)

def Get_coadd(fieldname,filt,newsky=False,diff_time=40.):

    if len(filt) == 0:
        return None

    #print(type(filt))
    
    diff_time=diff_time # in sec

    """
    if fieldname=='WFD':
        diff_time=0.
    """
    filt.sort(order=expMJD)
    """
    for val in filt:
        print val
    print np.mean(filt['expMJD'][:10])
    """
    #print filt['expMJD']
    
    diff=[io-jo for jo,io in zip(filt[expMJD][:-1], filt[expMJD][1:])]
    #print('hello',diff)
    sel=[i+1 for i in range(len(diff)) if 24.*3600.*diff[i]> diff_time]
    sel=[0]+sel+[len(filt)]
    r=[]

    #print(sel)

    #var_tot=['filter','fieldRA','fieldDec','ditheredRA','ditheredDec']
    var_tot=['filter','fieldRA','fieldDec']
    vars_mean=[expMJD,rawSeeing,'moonPhase',skyb,'airmass','fiveSigmaDepth']

    
    for i in range(len(sel)-1):
        ia=sel[i]
        ib=sel[i+1]
       
        restot=dict(zip([var for var in var_tot],[filt[var][ia] for var in var_tot]))

        theslice=filt[ia:ib]

        res=dict(zip([var for var in vars_mean],[np.mean(theslice[var]) for var in vars_mean]))
        res['visitExpTime']=np.sum(filt[visitTime][ia:ib])
        #res['finSeeing']=Get_mean_finSeeing(theslice)
        res['finSeeing']=np.mean(theslice['seeingFwhmEff'])
        res['kAtm']=params.kAtm[filt['filter'][ia]]
        #if fieldname == 'DD':
        res['fiveSigmaDepth']+=1.25*np.log10(res['visitExpTime']/30.)               
        res['Nvisits']=len(theslice)
        restot.update(res)
       
        r.append(tuple([restot[key] for key in restot.keys()]))
       
    resu=np.rec.fromrecords(r,names=[key for key in restot.keys()])   

    #print resu['visitExpTime'],len(resu)
    #print('Done')
    return resu


def Get_coadd_old(filt,newsky=False):
   
    dict_for_coadd={}
    filtc=filt.copy()
    filtc.reshape((filtc.size,1))

    #print 'alors',len(filt)
    if newsky:
        transmission=Throughputs()
        sm = sb.SkyModel(observatory='LSST', mags=False, preciseAltAz=True)
    if len(filtc) > 0:
        inum=0
        if newsky:
        #first : recompute the sky
            for i,val in enumerate(filtc):
                filtre=val['filter']
                ra_rad=val['fieldRA']
                dec_rad=val['fieldDec']

                mjd=val[expMJD]
                seeing=val[rawSeeing]
                visitExpTime=val['visitExpTime']
                airmass=val['airmass']

                #print 'airmass',airmass
                if airmass < 2.5:
                    transmission.Load_Atmosphere(airmass)
                    sm.setRaDecMjd(lon=[ra_rad], lat=[dec_rad], filterNames=[filtre],
                                   mjd=mjd, degrees=False, azAlt=False)
                    mag=sm.returnMags(bandpasses=transmission.lsst_atmos)[filtre][0]
                    filtc['filtSkyBrightness'][i]=mag
                    filtc['fiveSigmaDepth'][i]=m5_OpSim(filtre,seeing,mag,visitExpTime,airmass)
    
    
        dict_for_coadd[inum]=np.zeros((0,1),filtc.dtype)
        #print 'timediff',24.*60.*60.*(filtc['time']-filtc['time'][0])
                                
        iloop=0
        #print 'blablabla',dict_for_coadd[inum]
       
        dict_for_coadd[inum]=np.vstack([dict_for_coadd[inum],filtc[iloop]])
                                
        if len(filtc) > 1:
            while iloop < len(filtc)-1:   
                diff_time_sec=24.*60.*60.*(filtc[expMJD][iloop+1]-filtc[expMJD][iloop])
                #print 'alors ???',diff_time_sec,inum
                if diff_time_sec > 40.:
                    inum+=1
                    dict_for_coadd[inum]=np.zeros((0,1),filtc.dtype)
                
                dict_for_coadd[inum]=np.vstack([dict_for_coadd[inum],filtc[iloop]])
                    
                iloop+=1
        #print 'thedict',dict_for_coadd

    return dict_for_coadd

def Do_It(fieldname,band,toprocess,newsky=False,diff_time=40.,out_q=None):
        
    res=[]
    coadd=Get_coadd(fieldname,toprocess[toprocess['filter']==band],newsky=newsky,diff_time=diff_time)
    #print 'ee',band,len(coadd)
    if coadd is not None:
        for val in coadd:
            filtre=val['filter']
            ra_rad=val['fieldRA']
            dec_rad=val['fieldDec']
            toprint = 'LSSTPG::'+filtre+' '
            toprint+=str(format(val[expMJD],'.7f'))+' '
            toprint+=str(int(val['visitExpTime']))+' '
            #toprint+=str(format(np.mean(Get_mean_finSeeing(val)),'.7f'))+' '
            toprint+=str(format(val[rawSeeing],'.7f'))+' '
            #toprint+=str(format(np.mean(val['FWHMeff']),'.7f'))+' '
            toprint+=str(format(val['finSeeing'],'.7f'))+' '
            toprint+=str(format(val['moonPhase'],'.7f'))+' '
	    #print sm.getComputedVals()
        
            toprint+=str(format(val['skyBrightness'],'.7f'))+' '
            toprint+=str(val['kAtm'])+' '
            toprint+=str(format(val['airmass'],'.7f'))+' '
            #toprint+=str(format(np.mean(val['FWHMgeom']),'.7f'))+' '
            #toprint+=str(format(np.mean(val['FWHMeff']),'.7f'))+' '
            toprint+=str(format(val['fiveSigmaDepth'],'.7f'))+' '
            toprint+=str(val['Nvisits'])+' '
            toprint+=str(format(val['fieldRA'],'.7f'))+' '
            toprint+=str(format(val['fieldDec'],'.7f'))
            
            #toprint+=str(format(val['ditheredRA'],'.7f'))+' '
            #toprint+=str(format(val['ditheredDec'],'.7f'))
            #m5_coadd=Get_fiveSigmaDepth_coadd(len(val),np.mean(val['fiveSigmaDepth']))
            #m5_coadd=Get_fiveSigmaDepth_coadd(1.,np.mean(val['fiveSigmaDepth']))
            #print 'hello',np.mean(val['fiveSigmaDepth']),m5_coadd
            #outputfile.write(toprint+'\n')
            res.append(toprint+'\n')
            #print(toprint)

    if out_q is not None:
        out_q.put({band : res})
    else:
        return res

 

parser = OptionParser()
parser.add_option("-f", "--fieldname", type="string", default="DD", help="filter [%default]")
parser.add_option("-F", "--fieldid", type="int", default=290, help="filter [%default]")
parser.add_option("--fieldRA", type="float", default=0., help="filter [%default]")
parser.add_option("--fieldDec", type="float", default=0., help="filter [%default]")
parser.add_option("--season", type="int", default=1, help="filter [%default]")
parser.add_option("--simu", type="string", default="feature_baseline_10yrs", help="filter [%default]")
parser.add_option("--coadd", type="string", default="no", help="filter [%default]")
parser.add_option("--time_coadd", type="float", default=24.*3600., help="filter [%default]")

opts, args = parser.parse_args()

bands='ugrizy'
#bands='g'
prefix='Observations'
prefixb='Summary'
fieldname=opts.fieldname
fieldid=opts.fieldid
fieldRA=opts.fieldRA
fieldDec=opts.fieldDec
season=opts.season
coadd=''
if opts.coadd == 'yes':
    coadd='_coadd'
#name=prefix+'_'+opts.fieldname+'_'+str(opts.fieldid)+'.pkl'
#name=prefix+'_'+opts.fieldname+'_'+str(opts.fieldid)+'.npy'
name=prefix+'_'+str(opts.fieldid)+'.npy'
name_summary=prefixb+'_'+str(opts.fieldid)+'.npy'
#name=prefix+'_'+opts.fieldname+'_'+str(opts.fieldRA)+'_'+str(opts.fieldDec)+'.pkl'

suffix=opts.simu

thedir='/sps/lsst/data/dev/pgris/OpSimFiles_'+suffix+'/Season_'+str(season)+'/N_healpix_64'
#thedir='/sps/lsst/data/dev/pgris/Make_Cadence/Obs_minion_1016_DD_resh_new'
#pkl_file = open(thedir+'/'+name,'rb')
#thedict=pkl.load(pkl_file)

fi=thedir+'/'+name
fi_sum=thedir+'/'+name_summary
#print(thedict.keys(),thedict['dataSlice'].dtype.names)

#legend=['band','mjd','exptime','rawSeeing','FWHMeff','moon_frac','sky','kAtm','airmass','m5sigmadepth','Nexp','Ra','Dec','ditheredRA','ditheredDec']
legend=['band','mjd','exptime','rawSeeing','FWHMeff','moon_frac','sky','kAtm','airmass','m5sigmadepth','Nexp','Ra','Dec']

summa=np.load(fi_sum)

#print('Number of visits',np.asscalar(summa['Nvisits']))
toprocess=np.load(fi)
params=parameters()
   
outdir='OpSimLogs_'+suffix+'/'+opts.fieldname+'/Year_'+str(season)+'/N_healpix_64'+coadd
 
if not os.path.exists(outdir):
    os.makedirs(outdir)

outputfile  = open(outdir+'/'+prefix+'_'+opts.fieldname+'_'+str(opts.fieldid)+'.txt','wb')
#outputfile  = open(outdir+'/'+prefix+'_'+opts.fieldname+'_'+str(opts.fieldRA)+'_'+str(opts.fieldDec)+'.txt','w')

for leg in legend:
    lego=leg+' : '
    if leg == 'FWHMeff':
        lego='seeing : [was FWHMeff]'
    #print('writing',lego)
    outputfile.write('# '+lego+'\n')
outputfile.write('# end \n')
    
resultdict={}
for band in bands:
    resultdict[band]=Do_It(fieldname,band,toprocess,newsky=False,diff_time=opts.time_coadd)
       
#print('dumping')
for key in bands:
    for vals in resultdict[key]:
        outputfile.write(vals)
#print('done')


