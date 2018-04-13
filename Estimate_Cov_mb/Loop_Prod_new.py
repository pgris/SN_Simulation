import numpy as np
import os
from optparse import OptionParser

def batch(tab,inum,multiproc='yes'):
    cwd = os.getcwd()
    dirScript= cwd + "/scripts"

    if not os.path.isdir(dirScript) :
        os.makedirs(dirScript)
    
    dirLog = cwd + "/logs"
    if not os.path.isdir(dirLog) :
        os.makedirs(dirLog)    
    
    """
    name_id=fieldname+'_'+str(fieldid)+'_'+str(zmin)+'_'+str(zmax-zstep)+'_season_'+str(opts.season)+'_x1_'+stretch.replace(',','_')+'_c_'+color.replace(',','_')+'_T0min_'+str(T0min)+'_T0max_'+str(T0max)
    log = dirLog + '/'+name_id+'.log'
    """
    name_id='LC_Fit_'+str(inum)
    log = dirLog + '/'+name_id+'.log'

    if multiproc == 'yes':
        qsub = "qsub -P P_lsst -l os=sl6,sps=1,ct=30:00:00,h_vmem=16G -j y -o "+ log + " -pe multicores 8 <<EOF"
    else:
        qsub = "qsub -P P_lsst -l os=sl6,sps=1,ct=30:00:00,h_vmem=16G -j y -o "+ log + " <<EOF"


    scriptName = dirScript+'/'+name_id+'.sh'

    script = open(scriptName,"w")
    script.write(qsub + "\n")
    script.write("#!/usr/local/bin/bash\n")
    script.write(" cd " + cwd + "\n")
    script.write("bash " + "\n")
            #script.write("ls /usr/local/grid/emi-3/WN/SL6_64/3.10.0-1.2/usr/lib64/" + "\n")
    script.write(" source setups_cosmomaf.sh\n")
#for z in np.arange(zmin,zmax-zstep,zstep):

    for val in tab:
        cmd='python estimate_cov_mb.py --z '+str(val['z'])+' --fieldname '+val['fieldname']+' --fieldid '+str(val['fieldid'])+' --season '+str(val['season'])+' --sntype '+val['sntype']+' --stretch '+str(val['x1'])+' --color '+str(val['c'])+' --dirmeas '+val['dirmeas']+' --dirout '+val['dirout']+' --multiproc '+multiproc+' --DayMax '+str(val['DayMax'])
       
        script.write(cmd+" \n")
    script.write("EOF" + "\n")
    script.close()
    os.system("sh "+scriptName)

parser = OptionParser()
parser.add_option("--paramfile", type="string", default='', help="filter [%default]")
parser.add_option("--min_rf_phase", type="float", default=-20.0, help="filter [%default]")
parser.add_option("--max_rf_phase", type="float", default=60., help="filter [%default]")

opts, args = parser.parse_args()

paramfile=opts.paramfile
min_rf_phase=opts.min_rf_phase
max_rf_phase=opts.max_rf_phase

if paramfile == '':
    print('You have to give a parameter filename')
else:
    paramsb=np.load(paramfile)
    
    names=[name for name in paramsb.dtype.names if name != 'dirout' and name!= 'dirmeas']
    r=[]
    for par in paramsb:
        lo=[par[name] for name in names]
        if par['T0step'] >0. :
            ddout=par['dirout']+'_'+str(par['T0step']).replace('.','_')+'_'+str(min_rf_phase).replace('-','m')+'_'+str(max_rf_phase)
        else:
            ddout=par['dirout']+'_random_'+str(min_rf_phase).replace('-','m')+'_'+str(max_rf_phase)
        dirout=ddout
        dirmeas='Fitted_'+ddout
        dirout=dirmeas+'_covmb'
        lo.append(dirout)
        lo.append(dirmeas)
        r.append(lo)
    names+=['dirout','dirmeas']
    
    params=np.rec.fromrecords(r,names=names)
    print params
    
    
    n_per_batch=40
    multiproc='yes'
    #n_per_batch=400
    #multiproc='no'
    ivals=range(0,len(params),n_per_batch)
    ivals=np.append(ivals,len(params))

    print params.dtype
    for i in range(len(ivals)-1):
        batch(params[ivals[i]:ivals[i+1]],i,multiproc=multiproc)
    
    print len(params)
    
