import os

def batch(simu,fieldname,season,coadd,time_coadd):
    
    cwd = os.getcwd()
    dirScript= cwd + "/scripts"

    if not os.path.isdir(dirScript) :
        os.makedirs(dirScript)
    
    dirLog = cwd + "/logs"
    if not os.path.isdir(dirLog) :
        os.makedirs(dirLog)    
        
    
    name_id=fieldname+'_'+simu+'_'+str(season)+'_'+coadd
    log = dirLog + '/'+name_id+'.log'

    qsub = "qsub -P P_lsst -l sps=1,ct=30:00:00,h_vmem=16G -j y -o "+ log + " -pe multicores 8 <<EOF"

    scriptName = dirScript+'/'+name_id+'.sh'
    
    script = open(scriptName,"w")
    script.write(qsub + "\n")
    script.write("#!/usr/local/bin/bash\n")
    script.write(" cd " + cwd + "\n")
    script.write(" source setups_cosmomaf.sh\n")

    cmd='python Loop_Dump_newcad.py --simu '+simu+' --fieldname '+fieldname+' --season '+str(season)+' --coadd '+coadd+' --time_coadd '+str(time_coadd)

    script.write(cmd+" \n")
    script.write("EOF" + "\n")
    script.close()
    os.system("sh "+scriptName)
    #time.sleep(1)







fieldname='WFD'
season=2
coadd='yes'
time_coadd=24.*3600.

for simu in ['feature_baseline_10yrs','feature_rolling_half_mask_10yrs','feature_rolling_twoThird_10yrs','alt_sched','alt_sched_rolling']:
    batch(simu,fieldname,season,coadd,time_coadd)
