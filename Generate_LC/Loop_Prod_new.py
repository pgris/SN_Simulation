import numpy as np
import os
from optparse import OptionParser

def batch(tab,inum):
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
    name_id='LC_Simulation_'+str(inum)
    log = dirLog + '/'+name_id+'.log'


    qsub = "qsub -P P_lsst -l sps=1,ct=03:00:00,h_vmem=16G -j y -o "+ log + " <<EOF"
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
        #print 'go',val,val.dtype
        cmd='python generate_lc.py --z '+str(val['z'])+' --fieldname '+val['fieldname']+' --fieldid '+str(val['fieldid'])+' --season '+str(val['season'])+' --sntype '+val['sntype']+' --stretch '+str(val['x1'])+' --color '+str(val['c'])+' --dirmeas '+val['dirmeas']+' --dirout '+val['dirout']+ ' --T0step '+str(val['T0step'])+' --stretch_weight '+str(val['weight_x1'])+' --color_weight '+str(val['weight_c'])+' --Opsimlog '+val['Opsimlog']+' --DayMax '+str(val['DayMax'])
       
        script.write(cmd+" \n")
    script.write("EOF" + "\n")
    script.close()
    os.system("sh "+scriptName)


parser = OptionParser()

parser.add_option("--simu_name", type="string", default="sncosmo", help="filter [%default]")
parser.add_option("--paramfile", type="string", default='', help="filter [%default]")
opts, args = parser.parse_args()
simu_name=opts.simu_name
params=np.load(opts.paramfile)
print params
        
n_per_batch=20
#n_per_batch=100
ivals=range(0,len(params),n_per_batch)
ivals=np.append(ivals,len(params))

print params.dtype
for i in range(len(ivals)-1):
    batch(params[ivals[i]:ivals[i+1]],i)

print len(params)
