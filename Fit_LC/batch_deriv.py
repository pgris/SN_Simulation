from optparse import OptionParser
import os
import numpy as np
import cPickle as pkl

def batch(atype,idmin,idmax):

    cwd = os.getcwd()
    dirScript= cwd + "/scripts_deriv_mb"

    if not os.path.isdir(dirScript) :
        os.makedirs(dirScript)
    
    dirLog = cwd + "/logs_deriv_mb"
    if not os.path.isdir(dirLog) :
        os.makedirs(dirLog)  

    name_id=atype.strip()+'_'+str(idmin)+'_'+str(idmax)
    log = dirLog + '/'+name_id+'.log'
    qsub = "qsub -P P_lsst -l sps=1,ct=4:00:00,h_vmem=16G -j y -o "+ log + " <<EOF"

    
    scriptName = dirScript+'/'+name_id+'.sh'

    script = open(scriptName,"w")
    script.write(qsub + "\n")
    script.write("#!/usr/local/bin/bash\n")
    script.write(" cd " + cwd + "\n")
    script.write("bash " + "\n")
            #script.write("ls /usr/local/grid/emi-3/WN/SL6_64/3.10.0-1.2/usr/lib64/" + "\n")
    script.write(" source setups_cosmomaf.sh\n")

    cmd='python Make_Derivatives_mb.py --type '+atype+' --idmin '+str(idmin)+' --idmax '+str(idmax)+' --save_file yes'
    script.write(cmd+" \n")
    script.write("EOF" + "\n")
    script.close()
    os.system("sh "+scriptName)

parser = OptionParser()
parser.add_option("--type", type="string", default='low_z', help="filter [%default]") 
opts, args = parser.parse_args()

atype=opts.type

"""
pkl_file = open('../Generate_LC/Map_X1_C.pkl','rb')
tab_X1_c=pkl.load(pkl_file)

r=[]
idx= tab_X1_c['type'] == atype

for z in np.arange(0.01,1.1,0.01):
    for val in tab_X1_c[idx]:
        r.append((z,val['X1'],val['Color']))
"""
pkl_file = open('../../SNSim_Simulation/Files_Params_'+atype+'.pkl','rb')
all_params=pkl.load(pkl_file)

nbatch=100
n_per_batch=len(all_params)/100

print 'hello',len(all_params)
for i in range(nbatch+1):
    idmin=i*n_per_batch
    idmax=idmin+n_per_batch
    if idmax > len(all_params)-1:
        idmax=len(all_params)-1

    batch(atype,idmin,idmax)
