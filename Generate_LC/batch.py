import numpy as np
from optparse import OptionParser
import os
import time

parser = OptionParser()
parser.add_option("--zmin", type="float", default=0.0, help="filter [%default]")
parser.add_option("--zmax", type="float", default=1.2, help="filter [%default]")
parser.add_option("--zstep", type="float", default=0.1, help="filter [%default]")
#parser.add_option("-N", "--nevts", type="int", default=10, help="filter [%default]")
parser.add_option("-m", "--model", type="string", default='salt2-extended', help="filter [%default]")
parser.add_option("-v", "--version", type="string", default='', help="filter [%default]")
parser.add_option("-f", "--fieldname", type="string", default='WFD', help="filter [%default]")
parser.add_option("-i", "--fieldid", type="int", default=309, help="filter [%default]")
parser.add_option("-s", "--season", type="int", default=-1, help="filter [%default]")
parser.add_option("-t", "--sntype", type="string", default='Ia', help="filter [%default]")
#parser.add_option("-d", "--dbFile", type="string",default='None', help="dbFile to process [%default]")
parser.add_option("-x", "--stretch", type="float", default=2.0, help="filter [%default]")
parser.add_option("-c", "--color", type="float", default=-0.2, help="filter [%default]")
parser.add_option("-d", "--dirmeas", type="string", default="None", help="filter [%default]")
#parser.add_option("-r", "--T0random", type="string", default="No", help="filter [%default]")
#parser.add_option("--zrandom", type="string", default="yes", help="filter [%default]")
parser.add_option("--T0min", type="int", default=0, help="filter [%default]")
parser.add_option("--T0max", type="int", default=10, help="filter [%default]")
parser.add_option("--dirout", type="string", default="Light_Curves_sncosmo", help="filter [%default]")

opts, args = parser.parse_args()

zmin=opts.zmin
zmax=opts.zmax
zstep=opts.zstep
#nevts=opts.nevts
fieldname=opts.fieldname
fieldid=opts.fieldid
season=opts.season
sntype=opts.sntype
stretch=opts.stretch
color=opts.color
dirmeas=opts.dirmeas
#T0random=opts.T0random
#zrandom=opts.zrandom
T0min=opts.T0min
T0max=opts.T0max
dirout=opts.dirout

cwd = os.getcwd()
dirScript= cwd + "/scripts"

if not os.path.isdir(dirScript) :
    os.makedirs(dirScript)
    
dirLog = cwd + "/logs"
if not os.path.isdir(dirLog) :
    os.makedirs(dirLog)    
    
name_id=fieldname+'_'+str(fieldid)+'_'+str(zmin)+'_'+str(zmax)+'_season_'+str(opts.season)+'_x1_'+str(stretch)+'_c_'+str(color)+'_T0min_'+str(T0min)+'_T0max_'+str(T0max)
log = dirLog + '/'+name_id+'.log'


qsub = "qsub -P P_lsst -l sps=1,ct=03:00:00,h_vmem=16G -j y -o "+ log + " -pe multicores 8 <<EOF"
scriptName = dirScript+'/'+name_id+'.sh'

script = open(scriptName,"w")
script.write(qsub + "\n")
script.write("#!/usr/local/bin/bash\n")
script.write(" cd " + cwd + "\n")
script.write("bash " + "\n")
            #script.write("ls /usr/local/grid/emi-3/WN/SL6_64/3.10.0-1.2/usr/lib64/" + "\n")
script.write(" source setups_cosmomaf.sh\n")
for z in np.arange(zmin,zmax+zstep,zstep):
    cmd='python generate_lc.py --z '+str(z)+' --fieldname '+fieldname+' --fieldid '+str(fieldid)+' --season '+str(season)+' --sntype '+sntype+' --stretch '+str(stretch)+' --color '+str(color)+' --dirmeas '+dirmeas+' --T0min '+str(T0min)+' --T0max '+str(T0max)+' --dirout '+dirout
    script.write(cmd+" \n")
script.write("EOF" + "\n")
script.close()
os.system("sh "+scriptName)
time.sleep(1)
