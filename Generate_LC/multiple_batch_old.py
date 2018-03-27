from optparse import OptionParser
import os
import numpy as np
from Observations import *

parser = OptionParser()
parser.add_option("-z", "--zmin", type="float", default=0.1, help="filter [%default]")
parser.add_option("-Z", "--zmax", type="float", default=0.2, help="filter [%default]")
parser.add_option("--zstep", type="float", default=0.025, help="filter [%default]")
parser.add_option("-f", "--fieldname", type="string", default='WFD', help="filter [%default]")
parser.add_option("-i", "--fieldid", type="int", default=309, help="filter [%default]")
parser.add_option("-s", "--season", type="int", default=-1, help="filter [%default]")
parser.add_option("-t", "--sntype", type="string", default='Ia', help="filter [%default]")
#parser.add_option("-d", "--dbFile", type="string",default='None', help="dbFile to process [%default]")
parser.add_option("-x", "--stretch", type="float", default=2.0, help="filter [%default]")
parser.add_option("-c", "--color", type="float", default=-0.2, help="filter [%default]")
parser.add_option("-d", "--dirmeas", type="string", default="None", help="filter [%default]")
parser.add_option("--NT0", type="int", default=10, help="filter [%default]")
parser.add_option("--dirout", type="string", default="None", help="filter [%default]")
parser.add_option("--X1C", type="string", default="", help="filter [%default]")
parser.add_option("--T0step", type="float", default=0.2, help="filter [%default]")

opts, args = parser.parse_args()

zmin=opts.zmin
zmax=opts.zmax
zstep=opts.zstep
fieldname=opts.fieldname
fieldid=opts.fieldid
season=opts.season
sntype=opts.sntype
stretch=opts.stretch
color=opts.color
dirmeas=opts.dirmeas
NT0=opts.NT0
dirout=opts.dirout
X1C=opts.X1C
T0step=opts.T0step

#T0random=opts.T0random

OpSim_Logs_dir=os.getenv('OPSIM_LOGS')
filename=OpSim_Logs_dir+'/'+opts.dirmeas+'/Observations_'+opts.fieldname+'_'+str(fieldid)+'.txt'

myobs=Observations(fieldid=fieldid, filename=filename)

myseason=myobs.seasons[season]

#print myseason.dtype,np.min(myseason['mjd']),np.max(myseason['mjd'])

iddx=myseason['band']!='LSSTPG::u'
mysel=myseason[iddx]

min_season=np.min(mysel['mjd'])
max_season=np.max(mysel['mjd'])

T0_vals=np.arange(min_season,max_season,T0step)
N_T0=NT0
if NT0==-1:
    N_T0=len(T0_vals)

T0min=0
T0max=T0min+N_T0
print 'alors',T0min,T0max,N_T0
while T0min < len(T0_vals):
    if T0max > len(T0_vals):
        T0max = len(T0_vals)-1
    cmd='python batch.py --zmin '+str(zmin)+' --zmax '+str(zmax)+' --zstep '+str(zstep)+' --fieldname '+fieldname+' --fieldid '+str(fieldid)+' --season '+str(season)+' --sntype '+sntype+' --stretch '+str(stretch)+' --color '+str(color)+' --dirmeas '+dirmeas+' --T0min '+str(T0min)+' --T0max '+str(T0max)+' --dirout '+dirout+' --T0step '+str(T0step)
    if X1C != '':
        cmd+=' --X1C '+X1C
    print 'executing',cmd
    os.system(cmd)
    T0min+=N_T0
    T0max=T0min+N_T0
