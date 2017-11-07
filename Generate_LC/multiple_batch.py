from optparse import OptionParser
import os
import numpy as np

parser = OptionParser()
parser.add_option("-z", "--zmin", type="float", default=0.0, help="filter [%default]")
parser.add_option("-Z", "--zmax", type="float", default=1.1, help="filter [%default]")
parser.add_option("--zstep", type="float", default=0.01, help="filter [%default]")
parser.add_option("-f", "--fieldname", type="string", default='WFD', help="filter [%default]")
parser.add_option("-i", "--fieldid", type="int", default=309, help="filter [%default]")
parser.add_option("-s", "--season", type="int", default=-1, help="filter [%default]")
parser.add_option("-t", "--sntype", type="string", default='Ia', help="filter [%default]")
#parser.add_option("-d", "--dbFile", type="string",default='None', help="dbFile to process [%default]")
parser.add_option("-x", "--stretch", type="float", default=2.0, help="filter [%default]")
parser.add_option("-c", "--color", type="float", default=-0.2, help="filter [%default]")
parser.add_option("-d", "--dirmeas", type="string", default="None", help="filter [%default]")
#parser.add_option("-r", "--T0random", type="string", default="No", help="filter [%default]")

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
#T0random=opts.T0random

for z in np.arange(zmin,zmax+zstep,zstep):
    cmd='python batch.py --z '+str(z)+' --fieldname '+fieldname+' --fieldid '+str(fieldid)+' --season '+str(season)+' --sntype '+sntype+' --stretch '+str(stretch)+' --color '+str(color)+' --dirmeas '+dirmeas
    print 'executing',cmd
    os.system(cmd)
