import os
from optparse import OptionParser

parser = OptionParser()

parser.add_option("-f", "--fieldname", type="string", default='DD', help="filter [%default]")
parser.add_option("-i", "--fieldid", type="int", default=744, help="filter [%default]")
parser.add_option("-s", "--season", type="int", default=0, help="filter [%default]")
parser.add_option("-t", "--sntype", type="string", default='Ia', help="filter [%default]")
#parser.add_option("-d", "--dbFile", type="string",default='None', help="dbFile to process [%default]")
parser.add_option("-x", "--stretch", type="float", default=2.0, help="filter [%default]")
parser.add_option("-c", "--color", type="float", default=-0.2, help="filter [%default]")
parser.add_option("-d", "--dirmeas", type="string", default="DD", help="filter [%default]")
parser.add_option("--NT0", type="int", default=-1, help="filter [%default]")
parser.add_option("--dirout", type="string", default="None", help="filter [%default]")
parser.add_option("--T0min", type="int", default=0, help="filter [%default]")
parser.add_option("--T0max", type="int", default=20, help="filter [%default]")

opts, args = parser.parse_args()

fieldname=opts.fieldname
fieldid=opts.fieldid
season=opts.season
sntype=opts.sntype
stretch=opts.stretch
color=opts.color
dirmeas=opts.dirmeas
NT0=opts.NT0
dirout=opts.dirout
T0min=opts.T0min
T0max=opts.T0max


zstep=0.01
nz=2

for j in range(10):
    T0min=j*100
    T0max=T0min+100

    for i in range(0,50):
        zmin=0.01+0.01*float(i*nz)
        zmax=zmin+nz*zstep
    #cmd='python multiple_batch.py --zmin '+str(zmin)+' --zmax '+str(zmax)+' --zstep '+str(zstep)+' --fieldname '+fieldname+' --fieldid '+str(fieldid)+' --season '+str(season)+' --stretch '+str(stretch)+' --color '+str(color)+' --dirmeas '+dirmeas+' --NT0 '+str(NT0)+' --dirout '+dirout
        cmd='python batch.py --zmin '+str(zmin)+' --zmax '+str(zmax)+' --zstep '+str(zstep)+' --fieldname '+fieldname+' --fieldid '+str(fieldid)+' --season '+str(season)+' --stretch '+str(stretch)+' --color '+str(color)+' --dirmeas '+dirmeas+' --dirout '+dirout +' --T0min '+str(T0min)+' --T0max '+str(T0max)
        print cmd
        os.system(cmd)

    #break
