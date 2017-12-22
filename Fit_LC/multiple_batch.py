from optparse import OptionParser
import os
import numpy as np
import glob
import time

def batch(tab, ibatch,simu_name,multiproc):

    cwd = os.getcwd()
    dirScript= cwd + "/scripts"

    if not os.path.isdir(dirScript) :
        os.makedirs(dirScript)
    
    dirLog = cwd + "/logs"
    if not os.path.isdir(dirLog) :
        os.makedirs(dirLog)    
    
    fieldid=tab['fieldid'][0]
    season=tab['season'][0]
    stretch=tab['X1'][0]
    color=tab['Color'][0]
    T0min=tab['T0min'][0]
    T0max=tab['T0max'][0]

    name_id=fieldname+'_'+str(fieldid)+'_season_'+str(opts.season)+'_x1_'+str(stretch)+'_c_'+str(color)+'_T0min_'+str(T0min)+'_T0max_'+str(T0max)+'_'+str(ibatch)+'_'+simu_name
    log = dirLog + '/'+name_id+'.log'

    if multiproc == 'yes':
        qsub = "qsub -P P_lsst -l sps=1,ct=30:00:00,h_vmem=16G -j y -o "+ log + " -pe multicores 8 <<EOF"
    else:
        qsub = "qsub -P P_lsst -l sps=1,ct=30:00:00,h_vmem=16G -j y -o "+ log + " <<EOF"
 
    scriptName = dirScript+'/'+name_id+'.sh'
    
    script = open(scriptName,"w")
    script.write(qsub + "\n")
    script.write("#!/usr/local/bin/bash\n")
    script.write(" cd " + cwd + "\n")
    script.write("bash " + "\n")
            #script.write("ls /usr/local/grid/emi-3/WN/SL6_64/3.10.0-1.2/usr/lib64/" + "\n")
    script.write(" source setups_cosmomaf.sh\n")

    for val in tab:
        cmd='python fit_lcs.py --z '+str(val['z'])+' --fieldname '+val['fieldname']+' --fieldid '+str(val['fieldid'])+' --season '+str(val['season'])+' --sntype '+val['sntype']+' --stretch '+str(val['X1'])+' --color '+str(val['Color'])+' --dirmeas '+val['dirmeas']+' --T0min '+str(val['T0min'])+' --T0max '+str(val['T0max'])+' --dirout '+val['dirout']+' --multiproc '+multiproc
        script.write(cmd+" \n")
    script.write("EOF" + "\n")
    script.close()
    os.system("sh "+scriptName)
    #time.sleep(1)



parser = OptionParser()
#parser.add_option("-z", "--z", type="float", default=0.0, help="filter [%default]")
parser.add_option("-f", "--fieldname", type="string", default='WFD', help="filter [%default]")
parser.add_option("-i", "--fieldid", type="int", default=309, help="filter [%default]")
parser.add_option("-s", "--season", type="int", default=0, help="filter [%default]")
parser.add_option("--z", type="float", default=-1.0, help="redshift [%default]")
parser.add_option("-t", "--sntype", type="string", default='Ia', help="filter [%default]")
#parser.add_option("-d", "--dbFile", type="string",default='None', help="dbFile to process [%default]")
parser.add_option("-x", "--stretch", type="float", default=2.0, help="filter [%default]")
parser.add_option("-c", "--color", type="float", default=-0.2, help="filter [%default]")
parser.add_option("-d", "--dirmeas", type="string", default="None", help="filter [%default]")
parser.add_option("--simulator", type="string", default="snscosmo", help="filter [%default]")
parser.add_option("--dirout", type="string", default="Fitted_Light_Curves", help="filter [%default]")
parser.add_option("--n_per_batch", type="int", default="10", help="filter [%default]")
#parser.add_option("-r", "--T0random", type="string", default="No", help="filter [%default]")
parser.add_option("--multiproc", type="string", default='yes', help="filter [%default]")

opts, args = parser.parse_args()


fieldname=opts.fieldname
fieldid=opts.fieldid
season=opts.season
stretch=opts.stretch
color=opts.color
dirmeas=opts.dirmeas
sntype=opts.sntype
simu_name=opts.simulator
dirout=opts.dirout
#T0random=opts.T0random
z=opts.z
nfiles_per_batch=opts.n_per_batch
multiproc=opts.multiproc

main_dir_in='/sps/lsst/data/dev/pgris/'
if z == -1.0:
    files=glob.glob(main_dir_in++dirmeas+'/'+fieldname+'/'+str(fieldid)+'/Season_'+str(season)+'/'+fieldname+'_'+str(fieldid)+'_*_X1_'+str(stretch)+'_C_'+str(color)+'*.pkl')
else:
    files=glob.glob(main_dir_in+dirmeas+'/'+fieldname+'/'+str(fieldid)+'/Season_'+str(season)+'/'+fieldname+'_'+str(fieldid)+'_'+str(z)+'_X1_'+str(stretch)+'_C_'+str(color)+'*.pkl')


r=[]

for fi in files:
    """
    fisplit=fi.split('_')
    T0min=fisplit[-2]
    T0max=fisplit[-1].split('.')[0]
    """
    
    spl=fi.split('/')
    
    splb=spl[-1].split('_')
    
    z=splb[2]
    stretch=float(splb[4])
    color=float(splb[6])
    T0min=splb[7]
    T0max=splb[8].split('.')[0]
    
    r.append((dirmeas,fieldname,fieldid,z,season,stretch,color,sntype,T0min,T0max,dirout))

params= np.rec.fromrecords(r,names=['dirmeas','fieldname','fieldid','z','season','X1','Color','sntype','T0min','T0max','dirout'])
print len(params)



nbatches=len(params)/nfiles_per_batch
if len(params)%nfiles_per_batch > 0:
    nbatches+=1

for i in range(nbatches):
    imin=i*nfiles_per_batch
    imax=imin+nfiles_per_batch
    if imax >= len(params):
        imax=len(params)
    batch(params[imin:imax],i,simu_name,multiproc)

"""
cmd='python batch.py --z '+str(z)+' --fieldname '+fieldname+' --fieldid '+str(fieldid)+' --season '+str(season)+' --sntype '+sntype+' --stretch '+str(stretch)+' --color '+str(color)+' --dirmeas '+dirmeas+' --T0min '+str(T0min)+' --T0max '+str(T0max)
    print 'executing',cmd
    os.system(cmd)
"""
