from Generate_LC import *
from Observations import *
import numpy as np
from Telescope import *
#from Fit_LC import *
from Generate_Single import Generate_Single_LC
import time
import multiprocessing
import cPickle as pkl
from optparse import OptionParser
import os
from astropy.table import vstack,Table
import h5py

parser = OptionParser()

#parser.add_option("-N", "--nevts", type="int", default=100, help="filter [%default]")
parser.add_option("-f", "--fieldname", type="string", default="WFD", help="filter [%default]")
parser.add_option("-n", "--fieldid", type="int", default=309, help="filter [%default]")
parser.add_option("-s", "--season", type="int", default=1, help="filter [%default]")
parser.add_option("-z", "--z", type="float", default=0.0, help="filter [%default]")
#parser.add_option("-Z", "--zmax", type="float", default=0.1, help="filter [%default]")
#parser.add_option("--zrandom", type="string", default="yes", help="filter [%default]")
parser.add_option("-x", "--stretch", type="float", default=2.0, help="filter [%default]")
parser.add_option("-c", "--color", type="float", default=-0.2, help="filter [%default]")
parser.add_option("--stretch_weight", type="float", default=1.0, help="filter [%default]")
parser.add_option("--color_weight", type="float", default=1.0, help="filter [%default]")
parser.add_option("-t", "--sntype", type="string", default="Ia", help="filter [%default]")
parser.add_option("-d", "--dirmeas", type="string", default="None", help="filter [%default]")
#parser.add_option("-r", "--T0random", type="string", default="yes", help="filter [%default]")
parser.add_option("--min_rf_phase", type="float", default=-20.0, help="filter [%default]")
parser.add_option("--max_rf_phase", type="float", default=60., help="filter [%default]")
#parser.add_option("--T0min", type="int", default=0, help="filter [%default]")
#parser.add_option("--T0max", type="int", default=10, help="filter [%default]")
parser.add_option("--dirout", type="string", default="Light_Curves_sncosmo", help="filter [%default]")
parser.add_option("--multiproc", type="string", default="no", help="filter [%default]")
parser.add_option("--T0step", type="float", default=0.2, help="filter [%default]")
parser.add_option("--Opsimlog", type="string", default="OpSimLogs", help="filter [%default]")
parser.add_option("--DayMax", type="float", default=-1., help="filter [%default]")
#parser.add_option("-r", "--rolling", type="int", default=0, help="filter [%default]")
#parser.add_option("--nrolling", type="int", default=0, help="filter [%default]")
#parser.add_option("--merge_factor", type="int", default=0, help="filter [%default]")

opts, args = parser.parse_args()

telescope=Telescope(atmos=True,aerosol=False,airmass=1.2)

fieldid=opts.fieldid
num_season=opts.season
#N_sn=opts.nevts
z=opts.z
#zmax=opts.zmax
#zrandom=opts.zrandom
#T0random=opts.T0random
X1=opts.stretch
Color=opts.color
X1_weight=opts.stretch_weight
Color_weight=opts.color_weight
min_rf_phase=opts.min_rf_phase
max_rf_phase=opts.max_rf_phase
#T0min=opts.T0min
#T0max=opts.T0max
dirout=opts.dirout
multiproc=opts.multiproc
T0step=opts.T0step
DayMax=opts.DayMax

OpSim_Logs_dir=opts.Opsimlog
#filename=OpSim_Logs_dir+'/'+opts.dirmeas+'/Observations_'+opts.fieldname+'_'+str(fieldid)+'.txt'
filename=OpSim_Logs_dir+'/Observations_'+opts.fieldname+'_'+str(fieldid)+'.txt' 

#filename='OpSimLogs/WFD/Observations_WFD_'+str(fieldid)+'.txt'

myobs=Observations(fieldid=fieldid, filename=filename)
if T0step > 0:
    dirout+='_'+str(T0step).replace('.','_')
else:
    dirout+='_random'
dirout+='_'+str(min_rf_phase).replace('-','m')+'_'+str(max_rf_phase)

outdir='/sps/lsst/data/dev/pgris/'+dirout+'/'+opts.dirmeas+'/'+str(fieldid)+'/Season_'+str(num_season)+'/z_'+str(z)

if not os.path.exists(outdir):
    os.makedirs(outdir)


print(len(myobs.seasons),filename)

#num_season=1

myseason=myobs.seasons[num_season]

print myseason.dtype,np.min(myseason['mjd']),np.max(myseason['mjd'])

iddx=myseason['band']!='LSSTPG::u'
mysel=myseason[iddx]

#remove some of the obs points depending on the redshift

min_season=np.min(mysel['mjd'])
max_season=np.max(mysel['mjd'])

duration=max_season-min_season

time_begin=time.time()
"""
T0_vals=np.arange(min_season,max_season,0.5)
z_vals=np.arange(zmin,zmax,0.01)
"""

date_obs=min_season+20.

n_multi=8

n_multi=1
tab_X1_c = np.rec.fromrecords([(X1,X1_weight,Color,Color_weight,'')],names=['X1','X1_weight','Color','Color_weight','type'])
    #print 'aah',tab_X1_c

#n_batch=N_sn/n_multi

#n_batch=len(T0_vals)
#n_multi=len(z_vals)

#n_multi=1
print 'alors tab',len(tab_X1_c)
nbatch=len(tab_X1_c)/n_multi
nbatch+=1

if T0step > 0:
    T0_vals=np.arange(min_season,max_season,T0step)
else:
    T0_vals=[DayMax]

#print 'ooo',len(T0_vals),X1_Color_npzfile['x1_vals'],X1_Color_npzfile['x1_weights'],len(X1_Color_npzfile['x1_vals']),X1_Color_npzfile['c_vals'],X1_Color_npzfile['c_weights']/np.min(X1_Color_npzfile['c_weights']),len(X1_Color_npzfile['c_weights'])

print 'Number of DayMax',len(T0_vals),nbatch
lcs=[]

name_for_output=opts.fieldname+'_'+str(fieldid)+'_'+str(z)+'_X1_'+str(X1)+'_C_'+str(Color)

if T0step < 0:
    name_for_output+='_DayMax_'+str(DayMax)

name_out=outdir+'/'+name_for_output+'.hdf5'

if os.path.isfile(name_out):
    print 'file',name_out,'exists'
    os.remove(name_out)

itot=-1
print 'Processing ',name_out
#for T0 in T0_vals[T0min:T0max]:
for T0 in T0_vals:   
    if multiproc == 'yes':
        ival=0
        n_multi_l=n_multi
        #print 'processing',T0,T0min,T0max
        for i in range(nbatch):
       
            result_queue = multiprocessing.Queue()
            if ival+n_multi_l > len(tab_X1_c):
                n_multi_l=len(tab_X1_c)-ival
        
            for j in range(0,n_multi_l):
                p=multiprocessing.Process(name='Subprocess-'+str(j),target=Generate_Single_LC,args=(z,T0,tab_X1_c['X1'][ival],tab_X1_c['X1_weight'][ival],tab_X1_c['Color'][ival],tab_X1_c['Color_weight'][ival],myseason,telescope,j,min_rf_phase,max_rf_phase,duration,date_obs,multiproc,result_queue))
                p.start()
                ival+=1
            resultdict = {}
            for j in range(0,n_multi_l):
                resultdict.update(result_queue.get())
                
            for p in multiprocessing.active_children():
                p.join()
                
            for j in range(0,n_multi_l):
                lcs.append(resultdict[j])

    else:
        for ival in range(len(tab_X1_c)):
            itot+=1
            lc=Generate_Single_LC(z,T0,tab_X1_c['X1'][ival],tab_X1_c['X1_weight'][ival],tab_X1_c['Color'][ival],tab_X1_c['Color_weight'][ival],myseason,telescope,0,min_rf_phase,max_rf_phase,duration,date_obs,multiproc).get_lc()
            #print 'meta',ival,T0,T0min,T0max
            if itot == 0:
                lc.write(name_out, path='lc_'+str(itot), compression=True)
            else:
                lc.write(name_out, path='lc_'+str(itot), append=True,compression=True)
        

        #print 'total elapse time',time.time()-time_beginb   
      
    #lcs+=[Generate_Single_LC(z,T0,tab_X1_c['X1'][i],tab_X1_c['X1_weight'][i],tab_X1_c['Color'][i],tab_X1_c['Color_weight'][i],myseason,telescope,0,min_rf_phase,max_rf_phase,duration,date_obs).lc for i in range(len(tab_X1_c))]
    #lcs+=[Generate_Single_LC(z,T0,tab_X1_c['X1'][i],tab_X1_c['X1_weight'][i],tab_X1_c['Color'][i],tab_X1_c['Color_weight'][i],myseason,telescope,0,min_rf_phase,max_rf_phase,duration,date_obs).lc for i in range(16)]




"""
for i in range(0,n_batch):
    result_queue = multiprocessing.Queue()
    date_obs=min_season+20.+float(i)
#process=[]
    #print 'processing main',i
    #T0=T0_vals[i]
    name_for_output=opts.fieldname+'_'+str(fieldid)+'_'+str(zmin)+'_'+str(zmax)+'_X1_'+str(X1)+'_C_'+str(Color)+'_'+str(i)
    for j in range(0,n_multi):
        
        if zrandom=='yes':
            z=np.random.uniform(zmin,zmax)
        else:
            z=zmin


        if T0random == 'yes':
            T0 = np.random.uniform(min_season,max_season)
        else:
            T0=0
            T0=date_obs-delta[j]
            T0=date_obs-max_rf_phase*(1.+z)

        window_max=date_obs
        window_min=date_obs-2.*np.abs(T0-date_obs)
        window_min=T0+min_rf_phase*(1.+z)

        #print 'hello',i,j,date_obs,min_season,window_min, window_max, T0,min_rf_phase,max_rf_phase

        #z=z_vals[j]

        X1_val=X1
        Color_val=Color
        if X1 == -999.:
            X1_val=np.random.choice(X1_Color_npzfile['x1_vals'],1,p=X1_Color_npzfile['x1_weights'])[0]
        if Color==-999.:
            Color_val=np.random.choice(X1_Color_npzfile['c_vals'],1,p=X1_Color_npzfile['c_weights'])[0]

        #print 'hello I will process',X1_val,Color_val
        #print 'Processing',j,T0,z
        p=multiprocessing.Process(name='Subprocess-'+str(i),target=Generate_Single_LC,args=(z,T0,X1_val,Color_val,myseason,telescope,j,min_rf_phase,max_rf_phase,duration,date_obs,result_queue))
    #process.append(p)
        p.start()
    
    resultdict = {}
    for j in range(0,n_multi):
        resultdict.update(result_queue.get())

    for p in multiprocessing.active_children():
        p.join()

    
    pkl_file = open(outdir_obs+'/'+name_for_output+'.pkl','wb')
    lcs=[]
    for j in range(0,n_multi):
        #print 'hello there',resultdict[j].dtype
        outdict={}
        outdict['metadata']=resultdict[j][0]
        outdict['lc']=resultdict[j][1]
        lcs.append(outdict)

    pkl.dump(lcs, pkl_file)
        
    pkl_file.close()
"""

print 'total elapse time',time.time()-time_begin
