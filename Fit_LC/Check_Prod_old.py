import numpy as np
import glob
from Load_X1_C import Load_X1_Color
import os
from Observations import *

def batch(tab, ibatch,multiproc):

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

    name_id=fieldname+'_'+str(fieldid)+'_season_'+str(season)+'_z_'+str(tab['z'][0])+'_x1_'+str(stretch)+'_c_'+str(color)+'_T0min_'+str(T0min)+'_T0max_'+str(T0max)+'_'+str(ibatch)+'_'+tab['simu'][0]
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


fieldname='DD'
fieldid=290
season=0
stretch=-999.0
color=-999.0
T0step=0.5
simu='sncosmo'

thedir='/sps/lsst/data/dev/pgris/Fitted_Light_Curves_'+simu+'_'+str(T0step).replace('.','_')+'/'+fieldname+'/'+str(fieldid)+'/Season_'+str(season)

z=[0.1]
z+=[p for p in np.arange(0.025,1.1,0.025)]


#tab_X1_C=Get_X1_Color_Distribution(stretch,color)
tab_X1_C=Load_X1_Color(stretch,color).tab
torepro=[]
for zval in z:
    files=glob.glob(thedir+'/z_'+str(zval)+'/*')
    #print zval,len(files)
    if zval < 0.1:
        sel=tab_X1_C['low_z']
    else:
       sel=tab_X1_C['high_z']

    if len(files) != len(sel):
        print 'There are missing events!!!',len(files),len(sel),zval
        for x1c in sel:
            found = False
            what=fieldname+'_'+str(fieldid)+'_'+str(zval)+'_X1_'+str(np.asscalar(x1c['x1']))+'_C_'+str(np.asscalar(x1c['c']))
            #print 'lokking for',what
            for fi in files:
                if fi.count(what) >0: 
                    found=True
                    break
            if found is False:
                print 'probleme here',T0step,zval,x1c
                """
                cmd='python multiple_batch.py --fieldname '+fieldname+' --fieldid '+str(fieldid)+' --z '+str(zval)
                cmd+=' --stretch '+str(np.asscalar(x1c['x1']))+' --color '+str(np.asscalar(x1c['c']))
                cmd+=' --dirmeas Light_Curves_'+simu+'_'+str(T0step).replace('.','_')+' --season '+str(season)+' --simulator '+simu+' --dirout Fitted_Light_Curves_'+simu+'_'+str(T0step).replace('.','_')+' --n_per_batch '+str(1)+' --multiproc no'
                print cmd
                """
                
                OpSim_Logs_dir=os.getenv('OPSIM_LOGS')
                filename=OpSim_Logs_dir+'/'+fieldname+'/Observations_'+fieldname+'_'+str(fieldid)+'.txt'
                myobs=Observations(fieldid=fieldid, filename=filename)
                myseason=myobs.seasons[season]

                #print myseason.dtype,np.min(myseason['mjd']),np.max(myseason['mjd'])

                iddx=myseason['band']!='LSSTPG::u'
                mysel=myseason[iddx]
                
                min_season=np.min(mysel['mjd'])
                max_season=np.max(mysel['mjd'])

                T0_vals=np.arange(min_season,max_season,T0step)

                N_T0=len(T0_vals)
                
                T0min=0
                T0max=T0min+N_T0
                torepro.append((fieldname,fieldid,zval,season,simu,T0step,np.asscalar(x1c['x1']),np.asscalar(x1c['c']),'Fitted_Light_Curves_'+simu+'_'+str(T0step).replace('.','_'),'Light_Curves_'+simu+'_'+str(T0step).replace('.','_'),T0min,T0max,'Ia'))
                

forbatch=np.rec.fromrecords(torepro,names=['fieldname','fieldid','z','season','simu','T0step','X1','Color','dirout','dirmeas','T0min','T0max','sntype'])

if len(forbatch) > 0:
    print forbatch
    
    n_for_batch=5
    
    ibatch=range(0,len(forbatch),5)
    ibatch=np.append(ibatch,len(forbatch))

    print ibatch
    for i in range(len(ibatch)-1):
        print forbatch[ibatch[i]:ibatch[i+1]]
        print 'done'
        batch(forbatch[ibatch[i]:ibatch[i+1]],i,multiproc='yes')

else:
    print 'Production ok'
