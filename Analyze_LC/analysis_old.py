import numpy as np
from LC_Ana import *
from optparse import OptionParser
from ID import *
import cPickle as pkl
import matplotlib.pyplot as plt
from scipy.spatial import distance

parser = OptionParser()

parser.add_option("-f", "--fieldname", type="string", default='DD', help="filter [%default]")
parser.add_option("-i", "--fieldid", type="int", default=744, help="filter [%default]")
parser.add_option("-x", "--stretch", type="float", default=-999., help="filter [%default]")
parser.add_option("-c", "--color", type="float", default=-999., help="filter [%default]")
parser.add_option("-d", "--dirmeas", type="string", default="DD", help="filter [%default]")
parser.add_option("--dirobs", type="string", default="DD", help="filter [%default]")
parser.add_option("--min_rf_phase", type="float", default=-20.0, help="filter [%default]")
parser.add_option("--max_rf_phase", type="float", default=60., help="filter [%default]")

opts, args = parser.parse_args()

dict_data={}

fieldid=opts.fieldid
fieldname=opts.fieldname
thedir=opts.dirmeas
thedir_obs=opts.dirobs
X1=opts.stretch
Color=opts.color
min_rf_phase=opts.min_rf_phase
max_rf_phase=opts.max_rf_phase
#zmax=0.7
zmax=1.1
z=[0.01]
step=0.025
z+=[val for val in np.arange(step,zmax,step)]
#z=[0.1,0.6]
print 'hello z',z
#z=[0.375]

add='_0_5'
add+='_'+str(min_rf_phase).replace('-','m')+'_'+str(max_rf_phase)

simul_name='snsim'
fieldids=[290,744,1427,2412,2786]
for fieldid in fieldids:
    for seas in range(10):
#for seas in range(2,3):
    #for T0step in [0.2,0.3,0.5,1.0]:
        for T0step in [0.5]:
            strs=str(T0step).replace('.','_')
            """
            for fieldid in [290,744,1427,2412,2786]:
        
                for fieldid in [100,101,102,103]:
                for fieldid in [6084]:
                if 101 in fieldid:
                add='_feature_baseline_10yrs_0_5'
                if fieldid == 290:
                add+='_b'
            """
            for (x1,c) in [(-2.0,0.2),(2.0,-0.2),(0.0,0.0)]:
            #for (x1,c) in [(-2.0,0.2)]:
                for zval in z:
                    
                    tag_name='DD_'+simul_name+'_'+str(fieldid)+'_'+str(seas+1)+'_'+str(x1)+'_'+str(c)+'_'+str(zval)
                    dict_data[tag_name]=Id(thedir='Fitted_Light_Curves_'+simul_name+add,thedir_obs=thedir_obs,fieldname=fieldname,fieldid=fieldid,X1=x1,Color=c,season=seas,z=zval,colorfig='k',lstyle='-',T0step=T0step,simul_name=simul_name)
                    #print('boo',zval,tag_name)
                    #dict_data[fieldname+'_sncosmo_'+str(fieldid)+'_'+str(seas+1)+'_'+str(x1)+'_'+str(c)+'_'+str(zval)]=Id(thedir='Fitted_Light_Curves_sncosmo_alt_sched_rolling/Year_2/N_healpix_64_coadd_0_5',thedir_obs=thedir_obs,fieldname=fieldname,fieldid=fieldid,X1=x1,Color=c,season=seas,z=zval,colorfig='k',lstyle='-',T0step=T0step) 
            #dict_data['DD_snsim_'+str(fieldid)+'_'+str(seas+1)+'_'+strs]=Id(thedir='Fitted_Light_Curves_snsim_test',thedir_obs=thedir_obs,fieldname=fieldname,fieldid=fieldid,X1=X1,Color=Color,season=seas,z=z,colorfig='k',lstyle='-',T0step=T0step)
    #dict_data['DD_sncosmo'+str(fieldid)+'_'+str(seas+1)]=Id(thedir='Fitted_Light_Curves_sncosmo_test_0_2_b',thedir_obs=thedir_obs,fieldname=fieldname,fieldid=fieldid,X1=X1,Color=Color,season=seas,z=z,colorfig='k',lstyle='-')
    #dict_data['DD_sncosmob'+str(fieldid)+'_'+str(seas+1)]=Id(thedir='Fitted_Light_Curves_sncosmo_testb_0_2_b',thedir_obs=thedir_obs,fieldname=fieldname,fieldid=fieldid,X1=X1,Color=Color,season=seas,z=z,colorfig='k',lstyle='-')
    """
    dict_data['DD_all'+str(fieldid)+'_'+str(seas+1)]=Id(thedir=thedir,thedir_obs=thedir_obs,fieldname=fieldname,fieldid=fieldid,X1=X1,Color=Color,season=seas,z=z,colorfig='k',lstyle='-')
    dict_data['DD_minus'+str(fieldid)+'_'+str(seas+1)]=Id(thedir=thedir+'_b',thedir_obs=thedir_obs,fieldname=fieldname,fieldid=fieldid,X1=X1,Color=Color,season=seas,z=z,colorfig='g',lstyle='-')
    dict_data['DD_third'+str(fieldid)+'_'+str(seas+1)]=Id(thedir='Fitted_Light_Curves_snsim_0_2_b',thedir_obs=thedir_obs,fieldname=fieldname,fieldid=fieldid,X1=X1,Color=Color,season=seas,z=z,colorfig='r',lstyle='--')
    dict_data['DD_fourth'+str(fieldid)+'_'+str(seas+1)]=Id(thedir='Fitted_Light_Curves_snsim_0_2',thedir_obs=thedir_obs,fieldname=fieldname,fieldid=fieldid,X1=X1,Color=Color,season=seas,z=z,colorfig='b',lstyle='--')
    """
    #dict_data['DD_medium'+str(fieldid)+'_'+str(seas+1)]=Id(thedir=thedir,thedir_obs=thedir_obs,fieldname=fieldname,fieldid=fieldid,X1=0.0,Color=0.0,season=seas,colorfig='k') 
    #dict_data['DD_faint'+str(fieldid)+'_'+str(seas+1)]=Id(thedir=thedir,thedir_obs=thedir_obs,fieldname=fieldname,fieldid=fieldid,X1=-2.0,Color=0.2,season=seas,colorfig='k') 

LC_Ana(dict_data,zmin=0.0,zmax=1.2,bin_z=0.025,data_dir='/sps/lsst/data/dev/pgris',fieldids=fieldids)


   
