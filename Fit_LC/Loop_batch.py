import os

for params in [(-2.0,0.2),(0.0,0.0)]:
    for simu in ['sncosmo','snsim']:
        cmd='python multiple_batch.py --fieldname DD --fieldid 744 --stretch '+str(params[0])+' --color '+str(params[1])+' --dirmeas ../Light_Curves_'+simu+'_noairmass --season 0 --simulator '+simu+' --dirout ../Fitted_Light_Curves_'+simu+'_noairmass'
        print cmd
        os.system(cmd)
