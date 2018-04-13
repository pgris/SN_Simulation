import os
import numpy as np

fieldname='DD'

thefile='fieldIDs_minion_1016_'+fieldname+'.txt'

fields=np.loadtxt(thefile,dtype={'names': ('name','fieldid'),'formats': ('S8','i4')})

print fields['fieldid']

for fieldid in fields['fieldid']:
    cmd = 'python Dump_OpSim_in_File.py --fieldname '+fieldname+' --fieldid '+str(fieldid)
    os.system(cmd)
    
