import numpy as np
import glob
from Observations import *

fieldname='WFD'
thedir='../../Ana_Cadence/OpSimLogs/'+fieldname
files=glob.glob(thedir+'/*.txt')

#print files

for fi in files:
    fieldid=int(fi.split('_')[-1].split('.')[0])
    myobs=Observations(fieldid=fieldid, filename=fi)
    print fieldid,myobs.seasons[0]['Ra'][0],myobs.seasons[0]['Dec'][0]
