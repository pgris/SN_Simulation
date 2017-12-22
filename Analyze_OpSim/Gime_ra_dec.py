import numpy as np
import glob
from Observations import *

fieldname='GalacticPlane'
thedir='OpSimLogs/'+fieldname
files=glob.glob(thedir+'/*.txt')

#print files

for fi in files:
    if fi.find('orig') == -1 and fi.find('ref') == -1:
        #print fi
        fieldid=int(fi.split('_')[-1].split('.')[0])
        myobs=Observations(fieldid=fieldid, filename=fi)
        print fieldid,myobs.seasons[0]['Ra'][0],myobs.seasons[0]['Dec'][0]
