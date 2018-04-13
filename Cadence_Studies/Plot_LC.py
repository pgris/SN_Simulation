import numpy as np
import cPickle as pkl
import matplotlib.pyplot as plt
from optparse import OptionParser

parser = OptionParser()

parser.add_option("--fieldid", type="int", default=101, help="filter [%default]")
parser.add_option("--cadence", type="int", default=6, help="filter [%default]")
parser.add_option("--X1", type="float", default=0.0, help="filter [%default]")
parser.add_option("--Color", type="float", default=0.0, help="filter [%default]")

opts, args = parser.parse_args()

fieldid=opts.fieldid
cadence=opts.cadence
X1=opts.X1
Color=opts.Color


fname='Light_Curves/Simul_'+str(fieldid)+'_Cadence_'+str(cadence)+'_X1_'+str(X1)+'_Color_'+str(Color)+'.pkl'
objs = []

f = open(fname, "r")
while 1:
    try:
        #print 'trying'
        objs.append(pkl.load(f))
    except EOFError:
        break
print(len(objs))

band='r'
for obj in objs:
    for lc in obj:
        phase=(lc['time']-lc.meta['DayMax'])/(1.+lc.meta['z'])
        idxa=lc['snr_m5']>= 5.0
        phase_5=(lc[idxa]['time']-lc.meta['DayMax'])/(1.+lc.meta['z'])
        print(lc.meta['DayMax'],lc.meta['z'],len(lc),np.min(phase),np.max(phase),np.min(phase_5),np.max(phase_5))
        idx= lc['band']=='LSST::'+band
        sel=lc[idx]
        plt.errorbar((sel['time']-sel.meta['DayMax'])/(1.+sel.meta['z']),sel['flux'],yerr=sel['fluxerr'],ls='None',color='k')
        plt.show()

