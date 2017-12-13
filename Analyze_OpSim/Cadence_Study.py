import numpy as np
from Observations import *
import pylab as plt


def Cadence(obs):

    oob=obs['mjd'][1:]-obs['mjd'][:-1]

    return np.median(oob)

thefile='fieldIDs_minion_1016_WFD.txt'

fieldids=np.loadtxt(thefile,dtype={'names': ('name','fieldid'),'formats': ('S8','i4')})

myobs={}

season=0
fieldname='WFD'
thedir='OpSimLogs/'+fieldname

conv=dict(zip(['LSSTPG::'+band for band in 'ugrizy'],[1,2,3,4,5,6]))
class_map = lambda x: conv[x]
for i,fieldid in enumerate(fieldids):
    name='Observations_'+fieldname+'_'+str(fieldid[1])+'.txt'
    obs=Observations(fieldid=fieldid[1], filename=thedir+'/'+name)
    myobs[fieldid[1]]=obs.seasons[season]

    if i>100: 
        break

for fieldid in myobs.keys():
    print fieldid,Cadence(myobs[fieldid])
    #plt.plot(myobs[fieldid]['mjd'],[fieldid]*len(myobs[fieldid]['mjd']))
    plt.plot(myobs[fieldid]['mjd'],map(class_map,myobs[fieldid]['band']))
   
plt.show()
