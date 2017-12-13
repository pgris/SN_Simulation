from SN_Utils import *
import cPickle as pkl
from Telescope import *
from Generate_LC import *
import numpy as np
from optparse import OptionParser

parser = OptionParser()
parser.add_option("--type", type="string", default='low_z', help="filter [%default]")
parser.add_option("--idmin", type="int", default=0, help="filter [%default]")
parser.add_option("--idmax", type="int", default=10, help="filter [%default]")
parser.add_option("--save_file", type="string", default="no", help="filter [%default]")

opts, args = parser.parse_args()

atype=opts.type
idmin=opts.idmin
idmax=opts.idmax
save_file=opts.save_file

snutils=SN_Utils()

"""
pkl_file = open('../Generate_LC/Map_X1_C.pkl','rb')
tab_X1_c=pkl.load(pkl_file)

#print tab_X1_c

r=[]
idx= tab_X1_c['type'] == atype

for z in np.arange(0.01,1.1,0.01):
    for val in tab_X1_c[idx]:
        r.append((z,val['X1'],val['Color']))
"""

pkl_file = open('../../SNSim_Simulation/Files_Params_'+atype+'.pkl','rb')
all_params=pkl.load(pkl_file)

#print len(r)

telescope=Telescope(atmos=True,aerosol=False,airmass=1.2)

resu=[]
for val in all_params[idmin:idmax]:
    """
    params=dict(zip(['z','DayMax','X1','Color'],[val['z'],0.,val['X1'],val['Color']]))

    mysn=Generate_LC(params,telescope=telescope,airmass=1.2,X0=val['X0'],dL=val['dL'])

    print mysn.X0,mysn.dL
    """
    params=dict(zip(['x0','x1','c'],[val['X0'],val['X1'],val['Color']]))
    der_mb=snutils.Deriv_mb(params,['x0','x1','c'])

    #print der_mb
    resu.append((val['z'],val['X0'],val['X1'],val['Color'],np.asscalar(der_mb['x0']),np.asscalar(der_mb['x1']),np.asscalar(der_mb['c'])))
    #break

recres= np.rec.fromrecords(resu,names=['z','X0','X1','Color','dmb_dx0','dmb_dx1','dmb_dc'])
#names=['z','X0','X1','Color','dmb_dx0','dmb_dx1','dmb_dc']
#recres= np.array(resu)#,dtype=[(name, float) for name in names])

#print recres.shape
#print recres

if save_file == 'yes':
    name_out='Files_Deriv_mb/res_'+atype.strip()+'_'+str(idmin)+'_'+str(idmax)+'.pkl'
    pkl_file = open(name_out,'wb')
    pkl.dump(recres, pkl_file)
    pkl_file.close() 
