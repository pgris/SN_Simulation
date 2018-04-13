import numpy as np
import h5py
from astropy.table import Table,vstack,Column
from SN_Utils import SN_Utils
import time
import multiprocessing
from optparse import OptionParser
import os

def Calc_cov_mb(vals,j,output_q):

    #print(vals['status'],vals['fit_status'])
    #print(vals.dtype,len(vals))
    if vals['fit_status']=='fitok':
        cov=np.ndarray(shape=(3,3), dtype=float, order='F')
        cov[0,0]=vals['salt2.CovX0X0']
        cov[1,1]=vals['salt2.CovX1X1']
        cov[2,2]=vals['salt2.CovColorColor']
        cov[0,1]=vals['salt2.CovX1X0']
        cov[0,2]=vals['salt2.CovX0Color']
        cov[1,2]=vals['salt2.CovX1Color']
        cov[2,1]=cov[1,2]
        cov[1,0]=cov[0,1]
        cov[2,0]=cov[0,2]
    #print(cov) 
        #print('ooooo',len(vals))
        params=dict(zip(param_names,[vals['salt2.X0'], vals['salt2.X1'],vals['salt2.Color']] ))
        
        covar_mb=snutils.Covar(params,cov,param_names)
        #print('params',params,covar_mb)
    else:
        covar_mb = None
    """
    for key, res in covar_mb.items():
        outdict[key].append(res)
        print('hello',key,len(outdict[key]),outdict[key])
    """
    if output_q is not None:
        output_q.put({j:covar_mb})
    else:
        return covar_mb


def Multiproc(tot_data,nbatch,outdict,outvars):
    nbatch=8
    batch=np.arange(0,len(tot_data),nbatch)
    batch=np.append(batch,len(tot_data))

    #print('oulala',len(tot_data),tot_data)
    for i in range(len(batch)-1):
        result_queue = multiprocessing.Queue()
        if i > 0 and (i%10) == 0:
            print 'Processing',i

        ida=batch[i]
        idb=batch[i+1]

        for j in range(ida,idb):
            #print('pppoo',j,tot_data[j])
            p=multiprocessing.Process(name='Subprocess-'+str(j),target=Calc_cov_mb,args=(tot_data[j],j,result_queue))
            p.start()
    
        resultdict = {}
        for j in range(ida,idb):
            resultdict.update(result_queue.get())
        

        for p in multiprocessing.active_children():
            p.join()

        for j in range(ida,idb):
            if resultdict[j] is not None:
                for key, val in resultdict[j].items():
                    outdict[key].append(val)
            else:
                for key in outvars:
                    outdict[key].append(-999.)
                #print(key,len(outdict[key]),outdict[key])

parser = OptionParser()
parser.add_option("-z", "--z", type="float", default=0.0, help="filter [%default]")
parser.add_option("-f", "--fieldname", type="string", default='WFD', help="filter [%default]")
parser.add_option("-i", "--fieldid", type="int", default=309, help="filter [%default]")
parser.add_option("-x", "--stretch", type="float", default=2.0, help="filter [%default]")
parser.add_option("-c", "--color", type="float", default=-0.2, help="filter [%default]")
parser.add_option("-d", "--dirmeas", type="string", default="None", help="filter [%default]")
#parser.add_option("--numfile", type="int", default=0, help="filter [%default]")
parser.add_option("-s", "--season", type="int", default=-1, help="filter [%default]")
#parser.add_option("--T0min", type="int", default=0, help="filter [%default]")
#parser.add_option("--T0max", type="int", default=10, help="filter [%default]")
parser.add_option("-t", "--sntype", type="string", default='Ia', help="filter [%default]")
parser.add_option("--dirout", type="string", default='Fitted_Light_Curves', help="filter [%default]")
parser.add_option("--multiproc", type="string", default='no', help="filter [%default]")
parser.add_option("--DayMax", type="float", default=-1., help="filter [%default]")

opts, args = parser.parse_args()

fieldname=opts.fieldname
fieldid=opts.fieldid
season=opts.season
X1=opts.stretch
Color=opts.color
dirmeas=opts.dirmeas
dirout=opts.dirout
z=opts.z
multiproc=opts.multiproc
DayMax=opts.DayMax

print('DayMax?',DayMax)
dirmain='/sps/lsst/data/dev/pgris'
fi=dirmain+'/'+dirmeas
dirsec=fieldname+'/'+str(fieldid)+'/Season_'+str(season)+'/z_'+str(z)
dirout+='/'+dirsec
fname=fieldname+'_'+str(fieldid)+'_'+str(z)+'_X1_'+str(X1)+'_C_'+str(Color)
if DayMax > 0:
    fname+='_DayMax_'+str(DayMax)
fname+='.hdf5'
fi+='/'+dirsec+'/'+fname

outfile=dirmain+'/'+dirout+'/'+fname
if not os.path.isdir(dirmain+'/'+dirout):
    os.makedirs(dirmain+'/'+dirout)

tot_data=None
f = h5py.File(fi,'r')
print('try loading',fi,len(f.keys()))
n_lc=0
for i,keyb in enumerate(f.keys()):
    tab=Table.read(fi, path=keyb)
    #print(i,tab)
    n_lc+=len(tab)
    if tot_data is None:
        tot_data=tab
    else:
        tot_data=vstack([tot_data,tab])
print('number of lc',n_lc)

#print(tot_data.colnames)
"""
for name in tot_data.colnames:
    if 'Cov' in name:
        print name
"""
param_names=['x0','x1','c']
outvars=['mb_recalc','salt2.CovColormb','salt2.CovX1mb','salt2.CovX0mb','salt2.Covmbmb']



snutils=SN_Utils()
time_begin=time.time()

outdict={}
for var in outvars:
        outdict[var]=[]
if multiproc == 'yes':
    print('multiproc')
    Multiproc(tot_data,8,outdict,outvars)
else:
    print('No multiproc')
    res=Calc_cov_mb(tot_data,0,None)
    if res is None:
        for key in outvars:
            outdict[key].append(-999.)
    else:
        for key,val in res.items():
            outdict[key].append(val)

#print('hello',outdict)
for var in outvars:
    tot_data.add_column(Column(outdict[var],name=var))

tot_data.write(outfile, path='fit_1', compression=True)
print('time',tot_data['mbfit'],tot_data['mb_recalc'],time.time()-time_begin)
