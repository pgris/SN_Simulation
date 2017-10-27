from Fit_Single import *
import cPickle as pkl
from Telescope import *
import multiprocessing
import time

thefile='../Light_Curves/DD/290/Season_0/DD_290_0.1_0.3_X1_-999.0_C_-999.0_0.pkl'

list_lc=pkl.load(open(thefile,'rb'))

telescope=Telescope(atmos=True,airmass=1.2)

n_multi=2
n_batch=len(list_lc)/n_multi

time_begin=time.time()

print len(list_lc)

fit_list=[]
#for lc in list_lc:
for i in range(n_batch):
    result_queue = multiprocessing.Queue()
    for j in range(0,n_multi):
        num=j+n_multi*i
        p=multiprocessing.Process(name='Subprocess-'+str(j),target=Fit_Single_LC,args=(list_lc[num],telescope,j,result_queue))
            
        p.start()

    resultdict = {}
    for j in range(0,n_multi):
        resultdict.update(result_queue.get())

    for p in multiprocessing.active_children():
        p.join()

    for j in range(0,n_multi):
        #print 'hello there',resultdict[j].dtype
        fit_list.append(resultdict[j])
    break

for fitval in fit_list:
    print fitval['salt2.X1'],fitval['salt2.Color']

print 'total elapse time',time.time()-time_begin
