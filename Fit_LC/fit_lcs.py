from Fit_Single import *
import cPickle as pkl
from Telescope import *
import multiprocessing
import time
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-z", "--z", type="float", default=0.0, help="filter [%default]")
parser.add_option("-f", "--fieldname", type="string", default='WFD', help="filter [%default]")
parser.add_option("-i", "--fieldid", type="int", default=309, help="filter [%default]")
parser.add_option("-x", "--stretch", type="float", default=2.0, help="filter [%default]")
parser.add_option("-c", "--color", type="float", default=-0.2, help="filter [%default]")
parser.add_option("-d", "--dirmeas", type="string", default="None", help="filter [%default]")
#parser.add_option("--numfile", type="int", default=0, help="filter [%default]")
parser.add_option("-s", "--season", type="int", default=-1, help="filter [%default]")
parser.add_option("--T0min", type="int", default=0, help="filter [%default]")
parser.add_option("--T0max", type="int", default=10, help="filter [%default]")
parser.add_option("-t", "--sntype", type="string", default='Ia', help="filter [%default]")
parser.add_option("--dirout", type="string", default='Fitted_Light_Curves', help="filter [%default]")

opts, args = parser.parse_args()

z=opts.z
fieldname=opts.fieldname
fieldid=opts.fieldid
season=opts.season
stretch=opts.stretch
color=opts.color
dirmeas=opts.dirmeas
#numfile=opts.numfile
T0min=opts.T0min
T0max=opts.T0max
sntype=opts.sntype
dir_in='../'+dirmeas+'/'+fieldname+'/'+str(fieldid)+'/Season_'+str(season)
filename=fieldname+'_'+str(fieldid)+'_'+str(z)+'_X1_'+str(stretch)+'_C_'+str(color)
filename+='_'+str(T0min)+'_'+str(T0max)+'.pkl'
thefile=dir_in+'/'+filename
dirout=opts.dirout

dir_out=dir_in.replace(dirmeas,dirout)

print dir_out

if not os.path.isdir(dir_out) :
    os.makedirs(dir_out)

name_out=dir_out+'/'+filename

list_lc=pkl.load(open(thefile,'rb'))

telescope=Telescope(atmos=True,airmass=1.2)

n_multi=10
n_batch=len(list_lc)/n_multi

#n_batch=10
time_begin=time.time()

print len(list_lc)

fit_stack=None
#for lc in list_lc:
for i in range(n_batch):
    result_queue = multiprocessing.Queue()
    if (i%10) == 0:
        print 'Fitting',i
    time_begin=time.time()
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
        if fit_stack is None:
            #print 'hello there',resultdict[j].dtype
            fit_stack=resultdict[j]
        else:
            fit_stack=vstack([fit_stack,resultdict[j]])
    #break
    print 'total elapse time',time.time()-time_begin
"""
for fitval in fit_stack:
    print fitval['salt2.X1'],fitval['X1'],fitval['salt2.Color'],fitval['Color'],fitval['salt2.T0'],fitval['DayMax']
"""   

pkl_file = open(name_out,'wb')
pkl.dump(fit_stack, pkl_file)
pkl_file.close() 

print 'total elapse time',time.time()-time_begin
