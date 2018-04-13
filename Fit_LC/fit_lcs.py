from Fit_Single import *
import cPickle as pkl
from Telescope import *
import multiprocessing
import time
from optparse import OptionParser
#from SN_Utils import *
import sncosmo
import h5py
from astropy.table import Table,vstack

def Dump_in_File(name,fit_lc):

    val_tot=None

    for key, val in fit_lc.items():
        if val_tot is None:
            val_tot=val
        else:
            val_tot=vstack([val_tot,val])

    print('dumping in ',name)
    if not os.path.isfile(name):
        val_tot.write(name, path='fit_1', compression=True)
    else:
        f = h5py.File(name,'r')
        key_num=len(f.keys())
        print 'hello',key_num
        f.close()
        val_tot.write(name,path='fit_'+str(key_num+1), append=True,compression=True)
    

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
parser.add_option("--multiproc", type="string", default='yes', help="filter [%default]")
parser.add_option("--DayMax", type="float", default=-1., help="filter [%default]")

opts, args = parser.parse_args()

z=opts.z
fieldname=opts.fieldname
fieldid=opts.fieldid
season=opts.season
stretch=opts.stretch
color=opts.color
dirmeas=opts.dirmeas
DayMax=opts.DayMax
#numfile=opts.numfile
#T0min=opts.T0min
#T0max=opts.T0max
sntype=opts.sntype
main_dir_in='/sps/lsst/data/dev/pgris/'
dir_in=main_dir_in+dirmeas+'/'+fieldname+'/'+str(fieldid)+'/Season_'+str(season)+'/z_'+str(z)
filename=fieldname+'_'+str(fieldid)+'_'+str(z)+'_X1_'+str(stretch)+'_C_'+str(color)
if DayMax > 0.:
    filename+='_DayMax_'+str(DayMax)
#filename+='_'+str(T0min)+'_'+str(T0max)+'.hdf5'
filename+='.hdf5'
thefile=dir_in+'/'+filename
dirout=opts.dirout
multiproc=opts.multiproc
DayMax=opts.DayMax

dir_out=dir_in.replace(dirmeas,dirout)

print 'Processing',filename

"""
if z < 1:
    deriv_file='Files_Deriv_mb/Derib_mb_low_z_all.pkl'
else:
    deriv_file='Files_Deriv_mb/Derib_mb_high_z_all.pkl'

deriv_mb=pkl.load(open(deriv_file,'rb'))
"""
print dir_out

if not os.path.isdir(dir_out) :
    os.makedirs(dir_out)

name_out=dir_out+'/'+filename

if os.path.isfile(name_out):
    print 'file',name_out,'exists'
    os.remove(name_out)

print 'loading',thefile

#list_lc=pkl.load(open(thefile,'rb'))

f = h5py.File(thefile,'r')

#print f.keys()

telescope=Telescope(atmos=True,airmass=1.2)
transmission=telescope.throughputs

for filtre in 'ugrizy':
    if telescope.airmass > 0:
        band=sncosmo.Bandpass(transmission.atmosphere[filtre].wavelen,transmission.atmosphere[filtre].sb, name='LSST::'+filtre,wave_unit=u.nm)
    else:
        band=sncosmo.Bandpass(transmission.system[filtre].wavelen,transmission.system[filtre].sb, name='LSST::'+filtre,wave_unit=u.nm) 
    sncosmo.registry.register(band, force=True)

nlc=len(f.keys())

n_multi=8
#n_multi=1
n_batch=nlc/n_multi

#n_batch=10
time_begin=time.time()

print nlc

fit_stack=None
#for lc in list_lc:

fitted_lc={}
#snutils=SN_Utils()

if multiproc == 'yes':
    num_tot=0
    
    batch=np.arange(0,nlc,n_multi)
    batch=np.append(batch,nlc)

    print(batch,nlc%n_multi)

    for i in range(len(batch)-1):
        result_queue = multiprocessing.Queue()
        if (i%10) == 0:
            print('Fitting',i,n_multi)
           
        ida=batch[i]
        idb=batch[i+1]

        for j in range(ida,idb):
            num_tot+=1
            #lc=list_lc[num]
            lc=Table.read(thefile, path='lc_'+str(j))
            p=multiprocessing.Process(name='Subprocess-'+str(j),target=Fit_Single_LC,args=(lc,telescope,j,result_queue))
            
            p.start()

        
        resultdict = {}
        for j in range(ida,idb):
            resultdict.update(result_queue.get())
        

        for p in multiprocessing.active_children():
            p.join()

        
        for j in range(ida,idb):
            fitted_lc[j]=resultdict[j]

        if i> 0 and (i%80)==0:
            print 'Dumping',num_tot,i
            Dump_in_File(name_out,fitted_lc)
            fitted_lc={}

            """
            if not os.path.isfile(name_out):
                resultdict[j].write(name_out, path='fit_'+str(num), compression=True)
            else:
                resultdict[j].write(name_out,path='fit_'+str(num), append=True,compression=True)
            """
        
else:
    
    #for i,lc in enumerate(list_lc):

    fitted_lc={}
    for i,key in enumerate(f.keys()):
        lc=Table.read(thefile, path=key)
        
        print('Fitting',i,len(lc))
        #time_begin=time.time()
        fitted_lc[key]=Fit_Single_LC(lc,telescope,-1).Summary()
        #print 'total elapse time fit',time.time()-time_begin
        """
        if fit_stack is None:
            #print 'hello there',resultdict[j].dtype
            fit_stack=fitted_lc
        else:
            fit_stack=vstack([fit_stack,fitted_lc])
        """
        """
        if not os.path.isfile(name_out):
             fitted_lc.write(name_out, path=key.replace('lc','fit'), compression=True)
        else:
            fitted_lc.write(name_out,path=key.replace('lc','fit'), append=True,compression=True)
        
        if i >=50:
            break
        """
    
        #if (i%500)==0:
         #   Dump_in_File(name_out,fitted_lc)

    #break
#print 'total elapse time',time.time()-time_begin
"""
for fitval in fit_stack:
    print fitval['salt2.X1'],fitval['X1'],fitval['salt2.Color'],fitval['Color'],fitval['salt2.T0'],fitval['DayMax']
"""   

"""
pkl_file = open(name_out,'wb')
pkl.dump(fit_stack, pkl_file)
pkl_file.close() 
"""
Dump_in_File(name_out,fitted_lc)
print 'total elapse time',time.time()-time_begin
