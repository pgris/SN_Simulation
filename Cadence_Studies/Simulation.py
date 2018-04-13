from Observations import *
from Generate_LC import Generate_LC
from Telescope import *
import cPickle as pkl
import multiprocessing

blue_cutoff=300.
red_cutoff=800.
min_rf_phase=-20.
max_rf_phase=60.

class Simulation:
    def __init__(self,fieldid,cadence,X1,Color,telescope):
        
        X1_Color=[(X1,Color)]
        zrange=[0.01]
        zrange+=[val for val in np.arange(0.1,1.5,0.1)]
        DayMax_step=0.5

        filename='Observations/Obs_'+str(fieldid)+'_'+str(cadence)+'.txt'
        obs_tot=Observations(fieldid=fieldid, filename=filename)
        season=0
        
        obs=obs_tot.seasons[season]

        r=[]
        mjds=np.arange(min(obs['mjd']),max(obs['mjd']),DayMax_step)

        for z in zrange:
            for mjd in mjds:
                for (x1,c) in X1_Color:
                    r.append((z,x1,c,mjd))

        params=np.rec.array(r, dtype=[('z', 'f8'),('X1', 'f8'), ('Color', 'f8'),('DayMax','f8')])
        
        name='Light_Curves/Simul_'+str(fieldid)+'_Cadence_'+str(cadence)+'_X1_'+str(X1_Color[0][0])+'_Color_'+str(X1_Color[0][1])+'.pkl'

        #simulate
        self.Simul_LCs(params,obs,telescope,fname=name)

    def Simul_LCs(self,params,obs,telescope,fname=''):

        fout = open(fname, "w")
        print('Number of LC to simulate:',len(params))
        lc_list=[]


        n_per_batch=10
        
        inter=range(0,len(params),n_per_batch)
        inter=np.append(inter,len(params))

        for jo in range(len(inter)-1):
            ida=inter[jo]
            idb=inter[jo+1]
            result_queue = multiprocessing.Queue()
        
            for j in range(ida,idb):
                p=multiprocessing.Process(name='Subprocess-'+str(j),target=self.Simul_LC,args=(params[j],obs,telescope,j,result_queue))
                p.start()
 
            resultdict = {}

            for j in range(ida,idb):
                resultdict.update(result_queue.get())

            for p in multiprocessing.active_children():
                p.join()

            for j in range(ida,idb):
                lc_list.append(resultdict[j])
            
            if len(lc_list) >= 100:
                print('Dumping:',len(lc_list))
                pkl.dump(lc_list,fout)
                lc_list=[]

        if lc_list:
            pkl.dump(lc_list,fout)
        fout.close()

    def Simul_LC(self,param,obs,telescope,j,output_q=None):
   
        lc_sncosmo=Generate_LC(param,telescope=telescope)
        z=param['z']
        T0=param['DayMax']
        mean_restframe_wavelength = np.asarray([telescope.throughputs.mean_wavelength[obser['band'][-1]]/ (1. + z) for obser in obs])
    
        p=(obs['mjd']-T0)/(1.+z)
        idx = (p >= min_rf_phase)&(p<=max_rf_phase)&(mean_restframe_wavelength>blue_cutoff) & (mean_restframe_wavelength<red_cutoff)
        obs=obs[idx]
    
        lc=lc_sncosmo(obs)
        lc.meta=dict(zip(['z','DayMax','X1','Color','min_rf_phase','max_rf_phase'],[z,T0,param['X1'],param['Color'],min_rf_phase,max_rf_phase]))

        if output_q is not None:
            output_q.put({j : lc})
        else:
            return lc
