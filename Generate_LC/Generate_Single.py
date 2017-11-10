from Generate_LC import *
import multiprocessing
from astropy.table import vstack,Table
from Fit_LC import *
import cPickle as pkl


class Generate_Single_LC:
    def __init__(self,z,T0,X1,X1_weight,Color,Color_weight,obs,telescope,inum,min_rf_phase,max_rf_phase,duration,date_obs,output_q):
        
        params={}
        params['z']=z
        params['DayMax']=T0
        params['X1']=X1
        params['Color']=Color

        
        self.telescope=telescope
        self.z=z
        self.T0=T0
        self.X1=X1
        self.X1_weight=X1_weight
        self.Color=Color
        self.Color_weight=Color_weight
        """
        self.outdict={}
        self.outdict=params
        self.outdict['status']='unkown'
        self.outdict['fit']=None
        self.outdict['mbsim']=-999.
        self.outdict['observations']=None
        """
        self.duration=duration
        self.date_obs=date_obs
        self.min_rf_phase=min_rf_phase
        self.max_rf_phase=max_rf_phase
        self.bands_rest = 'grizy'
        self.bands_rest = 'g'

        p=(obs['mjd']-T0)/(1.+z)
        idx = (p >= min_rf_phase)&(p<=max_rf_phase)
        self.obs=obs[idx]

        if len(self.obs) >=5:
            self.mysn=Generate_LC(params,telescope=self.telescope)
            self.Gen_LC()

        #output_q.put({inum : (self.mysn.metadata,self.tot_obs)})
       
        self.tot_obs.meta['X0']=self.mysn.X0
        self.tot_obs.meta['z']=z
        self.tot_obs.meta['DayMax']=T0
        self.tot_obs.meta['X1']=X1
        self.tot_obs.meta['X1_weight']=X1_weight
        self.tot_obs.meta['Color']=Color
        self.tot_obs.meta['Color_weight']=Color_weight
        self.tot_obs.meta['mbsim']=self.mysn.mbsim
        self.tot_obs.meta['duration']=self.duration
        self.tot_obs.meta['min_rf_phase']=min_rf_phase
        self.tot_obs.meta['max_rf_phase']=max_rf_phase
        
        output_q.put({inum : self.tot_obs})


    def Gen_LC(self):

        """
        result_queue = multiprocessing.Queue()
        process=[]

        for b in self.bands_rest:
            idx= self.obs['band']=='LSSTPG::'+b
            sel=self.obs[idx]
            print 'aiaiai',b,np.median(sel['m5sigmadepth'])
            p=multiprocessing.Process(name='Subprocess-'+b,target=self.mysn,args=(sel['mjd'],sel['airmass'],sel['m5sigmadepth'],b,sel['exptime'],result_queue))
            process.append(p)
            p.start()

        resultdict = {}
        for b in self.bands_rest:
            resultdict.update(result_queue.get())

        for p in multiprocessing.active_children():
            p.join()
    
        self.tot_obs=None
    
        for b in self.bands_rest:
            if resultdict[b] is not None:
                if self.tot_obs is None:
                    self.tot_obs=resultdict[b]
                else:
                    self.tot_obs=vstack([self.tot_obs,resultdict[b]])
        #self.outdict['observations']=self.tot_obs
        """
        self.tot_obs=self.mysn(self.obs['mjd'],self.obs['airmass'],self.obs['m5sigmadepth'],self.obs['band'],self.obs['exptime'])

       
