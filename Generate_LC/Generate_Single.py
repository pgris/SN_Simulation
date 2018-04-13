from Generate_LC import *
import multiprocessing
from astropy.table import vstack,Table
from Fit_LC import *
import cPickle as pkl


class Generate_Single_LC:
    def __init__(self,z,T0,X1,X1_weight,Color,Color_weight,obs,telescope,inum,min_rf_phase,max_rf_phase,duration,date_obs,multiproc,output_q=None):
        
        self.multiproc=multiproc

        params=dict(zip(['z','DayMax','X1','Color'],[z,T0,X1,Color]))
        
        blue_cutoff=300.
        red_cutoff=800.

        """
        print obs.dtype
        for ban in 'ugrizy':
            print ban,telescope.throughputs.mean_wavelength[ban]
        """
        mean_restframe_wavelength = np.asarray([telescope.throughputs.mean_wavelength[obser['band'][-1]]/ (1. + z) for obser in obs])

        """
        for band in 'ugrizy':
            print(band,telescope.throughputs.mean_wavelength[band])
        print(mean_restframe_wavelength,[telescope.throughputs.mean_wavelength[obser['band'][-1]] for obser in obs])
        """
        #print (mean_restframe_wavelength>blue_cutoff) & (mean_restframe_wavelength<red_cutoff)
        p=(obs['mjd']-T0)/(1.+z)
        idx = (p >= min_rf_phase)&(p<=max_rf_phase)&(mean_restframe_wavelength>blue_cutoff) & (mean_restframe_wavelength<red_cutoff)
        #idx = (p >= min_rf_phase)&(p<=max_rf_phase)
        #print 'before',len(obs)
        obs=obs[idx]
        #print 'after',len(obs),np.sort(obs,order='band')
        
        self.lc=Table()

        if len(obs) >=5:
            #telescope.throughputs.Load_Atmosphere(np.median(obs['airmass']))
            #telescope.throughputs.Load_Atmosphere(1.)
            mysn=Generate_LC(params,telescope=telescope,airmass=np.median(obs['airmass']))
            #self.lc=mysn(obs['mjd'],obs['airmass'],obs['m5sigmadepth'],obs['band'],obs['exptime'])
            self.Gen_LC(mysn,obs)
            self.lc.meta=dict(zip(['X0','z','DayMax','X1','X1_weight','Color','Color_weight','mbsim','duration','min_rf_phase','max_rf_phase'],[mysn.X0,z,T0,X1,X1_weight,Color,Color_weight,mysn.mbsim,duration,min_rf_phase,max_rf_phase]))
        else:
            self.lc.meta=dict(zip(['X0','z','DayMax','X1','X1_weight','Color','Color_weight','mbsim','duration','min_rf_phase','max_rf_phase'],[-1.,z,T0,X1,X1_weight,Color,Color_weight,-1.,duration,min_rf_phase,max_rf_phase]))  

        

        if output_q is not None:
            output_q.put({inum : self.lc})

    def get_lc(self):
        return self.lc

    def Gen_LC(self,mysn,obs):

        if self.multiproc == 'yes':
            result_queue = multiprocessing.Queue()
            process=[]

            bands='grizy'

            for b in bands:
                idx= obs['band']=='LSSTPG::'+b
                sel=obs[idx]
            #print 'aiaiai',b,np.median(sel['m5sigmadepth'])
                p=multiprocessing.Process(name='Subprocess-'+b,target=mysn,args=(obs,result_queue))
                process.append(p)
                p.start()
                
            resultdict = {}
            for b in bands:
                resultdict.update(result_queue.get())

            for p in multiprocessing.active_children():
                p.join()
    
            self.lc=None
    
            for b in bands:
                if resultdict[b] is not None:
                    if self.lc is None:
                        self.lc=resultdict[b]
                    else:
                        self.lc=vstack([self.lc,resultdict[b]])
        #self.outdict['observations']=self.tot_obs
    
        #self.tot_obs=self.mysn(self.obs['mjd'],self.obs['airmass'],self.obs['m5sigmadepth'],self.obs['band'],self.obs['exptime'])
        #return mysn(obs['mjd'],obs['airmass'],obs['m5sigmadepth'],obs['band'],obs['exptime'])
       
        else:
            self.lc=mysn(obs)  
