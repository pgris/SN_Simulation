from Telescope import *
import cPickle as pkl
import matplotlib.pyplot as plt
import time
import operator
from optparse import OptionParser
from Simulation import *
from Fit_Ana import *


blue_cutoff=300.
red_cutoff=800.
min_rf_phase=-20.
max_rf_phase=60.

def Simul_LCs(params,obs,telescope,fname=''):

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
            p=multiprocessing.Process(name='Subprocess-'+str(j),target=Simul_LC,args=(params[j],obs,telescope,j,result_queue))
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

def Simul_LC(param,obs,telescope,j,output_q=None):
   
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

def Ana_LCs(fname,telescope):

    time_ref=time.time()

    transmission=telescope.throughputs
    for filtre in 'ugrizy':
        if telescope.airmass > 0:
            band=sncosmo.Bandpass(transmission.atmosphere[filtre].wavelen,transmission.atmosphere[filtre].sb, name='LSST::'+filtre,wave_unit=u.nm)
        else:
            band=sncosmo.Bandpass(transmission.system[filtre].wavelen,transmission.system[filtre].sb, name='LSST::'+filtre,wave_unit=u.nm) 
        sncosmo.registry.register(band, force=True)

    f = open(fname, "r")
    objs = []
    r=[]
    n_per_batch=10
    fit_lc=Fit_LC(telescope=telescope)

    while 1:
        try:
        #print 'trying'
            objs.append(pkl.load(f))
        except EOFError:
            break
    print(len(objs))
    for obj in objs:
        #print(len(obj),type(obj))
        inter=range(0,len(obj),n_per_batch)
        inter=np.append(inter,len(obj))
        
        for jo in range(len(inter)-1):
            ida=inter[jo]
            idb=inter[jo+1]
            result_queue = multiprocessing.Queue()
        
            for j in range(ida,idb):
                p=multiprocessing.Process(name='Subprocess-'+str(j),target=Ana_LC,args=(obj[j],fit_lc,j,result_queue))
                p.start()
 
            resultdict = {}

            for j in range(ida,idb):
                resultdict.update(result_queue.get())

            for p in multiprocessing.active_children():
                p.join()

            for j in range(ida,idb):
                r+=resultdict[j][0]
                names=resultdict[j][1]
        if len(r)%100 == 0:
            print('there',len(r),time.time()-time_ref)
                #print('result',names)
            #break
    """
    names=['DayMax','z']
    for val in ['all','g','r','i','z','y']:
        names+=['snr_'+val,'snr_5_'+val,'Nbef_'+val,'Naft_'+val,'Nbef_5_'+val,'Naft_5_'+val,'phase_min_'+val,'phase_max_'+val,'phase_min_5_'+val,'phase_max_5_'+val]
    #lc_ana=np.rec.fromrecords(r,names=['DayMax','z','band','snr','snr_5','Nbef','Naft','Nbef_5','Naft_5','phase_min','phase_max','phase_min_5','phase_max_5'])
    #print('yes',len(r),len(names))
    """
    lc_ana=np.rec.fromrecords(r,names=names)
    #print(lc_ana)
    fout = open('Ana_'+fname, "w")
    pkl.dump(lc_ana,fout)
    fout.close()
    #Plot(lc_ana)

def Ana_LC(tab,fit_lc,j,output_q=None):

    #print(tab)

   

    idxc=tab['flux']/tab['fluxerr']>5.
     #fit this LC
    res,fitted_model,mbfit,fit_status=fit_lc(tab[idxc])

    ro, names=Fill_Infos(res,fitted_model,mbfit,fit_status)
    #print(res,fitted_model,mbfit,fit_status)

    nb,na=Bef_Aft(tab,tab.meta['DayMax'])
    nb5,na5=Bef_Aft(tab[idxc],tab.meta['DayMax'])
    phase_min,phase_max= phase(tab)
    phase_min_5,phase_max_5= phase(tab[idxc])

    #print('before',ro,names)
    ro+=[tab.meta['DayMax'],tab.meta['z'],snr(tab),snr(tab[idxc]),nb,na,nb5,na5,phase_min,phase_max,phase_min_5,phase_max_5]
    names+=['DayMax','z','snr_all','snr_5_all','Nbef_all','Naft_all','Nbef_5_all','Naft_5_all','phase_min_all','phase_max_all','phase_min_5_all','phase_max_5_all']

    for band in 'grizy':
        #ro=[]
        idx = tab['band']=='LSST::'+band
        idxb= (tab['band']=='LSST::'+band)&(tab['flux']/tab['fluxerr']>5.)
        #print(band,snr(tab[idx]),snr(tab[idxb]),Bef_Aft(tab[idx],tab.meta['DayMax']))
        nb,na=Bef_Aft(tab[idx],tab.meta['DayMax'])
        nb5,na5=Bef_Aft(tab[idxb],tab.meta['DayMax'])
        phase_min,phase_max= phase(tab[idx])
        phase_min_5,phase_max_5= phase(tab[idxb])
        #r.append((tab.meta['DayMax'],tab.meta['z'],band,snr(tab[idx]),snr(tab[idxb]),nb,na,nb5,na5,phase_min,phase_max,phase_min_5,phase_max_5))
        ro+=[snr(tab[idx]),snr(tab[idxb]),nb,na,nb5,na5,phase_min,phase_max,phase_min_5,phase_max_5]
        names+=['snr_'+band,'snr_5_'+band,'Nbef_'+band,'Naft_'+band,'Nbef_5_'+band,'Naft_5_'+band,'phase_min_'+band,'phase_max_'+band,'phase_min_5_'+band,'phase_max_5_'+band]
    
    #r.append((tab.meta['DayMax'],tab.meta['z'],'all',snr(tab),snr(tab[idxc]),nb,na,nb5,na5,phase_min,phase_max,phase_min_5,phase_max_5))

    rsb, nsb= Var_Science_Book(tab,tab.meta['DayMax'],tab.meta['z'])

    ro+=rsb
    names+=nsb

    r=[]
    r.append(tuple(ro))

    if output_q is not None:
        output_q.put({j : (r,names)})
    else:
        return r

def Fill_Infos(res,fitted_model,mbfit,fit_status):

    r=[]
    names=[]
    vars_salt=['x1','x0','c']
    corresp=dict(zip(vars_salt,['X1','X0','Color']))
    corresp['t0']='t0'
    corresp['z']='z'
    r+=[fit_status]
    names+=['fit_status']
    if fit_status == 'fitok':
        r+=[mbfit]
        names+=['mbfit']

        corr={}
        for i,pal in enumerate(res['vparam_names']):
            corr[pal]=i
        for i,par in enumerate(fitted_model.param_names):
            r+=[fitted_model.parameters[i]]
            names+=['salt2.'+corresp[par]]
        if res['covariance'] is not None:
            for k,var in enumerate(vars_salt):
                for varb in vars_salt[k:]:
                    r+=[res['covariance'][corr[var]][corr[varb]]]
                    names+=['salt2.Cov'+corresp[var]+corresp[varb]]
            r+=[res.errors['t0']**2]
            names+=['salt2.Covt0t0']


    else: 
        r+=[-999.]
        names+=['mbfit']
        thenames=['salt2.z', 'salt2.t0', 'salt2.X0', 'salt2.X1', 'salt2.Color', 'salt2.CovX1X1', 'salt2.CovX1X0', 'salt2.CovX1Color', 'salt2.CovX0X0', 'salt2.CovX0Color', 'salt2.CovColorColor', 'salt2.Covt0t0']
        r+=[-999.]*len(thenames)
        names+=thenames
       
    return r, names


def Var_Science_Book(lc,T0,z):
       
    r=[]
    names=[]
 
    lc.sort('time')
    phases=(lc['time']-T0)/(1.+z)
        #print('phases',phases)
    idx = phases < -5.
    idxb = phases > 30.
    r+=[len(phases[idx]),len(phases[idxb])]
    names+=['N_Phase_m5','N_Phase_p30']
    

    diffa,seldiffa=Get_Diff(lc,phases,-20.,60.)
    

    r+=[len(seldiffa)-1]
    names+=['N_nights_m20_p_p30']

    diffb,seldiffb=Get_Diff(lc,phases,-5.,30.)
    nearpeak=-1
    if len(diffb) > 0:
        nearpeak=np.max(diffb)/(1.+z)
    r+=[nearpeak]
    names+=['Near_peak_gap']

    snr_refs=[10.,15.,20.]
    for band in 'grizy':
        idf = lc['band']=='LSST::'+band
        sel_idf= lc[idf]
        if len(sel_idf) > 0:
            #print('band',band,len(sel_idf))
            snr_max=np.max(sel_idf['flux']/sel_idf['fluxerr'])
            for valref in snr_refs:
                nres=0
                if snr_max > valref:
                    nres=1
                r+=[nres]
                names+=['N_'+band+'_snrmax_'+str(int(valref))]
                #self.dict_quality['N_'+band+'_snrmax_'+str(int(valref))]=nres
        
        else:
            for valref in snr_refs:
                r+=[0.]
                names+=['N_'+band+'_snrmax_'+str(int(valref))]
       

    return r,names

def Get_Diff(lc,phases,phase_min,phase_max):

    idxc=(phases > phase_min)&(phases < phase_max)
    
    sel=lc[idxc]
    
    diff_time=1.
    
    diff=[io-jo for jo,io in zip(sel['time'][:-1], sel['time'][1:])]

    seldiff=[i+1 for i in range(len(diff)) if diff[i]>= diff_time]
    seldiff=[0]+seldiff+[len(sel)]

    return diff,seldiff

def phase(tab):
    if len(tab) > 0:
        phases=(tab['time']-tab.meta['DayMax'])/(1.+tab.meta['z'])
        return np.min(phases),np.max(phases)
    else:
       return -999.,-999. 

def snr(sel):
    if len(sel) > 0:
        return np.sum(sel['flux'])/np.sqrt(np.sum(sel['fluxerr']*sel['fluxerr']))
    else:
        return 0.0

def Bef_Aft(lc, T0):
    if len(lc) > 0:
        diff=lc['time']-T0
        idxa=diff <= 0
        idxb=diff > 0
        return len(lc[idxa]),len(lc[idxb])
    else:
        return 0.0,0.0

def Plot(tab):
    bcolor = {'u':'b', 'g':'c', 'r':'g', 'i':'y', 'z':'r', 'y':'m'}

    for band in 'grizy':
        #figa, axa = plt.subplots()
        #figb, axb = plt.subplots()
        figc, axc = plt.subplots()
        figc.suptitle(band+' band')
        for z in np.unique(tab['z']):
            idx = (tab['z']==z)
            sel=tab[idx]
            #axa.plot(sel['DayMax'],sel['snr_'+band],color=bcolor[band],marker='.')
            #axa.plot(sel['DayMax'],sel['snr_5_'+band],color=bcolor[band],marker='o')

            #axb.plot(sel['Nbef_'+band],sel['snr_'+band])
            idxb = (sel['phase_min_'+band]>-90.)
            selb=sel[idxb]
            axc.plot(sel['phase_min_all'],sel['snr_'+band])

    figa, axa = plt.subplots()
    for z in np.unique(tab['z']):
    #for z in [0.3]:
        idx = np.abs(tab['z']-z)<1.e-5
        selb=tab[idx]
        #axa.plot(selb['phase_min_all'],selb['snr_all'])
        axa.plot(selb['phase_max_5_all'],selb['snr_5_all'])

    figb, axb = plt.subplots() 
    idd = tab['fit_status']=='fitok'
    seld=tab[idd]
    axb.plot(seld['z'],seld['salt2.CovColorColor'])

    plt.show()


def Presel(tab):
    idx=(tab['phase_max_all']>=10.)&(tab['phase_min_all']<=-5.)
    return tab[idx]

def Select_all(tot_data,selec):
        
     valscp=tot_data.copy()
     for ic,sel in enumerate(selec):
         valscp=Select_Cut(valscp,sel[0],sel[1],sel[2])
         
     return valscp

def Select_Cut(data,var,comp,val):

        tot_data=None
        
        for vv in var:
            if tot_data is None:
                tot_data=data[vv]
            else:
                tot_data+=data[vv]


        idx = comp(tot_data,val)

        #print idx
        return data[idx]

def Calc_Effi(tab,selection):
    
    ro=[]
    sel=Select_all(tab,selection)

    for z in np.unique(tab['z']):
    #for z in [0.1]:
        idx = np.abs(tab['z']-z)<1.e-5
        tab_z=tab[idx]
        idxb=np.abs(sel['z']-z)<1.e-5
        sel_z=sel[idxb]
        #print(sel['phase_max_all'],sel['phase_min_all'])
        #idxb=(sel['phase_max_all']>=10.)&(sel['phase_min_all']<=-5.)
        
        #print(selb,idxb)
        effi=float(len(sel_z))/float(len(tab_z))
        err=np.sqrt(1.-effi)*effi/np.sqrt(float(len(tab_z)))
        #print('hello',z,effi,err,len(sel))
        mjds_min=(1.+z)*sel['phase_min_all']+sel['DayMax']
        mjds_max=(1.+z)*sel['phase_max_all']+sel['DayMax']
        print(z,np.max(mjds_min),np.min(mjds_max))
        #print((1.+z)*sel['phase_min_all'],sel['DayMax'])
        ro.append((z,effi,err))

    return np.rec.fromrecords(ro,names=['z','effi','err_effi'])

def Plot_Effi(axa,tab):

    #figa, axa = plt.subplots()
    axa.errorbar(tab['z'],tab['effi'],yerr=tab['err_effi'])


def get_operator_fn(op):
    return {
        '+' : operator.add,
        '-' : operator.sub,
        '*' : operator.mul,
        '/' : operator.div,
        '%' : operator.mod,
        '^' : operator.xor,
        '>=': operator.ge,
        '>' : operator.gt,
        '<=': operator.le,
        '<' : operator.lt,
        '==': operator.eq,
        }[op]

def Load_Selection(selname):

    test=np.loadtxt(selname,dtype={'names': ('cut','comp','val','type'),'formats': ('S150','S2','S8','S5')})

    print(test)

    sel=[]
    for cut in test:
        thecut=[]
        for val in cut['cut'].split('+'):
            thecut.append(val)
        #print(cut['type'])
        if cut['type'] != 'str':
            sel.append((thecut,get_operator_fn(cut['comp']),eval(cut['type']+'('+cut['val']+')')))
        else:
            sel.append((thecut,get_operator_fn(cut['comp']),cut['val']))
            
    return sel


parser = OptionParser()
parser.add_option("--process_step", type="string", default='simulation', help="filter [%default]")
parser.add_option("--fieldid", type="int", default=100, help="filter [%default]")
parser.add_option("--cadence", type="int", default=3, help="filter [%default]")
parser.add_option("--X1", type="float", default=0.0, help="filter [%default]")
parser.add_option("--Color", type="float", default=0.0, help="filter [%default]")

opts, args = parser.parse_args()

process_step=opts.process_step
fieldid=opts.fieldid
cadence=opts.cadence
X1=opts.X1
Color=opts.Color

#Load instrument

atmos=True
telescope=Telescope(atmos=atmos,aerosol=False,airmass=1.2)

if process_step == 'simulation':
    """
    #supernova parameters

    X1_Color=[(X1,Color)]
    zrange=[0.01]
    zrange+=[val for val in np.arange(0.1,1.5,0.1)]
    DayMax_step=0.5

    print zrange
    #load observations


    filename='Observations/Obs_'+str(fieldid)+'_'+str(cadence)+'.txt'
    obs_tot=Observations(fieldid=fieldid, filename=filename)
    season=0

    obs=obs_tot.seasons[season]

    #Load parameters for simulation

    r=[]
    mjds=np.arange(min(obs['mjd']),max(obs['mjd']),DayMax_step)

    for z in zrange:
        for mjd in mjds:
            for (x1,c) in X1_Color:
                r.append((z,x1,c,mjd))

    params=np.rec.array(r, dtype=[('z', 'f8'),('X1', 'f8'), ('Color', 'f8'),('DayMax','f8')])

    name='Light_Curves/Simul_'+str(fieldid)+'_Cadence_'+str(cadence)+'_X1_'+str(X1_Color[0][0])+'_Color_'+str(X1_Color[0][1])+'.pkl'

    #simulate
    Simul_LCs(params,obs,telescope,fname=name)
    """
    Simulation(fieldid,cadence,X1,Color,telescope)
#Analyze LC
#Ana_LCs(name,telescope)

if process_step == 'fit_ana':
    fname='Light_Curves/Simul_'+str(fieldid)+'_Cadence_'+str(cadence)+'_X1_'+str(X1)+'_Color_'+str(Color)+'.pkl'
    Fit_Ana(fname,telescope)

"""
lc_ana={}
for (X1,Color) in [(0.0,0.0),(-2.0,0.2),(2.0,-0.2)]:
    fichname='Ana_Simul_100_Cadence_3X1_'+str(X1)+'_Color_'+str(Color)+'.pkl'
    lc_ana[(X1,Color)]=pkl.load(open(fichname,'r'))
#print(lc_ana['fit_status'])
#Plot(lc_ana)


#print(lc_ana.dtype)
selection_names=['Sela_cc','Selb_10_2','Selb_10_2_cc']
selection={}
for selname in selection_names:
    selection[selname]=Load_Selection(selname+'.txt')
    figa, axa = plt.subplots()
    for (X1,Color) in [(0.0,0.0),(-2.0,0.2),(2.0,-0.2)]:
        effi=Calc_Effi(Presel(lc_ana[(X1,Color)]),selection[selname])
        Plot_Effi(axa,effi)
plt.show()

"""
