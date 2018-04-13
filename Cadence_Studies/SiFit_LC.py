from Telescope import *
import cPickle as pkl
import matplotlib.pyplot as plt
import time
import operator
from optparse import OptionParser
from Simulation import *
from Fit_Ana import *


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
    idx=(tab['phase_max_5_all']>=10.)&(tab['phase_min_5_all']<=-5.)
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
        if len(sel_z) > 0.:
        #print('hello',z,effi,err,len(sel))
            sel_z.sort(order='DayMax')
            #idb = np.abs(sel['phase_min_5_all']-5.)< 0.01
            #print('hello',sel['phase_min_5_all'])
            #idb = np.searchsorted(sel['phase_min_5_all'],-5.0)
            #X = np.abs(sel_z['phase_min_5_all']+5.0)
            idb = sel_z['phase_min_5_all']<= -5.
            #Y = sel_z['phase_max_5_all']
            idc = sel_z['phase_max_5_all'] >= 10.
            mjds_min=(1.+z)*sel_z[idb]['phase_min_5_all']+sel_z[idb]['DayMax']
            #print('aaa',sel_z[idb]['DayMax'],sel_z[idb]['phase_min_5_all'])
            print('hello',z,np.min(sel_z[idb]['DayMax']),np.max(sel_z[idc]['DayMax']))
            #mjds_min=(1.+z)*sel[idb]['phase_min_5_all']+sel[idb]['DayMax']
            #mjds_max=(1.+z)*sel['phase_max_all']+sel['DayMax']
            #print(z,np.max(mjds_min),np.min(mjds_max))
            #print((1.+z)*sel['phase_min_all'],sel['DayMax'])
        ro.append((z,effi,err))

    return np.rec.fromrecords(ro,names=['z','effi','err_effi'])

def Plot_Effi(axa,tab,ls='-',color='k'):

    #figa, axa = plt.subplots()
    axa.errorbar(tab['z'],tab['effi'],yerr=tab['err_effi'],ls=ls,color=color)


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
parser.add_option('--fieldids', default='',type="string")
parser.add_option('--cadences', default='',type="string")

opts, args = parser.parse_args()

process_step=opts.process_step
fieldid=opts.fieldid
cadence=opts.cadence
X1=opts.X1
Color=opts.Color
fieldids=opts.fieldids
cadences=opts.cadences

#Load instrument

atmos=True
telescope=Telescope(atmos=atmos,aerosol=False,airmass=1.2)

if process_step == 'simulation':
    Simulation(fieldid,cadence,X1,Color,telescope)

if process_step == 'fit_ana':
    fname='Light_Curves/Simul_'+str(fieldid)+'_Cadence_'+str(cadence)+'_X1_'+str(X1)+'_Color_'+str(Color)+'.pkl'
    Fit_Ana(fname,telescope)

if process_step == 'analysis':
    lc_ana={}

    if fieldids != '':
        for i,field in enumerate(fieldids.split(',')):
            cad=cadences.split(',')
            lc_ana[(int(field),int(cad[i]))]={}
    else:
        lc_ana[(fieldid,cadence)]={}

    for key in lc_ana.keys():
        for (X1,Color) in [(0.0,0.0),(-2.0,0.2),(2.0,-0.2)]:
            fichname='SN_from_LC/Simul_'+str(key[0])+'_Cadence_'+str(key[1])+'_X1_'+str(X1)+'_Color_'+str(Color)+'.pkl'
            lc_ana[key][(X1,Color)]=pkl.load(open(fichname,'r'))
        #print(lc_ana['fit_status'])
        #Plot(lc_ana)


        #print(lc_ana.dtype)
    #selection_names=['Sela_cc','Selb_10_2','Selb_10_2_cc']
    selection_names=['Sela_cc']
    selection={}
    lsty=dict(zip([(100,3),(101,6)],['--',':']))
    colors=dict(zip([(0.0,0.0),(-2.0,0.2),(2.0,-0.2)],['k','r','b']))
    for selname in selection_names:
        selection[selname]=Load_Selection(selname+'.txt')
        figa, axa = plt.subplots()
        for key in lc_ana.keys():
            #for (X1,Color) in [(0.0,0.0),(-2.0,0.2),(2.0,-0.2)]:
            for (X1,Color) in [(-2.0,0.2)]:
                effi=Calc_Effi(Presel(lc_ana[key][(X1,Color)]),selection[selname])
                Plot_Effi(axa,effi,ls=lsty[key],color=colors[(X1,Color)])
    
    plt.show()
