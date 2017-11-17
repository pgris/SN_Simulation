import numpy as np
from LC_Ana import *
from optparse import OptionParser
from ID import *
import cPickle as pkl
import matplotlib.pyplot as plt
from scipy.spatial import distance

def Plot(vara,varb,name):
    nrow=3
    ncol=2
    fig, axes = plt.subplots(nrow, ncol, figsize=(10,9))
    
    bands=dict(zip(['00','01','10','11','20','21'],['g','r','i','z','y','None']))
    
    for axnum in range(ncol * nrow):
        row = axnum // ncol
        col = axnum % ncol
        ax = axes[row, col]
        band=bands[str(row)+str(col)]
        print 'yes',row,col,band
        idx = res['band'] == band
        sel = res[idx]
        if len(sel) > 0:
            ax.hist(sel[vara]/sel[varb])
            ax.set_ylabel('Number of Entries')
            ax.set_xlabel(name)


parser = OptionParser()

parser.add_option("-f", "--fieldname", type="string", default='DD', help="filter [%default]")
parser.add_option("-i", "--fieldid", type="int", default=744, help="filter [%default]")
parser.add_option("-x", "--stretch", type="float", default=-999., help="filter [%default]")
parser.add_option("-c", "--color", type="float", default=-999., help="filter [%default]")
parser.add_option("-d", "--dirmeas", type="string", default="DD", help="filter [%default]")
parser.add_option("--dirobs", type="string", default="DD", help="filter [%default]")

opts, args = parser.parse_args()

dict_data={}

fieldid=opts.fieldid
fieldname=opts.fieldname
thedir=opts.dirmeas
thedir_obs=opts.dirobs
X1=opts.stretch
Color=opts.color

for seas in [0]:
    dict_data['DD_'+str(fieldid)+'_'+str(seas+1)]=Id(thedir=thedir,thedir_obs=thedir_obs,fieldname=fieldname,fieldid=fieldid,X1=X1,Color=Color,season=seas,colorfig='k')

#LC_Ana(dict_data,zmin=0.01,zmax=1.,bin_z=0.02)

z=0.11

filename=fieldname+'/'+str(fieldid)+'/Season_0/'+fieldname+'_'+str(fieldid)+'_'+str(z)+'_X1_2.0_C_-0.2_0_949.pkl'

print filename
lc_sncosmo=pkl.load(open('../../Light_Curves_sncosmo_noairmass/'+filename,'rb'))
lc_snsim=pkl.load(open('../../Light_Curves_snsim_noairmass/'+filename,'rb'))

print len(lc_sncosmo),len(lc_snsim)

r=[]
for i in range(len(lc_sncosmo)):
    lc_cos=lc_sncosmo[i]

    points=[(lc_snsim[i].meta['DayMax'],lc_snsim[i].meta['X1'],lc_snsim[i].meta['Color']) for i in range(len(lc_snsim))]
    idx_nearest=distance.cdist([(lc_cos.meta['DayMax'],lc_cos.meta['X1'],lc_cos.meta['Color'])],points).argmin()

    lc_sim=lc_snsim[idx_nearest]
    
    diff={}
    vals=['DayMax','X1','Color','z']
    for val in vals:
        diff[val]=lc_snsim[idx_nearest].meta[val]-lc_cos.meta[val]

    for band in 'grizy':
        idx = lc_cos['band']=='LSST::'+band
        lc_cosf=lc_cos[idx]
        idx = lc_sim['band']=='LSST::'+band
        lc_simf=lc_sim[idx]
        for i, val in enumerate(lc_cosf):
            idf=lc_simf['time']==val['time']
            lc_simff=lc_simf[idf]
            r.append((band,lc_cos.meta['DayMax'],lc_cos.meta['X1'],lc_cos.meta['Color'],lc_cos.meta['z'],diff['DayMax'],diff['X1'],diff['Color'],diff['z'],val['flux_e_sec'],lc_simff['flux_e_sec'][0],val['flux_e_sec']*val['fluxerr']/val['flux'],lc_simff['err_flux_e_sec'][0],val['time']))
    """
    fig, axes = plt.subplots(ncols=1, nrows=1, figsize=(10,9))

    axes.plot(lc_cos['time'],lc_cos['flux_e_sec'],'bo')
    axes.plot(lc_sim['time'],lc_sim['flux_e_sec'],'ko')

    print 100.*(lc_cos['flux_e_sec']-lc_sim['flux_e_sec'])/lc_cos['flux_e_sec']
    plt.show()
    """
    #print r
    #break

res=np.rec.fromrecords(r,names=['band','DayMax','X1','Color','z','diff_DayMax','diff_X1','diff_Color','diff_z','flux_sncosmo','flux_snsim','err_flux_sncosmo','err_flux_snsim','time'])

#print res

Plot('flux_sncosmo','flux_snsim','flux$^{sncosmo}$/flux$^{snsim}$')
Plot('err_flux_sncosmo','err_flux_snsim','errflux$^{sncosmo}$/errflux$^{snsim}$')
plt.show()


   
