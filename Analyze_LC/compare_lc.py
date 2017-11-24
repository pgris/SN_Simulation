import numpy as np
from LC_Ana import *
from optparse import OptionParser
from ID import *
import cPickle as pkl
import matplotlib.pyplot as plt
from scipy.spatial import distance

def Index_nearest(lca,lcb):
    
    points=[(lcb[i].meta['DayMax'],lcb[i].meta['X1'],lcb[i].meta['Color']) for i in range(len(lcb))]
    idx_nearest=distance.cdist([(lca.meta['DayMax'],lca.meta['X1'],lca.meta['Color'])],points).argmin()

    return idx_nearest

def Plot(var,vara,varb,varc,name):
    nrow=3
    ncol=2
    fig, axes = plt.subplots(nrow, ncol, figsize=(10,9))
    fontsize=12
    
    bands=dict(zip(['00','01','10','11','20','21'],['g','r','i','z','y','None']))
    
    for axnum in range(ncol * nrow):
        row = axnum // ncol
        col = axnum % ncol
        ax = axes[row, col]
        band=bands[str(row)+str(col)]
        print 'yes',row,col,band
        idx = (res['band'] == band)&(res['time']>2.)
        sel = res[idx]
        if len(sel) > 0:
            ll=vara+'/'+varb
            ax.hist(sel[var+'_'+vara]/sel[var+'_'+varb],histtype='step',label=ll,normed=True,stacked=True,bins=15)
            if varc is not None:
                ll=vara+'/'+varc
                ax.hist(sel[var+'_'+vara]/sel[var+'_'+varc],histtype='step',label=ll,normed=True,stacked=True,bins=15)
                ll=varc+'/'+varb
                ax.hist(sel[var+'_'+varc]/sel[var+'_'+varb],histtype='step',label=ll,normed=True,stacked=True,bins=15)
            ax.set_ylabel('Number of Entries')
            ax.set_xlabel(name)
            
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            ax.text(xlim[0]+0.95*(xlim[1]-xlim[0]),ylim[0]+0.92*(ylim[1]-ylim[0]),band)
            if row == 0 and col == 0:
                ax.legend(loc='upper left',prop={'size':fontsize},frameon=False)

def Plot_vs(vara,varb,varc,name):
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
        idx = (res['band'] == band)&(res['time']>2.)
        sel = res[idx]
        if len(sel) > 0:
            ax.plot(sel[varc],sel[vara]/sel[varb],'k.')
            ax.set_ylabel(name)
            ax.set_xlabel(varc)
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            #print 'h',0.98*xlim[1],1.1*ylim[0]
            ax.text(xlim[0]+0.95*(xlim[1]-xlim[0]),ylim[0]+0.1*(ylim[1]-ylim[0]),band)

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

z=0.5

filename=fieldname+'/'+str(fieldid)+'/Season_0/'+fieldname+'_'+str(fieldid)+'_'+str(z)+'_X1_2.0_C_-0.2_0_949.pkl'

print filename
lc_sncosmo=pkl.load(open('../../Light_Curves_sncosmo_test/'+filename,'rb'))
lc_snsim=pkl.load(open('../../Light_Curves_snsim_new/'+filename,'rb'))
#lc_snsim_psf_gauss=pkl.load(open('../../Light_Curves_snsim_new/'+filename,'rb'))

print len(lc_sncosmo),len(lc_snsim)

"""
fig, ax = plt.subplots(2, 2, figsize=(10,9))

ax[0][0].hist([lc_snsim[i].meta['DayMax'] for i in range(len(lc_snsim))])
ax[0][0].hist([lc_sncosmo[i].meta['DayMax'] for i in range(len(lc_sncosmo))])

ax[0][1].hist([lc_snsim[i].meta['z'] for i in range(len(lc_snsim))])
ax[0][1].hist([lc_sncosmo[i].meta['z'] for i in range(len(lc_sncosmo))])

ax[1][0].hist([lc_snsim[i].meta['X1'] for i in range(len(lc_snsim))])
ax[1][0].hist([lc_sncosmo[i].meta['X1'] for i in range(len(lc_sncosmo))])

ax[1][1].hist([lc_snsim[i].meta['Color'] for i in range(len(lc_snsim))])
ax[1][1].hist([lc_sncosmo[i].meta['Color'] for i in range(len(lc_sncosmo))])

plt.show()
"""

r=[]
for i in range(len(lc_sncosmo)):
    lc_cos=lc_sncosmo[i]

    idx_nearest=Index_nearest(lc_cos,lc_snsim)
    #idx_nearestb=Index_nearest(lc_cos,lc_snsim_psf_gauss)

    lc_sim=lc_snsim[idx_nearest]
    #lc_sim_b=lc_snsim_psf_gauss[idx_nearestb]

    #print 'oh yes',idx_nearest,idx_nearestb
    diff={}
    vals=['DayMax','X1','Color','z']
    for val in vals:
        diff[val]=lc_snsim[idx_nearest].meta[val]-lc_cos.meta[val]
        #print val,diff[val]

    pb=False
    for band in 'grizy':
        idx = lc_cos['band']=='LSST::'+band
        lc_cosf=lc_cos[idx]
        idx = lc_sim['band']=='LSST::'+band
        lc_simf=lc_sim[idx]
        """
        idxb = lc_sim_b['band']=='LSST::'+band
        lc_simfb=lc_sim_b[idxb]
        """
        """
        if band == 'y':
            print len(lc_cosf),len(lc_simf),len(lc_simfb)
        """
        for i, val in enumerate(lc_cosf):
            idf=lc_simf['time']==val['time']
            lc_simff=lc_simf[idf]
            #idfb=lc_simfb['time']==val['time']
            #lc_simffb=lc_simfb[idfb]
            """
            snr_sncosmo=val['flux']/val['fluxerr']
            snr_snsim=lc_simff['flux_e_sec'][0]/lc_simff['err_flux_e_sec'][0]
            snr_snsim_psf_gauss=lc_simffb['flux_e_sec'][0]/lc_simffb['err_flux_e_sec'][0]
            """
            snr_sncosmo=0
            snr_snsim=0
            snr_snsim_psf_gauss=0
            #print lc_simff['flux_e_sec'][0]
            #print len(lc_simffb),len(lc_cosf),len(lc_simff)
            #if len(lc_simffb) > 0:
            #r.append((band,lc_cos.meta['DayMax'],lc_cos.meta['X1'],lc_cos.meta['Color'],lc_cos.meta['z'],diff['DayMax'],diff['X1'],diff['Color'],diff['z'],val['flux_e_sec'],lc_simff['flux_e_sec'][0],val['flux_e_sec']*val['fluxerr']/val['flux'],lc_simff['err_flux_e_sec'][0],val['time']-np.min(lc_cosf['time']),lc_simffb['flux_e_sec'][0],lc_simffb['err_flux_e_sec'][0],snr_sncosmo,snr_snsim,snr_snsim_psf_gauss))
            r.append((band,lc_cos.meta['DayMax'],lc_cos.meta['X1'],lc_cos.meta['Color'],lc_cos.meta['z'],diff['DayMax'],diff['X1'],diff['Color'],diff['z'],val['flux_e_sec'],lc_simff['flux_e_sec'][0],val['flux_e_sec']*val['fluxerr']/val['flux'],lc_simff['err_flux_e_sec'][0],val['time']-np.min(lc_cosf['time'])))
          
    #break

    """
    if pb is True:
            figb, axb = plt.subplots(1, 2, figsize=(10,9))
            axb[0].plot(lc_cosf['time']-np.min(lc_cosf['time']),lc_cosf['flux_e_sec'],'ko')
            axb[0].plot(lc_simf['time']-np.min(lc_simf['time']),lc_simf['flux_e_sec'],'ro')
            axb[1].plot(lc_cosf['time'],lc_cosf['flux_e_sec']/lc_simf['flux_e_sec'],'ko')
            print lc_simf['time'],lc_cosf['time'],lc_cosf['flux_e_sec'],lc_simf['flux_e_sec']
            plt.show()
    """
    """
    fig, axes = plt.subplots(ncols=1, nrows=1, figsize=(10,9))

    axes.plot(lc_cos['time'],lc_cos['flux_e_sec'],'bo')
    axes.plot(lc_sim['time'],lc_sim['flux_e_sec'],'ko')

    print 100.*(lc_cos['flux_e_sec']-lc_sim['flux_e_sec'])/lc_cos['flux_e_sec']
    plt.show()
    """
    #print r
    #break

#res=np.rec.fromrecords(r,names=['band','DayMax','X1','Color','z','diff_DayMax','diff_X1','diff_Color','diff_z','flux_sncosmo','flux_snsim','err_flux_sncosmo','err_flux_snsim','time','flux_snsim_gain','err_flux_snsim_gain','snr_sncosmo','snr_snsim','snr_snsim_gain'])
res=np.rec.fromrecords(r,names=['band','DayMax','X1','Color','z','diff_DayMax','diff_X1','diff_Color','diff_z','flux_sncosmo','flux_snsim','err_flux_sncosmo','err_flux_snsim','time','flux_snsim_gain'])


#print res

Plot('flux','sncosmo','snsim',varc=None,name='flux ratio')
Plot('err_flux','sncosmo','snsim',varc=None,name='errflux ratio')
Plot_vs('flux_sncosmo','flux_snsim','time','flux$^{sncosmo}$/flux$^{snsim}$')
Plot_vs('err_flux_sncosmo','err_flux_snsim','time','errflux$^{sncosmo}$/errflux$^{snsim}$')

plt.show()


   
