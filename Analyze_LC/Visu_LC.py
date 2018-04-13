import numpy as np
from optparse import OptionParser
import cPickle as pkl
import glob
from astropy.table import Table,vstack
from scipy.spatial import distance
import sncosmo
from Telescope import *
from astropy import (cosmology, units as u, constants as const)
import matplotlib.pyplot as plt
from matplotlib import cm
import h5py

class Visu_LC:
    def __init__(self,fieldname='DD',fieldid=744,X1=0.,Color=0.,season=0,z=0.1,DayMax=0.2,data_dir='../Fitted_Light_Curves',lc_dir='../Light_Curves'):

        self.data_dir=data_dir
        
        dirmeas=self.data_dir+'/'+fieldname+'/'+str(fieldid)+'/Season_'+str(season)
        #files = glob.glob(dirmeas+'/'+fieldname+'_'+str(fieldid)+'*_'+str(T0min)+'_'+str(T0max)+'.pkl')
        fichid=fieldname+'_'+str(fieldid)+'_'+str(z)+'_X1_'+str(X1)+'_C_'+str(Color)+'*.hdf5'
        print 'looking for',dirmeas+'/'+fichid,DayMax
        files = glob.glob(dirmeas+'/'+fichid)

        lc=None
        for fi in files:
            f = h5py.File(fi,'r')
            print 'loading',fi,len(f.keys())
            for i,key in enumerate(f.keys()):
                tab=Table.read(f, path=key)
                """
                for val in tab:
                    print val['DayMax'],float(DayMax),val['DayMax']-float(DayMax)
                """
                idxb = np.abs(tab['DayMax']-DayMax)<1.e-5
                sel=tab[idxb]
                if len(sel) > 0:
                    print 'yes found',len(tab)
                    lc=sel
                    idx=(sel['N_bef_all']>=4)&(sel['N_aft_all']>=10)
                    idx&=(sel['status']=='go_fit')&(sel['fit_status']=='fit_ok')
                    idx&=(sel['phase_first']<=-5)&(sel['phase_last']>=20)
                    print 'Selection',len(sel[idx]),sel['N_bef_all'],sel['N_aft_all'],sel['status'],sel['fit_status'],sel['phase_first'],sel['phase_last']
                    break

        #print len(lc),lc
        #grab corresponding LC

        dir_lc=dirmeas.replace(data_dir,lc_dir)

        tot_lc=[]
        files = glob.glob(dir_lc+'/'+fichid)
        lcpoints=None
        for fi in files:
            print 'loading',fi    
            f = h5py.File(fi,'r')
            for i,key in enumerate(f.keys()):
                tab=Table.read(f, path=key)
                #print tab
                if np.abs(tab.meta['DayMax']-DayMax)<1.e-5:
                    print 'yes found',len(tab)
                    lcpoints=tab
                    break 

        #print lcpoints

        
        telescope=Telescope(airmass=np.median(1.2))

        transmission=telescope.throughputs

        for filtre in 'ugrizy':
            band=sncosmo.Bandpass(transmission.atmosphere[filtre].wavelen,transmission.atmosphere[filtre].sb, name='LSST::'+filtre,wave_unit=u.nm)
            sncosmo.registry.register(band)

        self.Plot_LC_Points(data=lcpoints,flux_name='flux')
        self.Plot_LC(lcpoints,lc)
        plt.show()
        """

        for lc in tot_lc:
            if lc.meta['X1']==sn_id['X1'] and lc.meta['Color']==sn_id['Color'] and lc.meta['DayMax']==sn_id['DayMax']:
                print 'yes found'
                print lc
                self.Plot_LC_Points(data=lc,flux_name='flux_e_sec')
                self.Plot_LC(lc,sn_id)
                break

        plt.show()
        """

    def Plot_LC_Points(self,data=None, flux_name='flux',xfigsize=None, yfigsize=None, figtext=None,figtextsize=1.,ncol=2,color=None,cmap=None, cmap_lims=(3000., 10000.),zp=25., zpsys='ab',):
       
        

        toff=data.meta['DayMax']
        bands=set(data['band'])

        tmin, tmax = [], []
        if data is not None:
            tmin.append(np.min(data['time']) - 10.)
            tmax.append(np.max(data['time']) + 10.)

        # Calculate layout of figure (columns, rows, figure size). We have to
        # calculate these explicitly because plt.tight_layout() doesn't space the
        # subplots as we'd like them when only some of them have xlabels/xticks.
        wspace = 0.6  # All in inches.
        hspace = 0.3
        lspace = 1.0
        bspace = 0.7
        trspace = 0.2
        nrow = (len(bands) - 1) // ncol + 1
        if xfigsize is None and yfigsize is None:
            hpanel = 2.25
            wpanel = 3.
        elif xfigsize is None:
            hpanel = (yfigsize - figtextsize - bspace - trspace -
                      hspace * (nrow - 1)) / nrow
            wpanel = hpanel * 3. / 2.25
        elif yfigsize is None:
            wpanel = (xfigsize - lspace - trspace - wspace * (ncol - 1)) / ncol
            hpanel = wpanel * 2.25 / 3.
        else:
            raise ValueError('cannot specify both xfigsize and yfigsize')

        figsize = (lspace + wpanel * ncol + wspace * (ncol - 1) + trspace,
               bspace + hpanel * nrow + hspace * (nrow - 1) + trspace +
               figtextsize)  

        # Create the figure and axes.
        fig, axes = plt.subplots(nrow, ncol, figsize=figsize, squeeze=False)
        
        fig.subplots_adjust(left=lspace / figsize[0],
                            bottom=bspace / figsize[1],
                            right=1. - trspace / figsize[0],
                            top=1. - (figtextsize + trspace) / figsize[1],
                            wspace=wspace / wpanel,
                            hspace=hspace / hpanel)
        # Color options.
        if color is None:
            if cmap is None:
                cmap = cm.get_cmap('jet_r')
        # Loop over bands
        bands = list(bands)
        waves = [sncosmo.get_bandpass(b).wave_eff for b in bands]
        waves_and_bands = sorted(zip(waves, bands))

        for axnum in range(ncol * nrow):
            row = axnum // ncol
            col = axnum % ncol
            ax = axes[row, col]
 
            if axnum >= len(waves_and_bands):
                ax.set_visible(False)
                ax.set_frame_on(False)
                continue

            wave, band = waves_and_bands[axnum]

            bandname_coords = (0.92, 0.92)
            bandname_ha = 'right'
            if color is None:
                bandcolor = cmap((cmap_lims[1] - wave) /
                                 (cmap_lims[1] - cmap_lims[0]))
            else:
                bandcolor = color

            if data is not None:
                mask = data['band'] == band
                time = data['time'][mask]
                flux = data[flux_name][mask]
                fluxerr = data['fluxerr'][mask]/data['flux'][mask]*data[flux_name][mask]
                ax.errorbar(time - toff, flux, fluxerr, ls='None',
                            color=bandcolor, marker='.', markersize=10.)    
            # Band name in corner
            ax.text(bandname_coords[0], bandname_coords[1], band,
                    color='k', ha=bandname_ha, va='top', transform=ax.transAxes)

            ax.axhline(y=0., ls='--', c='k')  # horizontal line at flux = 0.
            ax.set_xlim((tmin-toff, tmax-toff))

            if (len(bands) - axnum - 1) < ncol:
                if toff == 0.:
                    ax.set_xlabel('time')
                else:
                    ax.set_xlabel('time - {0:.2f}'.format(toff))
            else:
                for l in ax.get_xticklabels():
                    l.set_visible(False)
            if col == 0:
                if flux_name == 'flux':
                    ax.set_ylabel('flux ($ZP_{{{0}}} = {1}$)'
                                  .format(sncosmo.get_magsystem(zpsys).name.upper(), zp))
                if flux_name == 'flux_e_sec':
                   ax.set_ylabel('flux (e/sec)') 


    def Plot_LC(self,lc,sn):

        dust = sncosmo.OD94Dust()
        fitted_model=sncosmo.Model(source='salt2-extended', effects=[dust, dust],
                                   effect_names=['host', 'mw'],
                                   effect_frames=['rest', 'obs'])
        fitted_model.set(z=sn['z'])
        fitted_model.set(t0=sn['salt2.T0'])
        fitted_model.set(x0=sn['salt2.X0'])
        fitted_model.set(x1=sn['salt2.X1'])
        fitted_model.set(c=sn['salt2.Color']) 
        
        errors={}
        errors['t0']=np.sqrt(sn['salt2.CovT0T0'])
        errors['x0']=np.sqrt(sn['salt2.CovX0X0'])
        errors['x1']=np.sqrt(sn['salt2.CovX1X1'])
        errors['c']=np.sqrt(sn['salt2.CovColorColor'])
        
        print 'Phases : first',(lc['time'][0]-sn['salt2.T0'])/(1.+sn['z']),'last',(lc['time'][-1]-sn['salt2.T0'])/(1.+sn['z'])
        """
        res, fitted_modelb = sncosmo.fit_lc(lc, fitted_model,['t0', 'x0', 'x1', 'c'],bounds={'z':(sn['z']-0.001, sn['z']+0.001)})

        print 'ooo',res.errors
        print 'bbb',errors
        """
        sncosmo.plot_lc(lc, model=fitted_model,pulls=True)
            
        plt.show()

parser = OptionParser()
parser.add_option("-f", "--fieldname", type="string", default='DD', help="filter [%default]")
parser.add_option("-i", "--fieldid", type="int", default=744, help="filter [%default]")
parser.add_option("-x", "--stretch", type="float", default=0., help="filter [%default]")
parser.add_option("-c", "--color", type="float", default=0., help="filter [%default]")
parser.add_option("-d", "--dirmeas", type="string", default="DD", help="filter [%default]")
parser.add_option("--T0step", type="float", default=0.2, help="filter [%default]")
parser.add_option("--dirobs", type="string", default="DD", help="filter [%default]")
parser.add_option("--DayMax", type="float", default=0.2, help="filter [%default]")
parser.add_option("--z", type="float", default=0., help="filter [%default]")

opts, args = parser.parse_args()

fieldid=opts.fieldid
fieldname=opts.fieldname
thedir=opts.dirmeas
thedir_obs=opts.dirobs
X1=opts.stretch
Color=opts.color
z=opts.z
T0step=opts.T0step
DayMax=opts.DayMax

prefix='/sps/lsst/data/dev/pgris/'
Visu_LC(fieldname=fieldname,fieldid=fieldid,X1=X1,Color=Color,z=z,DayMax=DayMax,data_dir=prefix+'Fitted_Light_Curves_sncosmo_testb_'+str(T0step).replace('.','_')+'_b',lc_dir=prefix+'Light_Curves_sncosmo_testb_'+str(T0step).replace('.','_')+'_b')
