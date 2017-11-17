import numpy as np
from optparse import OptionParser
import cPickle as pkl
import glob
from astropy.table import vstack
from scipy.spatial import distance
import sncosmo
from Telescope import *
from astropy import (cosmology, units as u, constants as const)
import matplotlib.pyplot as plt
from matplotlib import cm

class Visu_LC:
    def __init__(self,fieldname='DD',fieldid=744,X1=0.,Color=0.,season=0,z=0.1,T0min=900,T0max=920,data_dir='../Fitted_Light_Curves',lc_dir='../Light_Curves'):

        self.data_dir=data_dir

        dirmeas=self.data_dir+'/'+fieldname+'/'+str(fieldid)+'/Season_'+str(season)
        files = glob.glob(dirmeas+'/'+fieldname+'_'+str(fieldid)+'*_'+str(T0min)+'_'+str(T0max)+'.pkl')

        tot_fit=None
        for fi in files:
            pkl_file = open(fi,'rb')
            print 'loading',fi
            if tot_fit is None:
                tot_fit=pkl.load(pkl_file)
            else:
                tot_fit=vstack([tot_fit,pkl.load(pkl_file)])
            
        print tot_fit

        points=[(tot_fit['X1'][i],tot_fit['Color'][i]) for i in range(len(tot_fit))]
        
        idx=distance.cdist([(X1,Color)],points).argmin()

        print idx,tot_fit[idx]

        sn_id=tot_fit[idx]

        #grab corresponding LC
        dir_lc=dirmeas.replace(data_dir,lc_dir)

        tot_lc=[]
        files = glob.glob(dir_lc+'/'+fieldname+'_'+str(fieldid)+'*_'+str(T0min)+'_'+str(T0max)+'.pkl')
        for fi in files:
            pkl_file = open(fi,'rb')
            print 'loading',fi
            if not tot_lc:
                tot_lc=pkl.load(pkl_file)
            else:
                tot_lc=tot_lc+pkl.load(pkl_file)

        print 'hello',len(tot_lc)


        telescope=Telescope(airmass=np.median(1.2))

        transmission=telescope.throughputs

        for filtre in 'ugrizy':
            band=sncosmo.Bandpass(transmission.atmosphere[filtre].wavelen,transmission.atmosphere[filtre].sb, name='LSST::'+filtre,wave_unit=u.nm)
            sncosmo.registry.register(band)


        for lc in tot_lc:
            if lc.meta['X1']==sn_id['X1'] and lc.meta['Color']==sn_id['Color'] and lc.meta['DayMax']==sn_id['DayMax']:
                print 'yes found'
                print lc
                self.Plot_LC_Points(data=lc,flux_name='flux_e_sec')
                self.Plot_LC(lc,sn_id)
                break

        plt.show()

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
        
        """
        res, fitted_modelb = sncosmo.fit_lc(lc, fitted_model,['t0', 'x0', 'x1', 'c'],bounds={'z':(sn['z']-0.001, sn['z']+0.001)})

        print 'ooo',res.errors
        print 'bbb',errors
        """
        sncosmo.plot_lc(lc, model=fitted_model,pulls=True,errors=errors)
            
        plt.show()

parser = OptionParser()
parser.add_option("-f", "--fieldname", type="string", default='DD', help="filter [%default]")
parser.add_option("-i", "--fieldid", type="int", default=744, help="filter [%default]")
parser.add_option("-x", "--stretch", type="float", default=0., help="filter [%default]")
parser.add_option("-c", "--color", type="float", default=0., help="filter [%default]")
parser.add_option("-d", "--dirmeas", type="string", default="DD", help="filter [%default]")
parser.add_option("--dirobs", type="string", default="DD", help="filter [%default]")

opts, args = parser.parse_args()

fieldid=opts.fieldid
fieldname=opts.fieldname
thedir=opts.dirmeas
thedir_obs=opts.dirobs
X1=opts.stretch
Color=opts.color
z=0.1

Visu_LC(fieldname=fieldname,fieldid=fieldid,X1=0,Color=0.,z=z)
