import numpy as np
from Telescope import *
import matplotlib.pyplot as plt
from astropy.table import Table

t = Table.read('throughputs_snsim.txt',format='ascii')
print t

telescope=Telescope(atmos=True,aerosol=False,airmass=1.2)


style=[',',',',',',',']
for i,band in enumerate(['u','g','r','i','z','y']):
            #plt.plot(telescope.throughputs.lsst_std[band].wavelen,telescope.throughputs.lsst_std[band].sb,linestyle='-',color=telescope.throughputs.filtercolors[band], label='%s - std' %(band))
    plt.plot(telescope.throughputs.lsst_system[band].wavelen,telescope.throughputs.lsst_system[band].sb,linestyle='--',color=telescope.throughputs.filtercolors[band], label='%s - syst' %(band))
    plt.plot(telescope.throughputs.lsst_atmos[band].wavelen,telescope.throughputs.lsst_atmos[band].sb,linestyle='-.',color=telescope.throughputs.filtercolors[band], label='%s - syst+atm' %(band))
    """
    if len(telescope.throughputs.lsst_atmos_aerosol) > 0:
        plt.plot(telescope.throughputs.lsst_atmos_aerosol[band].wavelen,telescope.throughputs.lsst_atmos_aerosol[band].sb,linestyle='-',color=telescope.throughputs.filtercolors[band], label='%s - syst+atm+aero' %(band))
        plt.plot(telescope.throughputs.atmos.wavelen, telescope.throughputs.atmos.sb, 'k:', label='X =%.1f atmos' %(telescope.throughputs.airmass),linestyle='-')
        if len(telescope.throughputs.lsst_atmos_aerosol) > 0:
            plt.plot(telescope.throughputs.atmos_aerosol.wavelen, telescope.throughputs.atmos_aerosol.sb, 'k:', label='X =%.1f atm+aero' %(telescope.throughputs.airmass),linestyle='--')
     """
    idx=t['band']==band
    plt.plot(t[idx]['wave']/10.,t[idx]['trans'],linestyle='-',color=telescope.throughputs.filtercolors[band])
    plt.plot(t[idx]['wave']/10.,t[idx]['trans']*t[idx]['atmos'],linestyle='-',color=telescope.throughputs.filtercolors[band])

plt.legend(loc=(0.85, 0.1), fontsize='smaller', fancybox=True, numpoints=1)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Sb (0-1)')
plt.title('System throughput')


plt.show()
