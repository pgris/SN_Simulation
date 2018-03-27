import numpy as np
from Telescope import *
import matplotlib.pyplot as plt

telescope=Telescope(atmos=True,aerosol=False,airmass=1.2)
mean_waves={}
for band in 'ugrizy':
    print(band,telescope.throughputs.mean_wavelength[band])
    mean_waves[band]=telescope.throughputs.mean_wavelength[band]


zvals=np.arange(0.,1.501,0.01)
fontsize=12
tot_label=[]
for band in 'ugrizy':
    tot_label.append(plt.errorbar(zvals,mean_waves[band]/(1.+zvals),label=band+' band'))

plt.plot(zvals,[300.]*len(zvals),color='k')
plt.plot(zvals,[800.]*len(zvals),color='k')
plt.xlabel(r'z',{'fontsize': 2*fontsize})
plt.ylabel(r'$\frac{\lambda_{mean}}{(1+z)}$ [nm]',{'fontsize': 2*fontsize})
plt.xlim([0.,1.5])
labs = [l.get_label() for l in tot_label]
plt.legend(tot_label, labs, ncol=6,loc='upper right',prop={'size':1.5*fontsize},frameon=False)
plt.xticks(size=1.5*fontsize)
plt.yticks(size=1.5*fontsize)
plt.show()
