import numpy as np
import pickle as pkl
from mpl_toolkits.basemap import Basemap
import pylab as plt


def Plot_Ra_Dec():
    fieldnames=['WFD','NorthEclipticSpur-18c','SouthCelestialPole-18','GalacticPlane','DD']
    colors=['b','y','g','k','r']

    corresp=dict(zip(fieldnames,['WFD','NorthEclipticSpur','SouthCelestialPole','GalacticPlane','DD']))
    m = Basemap(projection='moll',lon_0=180)
    parallels = np.arange(-90.,90,30.)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
    meridians = np.arange(0.,360.,60.)
    m.drawmeridians(meridians,labels=[0,1,1,0],fontsize=10)
    for i in np.arange(len(meridians)):
        plt.annotate(np.str(int(meridians[i]))+'$^o$',xy=m(meridians[i],30),xycoords='data')

    for i,fieldname in enumerate(fieldnames):
        fich=fieldname+'_Ra_Dec.txt'
        res=np.loadtxt(fich,dtype={'names': ('fieldid','Ra','Dec'),'formats': ('i4','f8','f8')})
    #print res['Ra']
        x, y = m(np.rad2deg(res['Ra']),np.rad2deg(res['Dec']))
        m.scatter(x,y,marker='s',color=colors[i], label=corresp[fieldname])
    

    plt.legend(bbox_to_anchor=(0.05, -0.20), ncol=3,loc=2, borderaxespad=0.,fontsize=12.)   
    plt.gcf().savefig('Map_Plots/Ra_Dec.png')

def Plot_Hist(what,tab_all,season):

    conv=dict(zip([band for band in 'ugrizy'],[1,2,3,4,5,6]))
    class_map = lambda x: conv[x]
    figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
    for band in 'ugrizy':
        idx = (tab_all['band']== band)&(tab_all['season']== season)&(tab_all['median_cadence']>0.)
        sel=tab_all[idx]
        print len(sel),season,band,sel['median_cadence'],sel['fieldid'],np.median(sel['median_cadence'])
        axa.hist(sel['median_cadence'],histtype='step',bins=35)
        #axa.scatter(sel['median_cadence'],map(class_map,sel['band']))

    axa.set_xlim([0,35])

def Plot_Map(what,tab_all,band,season,thresh,legend):

#print tab['ra']


    idx = (tab_all['band']==band)&(tab_all['season']==season)&(tab_all[what]>thresh)
    tab=tab_all[idx]
    
    print tab['fieldid'],len(tab['fieldid'])

    lons = np.rad2deg(tab['ra'])
    lats = np.rad2deg(tab['dec'])

    m = Basemap(projection='moll',lon_0=180)

    x, y = m(lons,lats)
#x, y = m(*np.meshgrid(lons, lats))
#m.drawmapboundary(fill_color='#99ffff')
    parallels = np.arange(-90.,90,30.)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
    meridians = np.arange(0.,360.,60.)
    m.drawmeridians(meridians,labels=[0,1,1,0],fontsize=10)
    for i in np.arange(len(meridians)):
        plt.annotate(np.str(int(meridians[i]))+'$^o$',xy=m(meridians[i],10),xycoords='data')
#m.fillcontinents(color='#cc9966',lake_color='#99ffff')
    m.scatter(x,y,c=tab[what],marker='s',s=10,cmap=plt.cm.jet)
#print 'hello',len(x),tab['median_cadence']
#cs = m.contourf(x,y,tab['median_cadence'],cmap=plt.cm.jet)
#m.pcolormesh(x, y,tab['median_cadence'] , cmap='Paired')
    toprint=legend+' - '+band+' band - season '+str(season+1)
    plt.annotate(toprint, xy=(0.30, 1.1), xycoords='axes fraction')
    plt.colorbar(fraction=0.02, pad=0.04)
    plt.grid()

def Olp_Plot():
    tab_all=None

    for fieldname in ['WFD']:
        if tab_all is None:
            tab_all=pkl.load(open(fieldname+'.pkl','rb'))
        else:
            tab_all=np.concatenate((tab_all,pkl.load(open(fieldname+'.pkl','rb'))))


    what='median_cadence'
    for season in range(10):
    #Plot_Map('m5',tab_all,'r',season,10.,'5-$\sigma$ depth')
        """
        for band in 'ugrizy':
        Plot_Map('median_cadence',tab_all,band,season,0.,'Cadence')
        #Plot_Hist('median_cadence',tab_all,season)
        plt.show()
        """
        for band in 'ugrizy':
            idx = (tab_all['band']==band)&(tab_all['season']==season)&(tab_all[what]>-0.5)
            tab=tab_all[idx]
            print season,band,np.percentile(np.array(tab['median_cadence']),50.),np.median(np.array(tab['median_cadence']))
         
    idb = (tab_all['fieldid'] == 311)&(tab_all['season'] == 0)
    tab=tab_all[idb]
    print len(tab),tab[['fieldid','band','median_cadence']]

#Plot_Ra_Dec()
def Plot_Cadence(tab):

    radec=np.unique((tab[['Ra','Dec']]))
    print radec[:]
    
    norm=[]
    for val in radec:
        idd = (tab['Ra']==val[0])&(tab['Dec']==val[1])
        sel=tab[idd]
        #norm.append(np.median(sel['mjd'][1:]-sel['mjd'][:-1]))
        norm.append(np.median(sel['Nvisits']))

    lons = np.rad2deg([radec[i][0] for i in range(len(radec))])
    lats = np.rad2deg([radec[i][1] for i in range(len(radec))])
    
    m = Basemap(projection='moll',lon_0=180)
    
    x, y = m(lons,lats)
    
    m.scatter(x,y,c=norm,marker='s',s=10,cmap=plt.cm.jet)
    
   
    cbar = plt.colorbar(fraction=0.02, pad=0.04)
    cbar.ax.set_ylabel('Median number of Visits')

    plt.grid()

fieldname='WFD'

tab_all=pkl.load(open('Cadence_'+fieldname+'.pkl','rb'))
bands='g'

for season in range(1):
    for band in bands:
        idx = (tab_all['season']==season)&(tab_all['band']=='LSSTPG::'+band)
        tab=tab_all[idx]
        Plot_Cadence(tab)
        plt.show()
