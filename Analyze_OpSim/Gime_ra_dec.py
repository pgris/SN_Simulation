import numpy as np
import glob
from Observations import *
import pylab as plt
from shapely import geometry

simu_name='feature_baseline_10yrs'
fieldname='WFD'
season=2
thedir='OpSimLogs_'+simu_name+'/'+fieldname
thedir+='/Year_'+str(season)+'/N_healpix_64'
files=glob.glob(thedir+'/*.txt')
from descartes.patch import PolygonPatch
#print files

r=[]
for fi in files:
    if fi.find('orig') == -1 and fi.find('ref') == -1:
        #print fi
        fieldid=int(fi.split('_')[-1].split('.')[0])
        myobs=Observations(fieldid=fieldid, filename=fi)
        if len(myobs.seasons) >= 1:
            season=myobs.seasons[0]
            r.append((fieldid,len(season)))
            #print fieldid,season['Ra'][0],season['Dec'][0],len(season)
            if len(season) > 150.:
                vo=[]
                figb,axb = plt.subplots(ncols=1, nrows=1) 
                for val in season[['Ra','Dec']]:
                    vo.append(val)
                poly=geometry.Polygon(vo)
                x, y = poly.exterior.coords.xy
                print('oo',x,y)
                r=[]
                for (xa,ya) in zip(x,y):
                    r.append([xa,ya])
                poly_ext=geometry.Polygon(np.array(r)).convex_hull
                print('area',poly_ext.area)
                patch = PolygonPatch(poly_ext, facecolor='#6699cc', edgecolor='k')
                axb.add_patch(patch)
                axb.plot(season['Ra'],season['Dec'],'ko')
                plt.show()
            
    
res=np.rec.fromrecords(r,names=['fieldid','Nvisits'])
plt.hist(res['Nvisits'],histtype='step')
