import cPickle as pkl
import matplotlib.pyplot as plt
import numpy as np

def Histo_ratio(sela,selb,varname,zmin,zmax,bin_z):

        range=[zmin,zmax]
        bin_width=bin_z
        #print 'alors pal',range,bin_width
        num_bins=int((range[1]-range[0])/bin_width)
        #range=[0.0,1.2]
        
        #hista, bin_edgesa = np.histogram(sela[varname],bins=num_bins,range=range,weights=sela['X1_weight']*sela['Color_weight'])
        #histb, bin_edgesb = np.histogram(selb[varname],bins=num_bins,range=range,weights=selb['X1_weight']*selb['Color_weight'])
        hista, bin_edgesa = np.histogram(sela[varname],bins=num_bins,range=range)
        histb, bin_edgesb = np.histogram(selb[varname],bins=num_bins,range=range)
        #print hista
        #print histb
       
        #bin_center = (bin_edgesa[:-1] + bin_edgesa[1:]) / 2

        
        ratio=[]
        ratio_err=[]
        norm=[]
        norm_err=[]
        
        for a,b in zip(histb,hista):
            if b==0:
                ratio.append(a)
                ratio_err.append(0)
            else:
                #print 'hello',a,b
                effi=float(a)/float(b)
                ratio.append(effi)
                ratio_err.append(np.sqrt(1.-effi)*effi/np.sqrt(float(b)))
            eff=float(a)/float(np.sum(hista))
            norm.append(eff)
            norm_err.append(np.sqrt(1.-eff)*eff/np.sqrt(float(np.sum(hista))))
    
        return bin_edgesa[:-1], ratio, ratio_err,norm,norm_err

dirsel='Sel_Selb_10_2'
fname='DD_744_Season_0.pkl'

tab=pkl.load(open(dirsel+'/Sel_'+fname,'rb'))
tab_all=pkl.load(open(dirsel+'/Tot_'+fname,'rb'))
print(len(tab),len(tab_all))

"""
plt.hist(tab['z'],histtype='step',range=[0.,1.5],bins=30)
plt.hist(tab_all['z'],histtype='step',range=[0.,1.5],bins=30)
"""
"""
edges, ratio, ratio_err,norm,norm_err=Histo_ratio(tab_all,tab,'z',0.,1.5,0.1)

plt.errorbar(edges,ratio,yerr=ratio_err)
"""

names=['name','zcmb','zhel','dz','mb','dmb','x1','dx1','color','dcolor','3rdvar','d3rdvar','tmax','dtmax','cov_m_s','cov_m_c','cov_s_c','set','ra','dec','biascor']

jla_tab=np.loadtxt('jla_lcparams.txt',dtype={'names': tuple(names),'formats': tuple(['S15']+[np.float]*(len(names)-1))})
print(jla_tab)
figb, axb = plt.subplots(2, 2, figsize=(10,9))
idx = tab_all['z']<0.1
idxb = jla_tab['zcmb']<0.1
axb[0][0].hist(tab_all[idx]['X1'],histtype='step',color='k')
axb[0][0].hist(jla_tab[idxb]['x1'],histtype='step',color='r')
axb[0][1].hist(tab_all[idx]['Color'],histtype='step',color='k')
axb[0][1].hist(jla_tab[idxb]['color'],histtype='step',color='r')
idxa = tab_all['z']>=0.1
idxab = jla_tab['zcmb']>=0.1
axb[1][0].hist(tab_all[idxa]['X1'],histtype='step',color='k')
axb[1][0].hist(jla_tab[idxab]['x1'],histtype='step',color='r')
axb[1][1].hist(tab_all[idxa]['Color'],histtype='step',color='k')
axb[1][1].hist(jla_tab[idxab]['color'],histtype='step',color='r')
plt.show()
