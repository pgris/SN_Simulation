import cPickle as pkl
import matplotlib.pyplot as plt
from scipy.spatial import distance

fich='DD_744_0.1_X1_-999.0_C_-999.0_900_901.pkl'

dira='../../Light_Curves_sncosmo_test_multi/DD/744/Season_0/'
dirb='../../Light_Curves_sncosmo_test/DD/744/Season_0/'

lca=pkl.load(open(dira+fich,'rb'))
lcb=pkl.load(open(dirb+fich,'rb'))

print len(lca),len(lcb)
print lca[1207].meta
print lcb[1207].meta


for i in range(len(lca)):
    lc_cos=lca[i]

    points=[(lcb[i].meta['DayMax'],lcb[i].meta['X1'],lcb[i].meta['Color']) for i in range(len(lcb))]
    
    idx_nearest=distance.cdist([(lc_cos.meta['DayMax'],lc_cos.meta['X1'],lc_cos.meta['Color'])],points).argmin()
    
    lc_sim=lcb[idx_nearest]

    print len(lc_cos),len(lc_sim)
    print lc_cos
    print lc_sim
    break


"""
figb, axb = plt.subplots(ncols=2, nrows=2, figsize=(10,9))
axb[0][0].hist([lca[i].meta['X1'] for i in range(len(lca))])
axb[0][1].hist([lca[i].meta['Color'] for i in range(len(lca))])
axb[1][0].hist([lcb[i].meta['X1'] for i in range(len(lcb))])
axb[1][1].hist([lcb[i].meta['Color'] for i in range(len(lcb))])

plt.show()
"""
