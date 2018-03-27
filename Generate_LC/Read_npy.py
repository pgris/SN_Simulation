import numpy as np

tab=np.load('prod_lc_sncosmo_random.npy')

print(len(tab))

idx = tab['season'] == 0


test=np.unique(tab[idx][['z','x1','c','DayMax']])
print(len(tab[idx]),len(test),test['DayMax'])
