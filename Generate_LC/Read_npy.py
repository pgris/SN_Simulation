import numpy as np
import matplotlib.pyplot as plt

tab=np.load('prod_lc_sncosmo_random_744_test.npy')

print(len(tab))

idx = tab['season'] == 0


test=np.unique(tab[idx][['z','x1','c','DayMax']])
print(len(tab[idx]),len(test),test['DayMax'],tab['x1'])

plt.hist(tab['c'],range=[-6.,6.])
plt.show()
