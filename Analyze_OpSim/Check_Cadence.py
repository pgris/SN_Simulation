import numpy as np
import pickle as pkl

fichname='Cadence_WFD_all_Season_2.pkl'
tab=pkl.load(open(fichname,'rb'))

print(tab.dtype)
print(len(tab['Nvisits_g'][tab['Nvisits_g']>00.]),tab['mjd_z'])
