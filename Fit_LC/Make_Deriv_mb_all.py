from scipy.interpolate import griddata
import glob
import numpy as np
import cPickle as pkl


atype='high_z'
files=glob.glob('Files_Deriv_mb/res_'+str(atype)+'*.pkl')

tot_pars=None
print 'alors',len(files)

for fi in files:
    list_pars=pkl.load(open(fi,'rb'))
    print 'loading',fi
    if tot_pars is None:
        #print list_pars.dtype
        tot_pars=list_pars
    else:
        #print tot_pars.dtype,len(tot_pars)
        #print list_pars.dtype,len(list_pars)
        tot_pars=np.concatenate((tot_pars,list_pars))

pkl_file = open('Files_Deriv_mb/Derib_mb_'+atype+'_all.pkl','wb')
pkl.dump(tot_pars, pkl_file)       
pkl_file.close()

print 'hello',len(tot_pars),tot_pars[['z','X0','X1','Color']][:10]


m_z=np.mean(tot_pars['z'])
m_x1=np.mean(tot_pars['X1'])
m_c=np.mean(tot_pars['Color'])
m_x0=np.mean(tot_pars['X0'])

z0 = griddata((tot_pars['z'],tot_pars['X0'],tot_pars['X1'],tot_pars['Color']), tot_pars['dmb_dx0'], (m_z,m_x0,m_x1,m_c), method='nearest')


print 'hello',z0,np.mean(tot_pars['dmb_dx0'])
