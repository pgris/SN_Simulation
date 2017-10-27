import numpy as np
from tempfile import TemporaryFile

def Save_npz(name,x1_mean,sig_m_x1,sig_p_x1,c_mean,sig_m_c,sig_p_c):

    x1_vals, x1_weights = gauss_asym_distrib(x1_mean,sig_m_x1,sig_p_x1)
    c_vals, c_weights = gauss_asym_distrib(c_mean,sig_m_c,sig_p_c)
    
    fichname='Dist_X1_Color_'+name
    np.savez(fichname,x1_vals=x1_vals, x1_weights=x1_weights,c_vals=c_vals, c_weights=c_weights)

def gauss_asym_distrib(mean,sigma_minus,sigma_plus):
    xmin=mean-5.*sigma_minus
    xmax=mean+5.*sigma_plus

    pas=1.e-4

    nsteps=int((xmax-xmin)/pas)
    
    xvals=[]
    weights=[]
    
    for i in range(nsteps):
        x=xmin+float(i)*pas
        if x < mean:
            res=np.exp(-np.power(x-mean,2)/(2*np.power(sigma_minus,2.)))
        else:
            res=np.exp(-np.power(x-mean,2)/(2*np.power(sigma_plus,2.)))
            
        xvals.append(x)
        weights.append(res)


    return xvals,weights/np.sum(weights)


#Load x1 and c asymmetric distributions
# values from Scolnic & Kessler may 2016 arXiv:1603.01559v2

x1_mean=0.964
sig_m_x1=1.467
sig_p_x1=0.235
c_mean=-0.099
sig_m_c=0.003
sig_p_c=0.119

Save_npz('high_z',x1_mean=0.964,sig_m_x1=1.467,sig_p_x1=0.235,c_mean=-0.099,sig_m_c=0.003,sig_p_c=0.119)
Save_npz('low_z',x1_mean=0.419,sig_m_x1=3.024,sig_p_x1=0.742,c_mean=-0.069,sig_m_c=0.003,sig_p_c=0.148)



    #npzfile = np.load('test.npz','r')
    #print npzfile['x1_vals']
