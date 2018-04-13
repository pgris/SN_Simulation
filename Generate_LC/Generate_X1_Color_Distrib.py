import numpy as np
from tempfile import TemporaryFile
import matplotlib.pyplot as plt
import cPickle as pkl

def Save_npz(name,x1_mean,sig_m_x1,sig_p_x1,c_mean,sig_m_c,sig_p_c,min_x1,max_x1,min_c,max_c,x1_step,c_step):

    x1_vals, x1_weights = gauss_asym_distrib(x1_mean,sig_m_x1,sig_p_x1,x1_step,min_x1,max_x1+x1_step)
    c_vals, c_weights = gauss_asym_distrib(c_mean,sig_m_c,sig_p_c,c_step,min_c,max_c+c_step)
    
    #print x1_vals, x1_weights

    fichname='Dist_X1_Color_'+name
    np.savez(fichname,x1_vals=x1_vals, x1_weights=x1_weights,c_vals=c_vals, c_weights=c_weights)

def gauss_asym_distrib(mean,sigma_minus,sigma_plus,pas,min_val,max_val):
    xmin=mean-5.*sigma_minus
    xmax=mean+5.*sigma_plus

    xmin=min_val
    xmax=max_val
    #pas=0.1

    nsteps=int((xmax-xmin)/pas)
    
    xvals=[]
    weights=[]
    
    for i in range(nsteps):
        x=xmin+float(i)*pas
        if x < mean:
            res=np.exp(-np.power(x-mean,2)/(2*np.power(sigma_minus,2.)))
        else:
            res=np.exp(-np.power(x-mean,2)/(2*np.power(sigma_plus,2.)))
        #if res>0.0000001:
        xvals.append(x)
        weights.append(res)

    print 'alors',xvals
    return xvals,weights/np.sum(weights)

def Gime_Numbers(filename):

    X1_Color_npzfile = np.load(filename,'r')

    n_x1=len(X1_Color_npzfile['x1_vals'])
    n_color=len(X1_Color_npzfile['c_vals']) 

    thestr='N_x1 = '+str(n_x1)+' N_col = '+str(n_color)
    return thestr

def Show_Results(filename):

    X1_Color_npzfile = np.load(filename,'r')
    
    print 'alors',filename
    ftxt = open(filename.replace('.npz','.txt'),"w")
    ftxt.write('# x1 c weight_x1 weight_c weight_tot\n')
    for i in range(len(X1_Color_npzfile['x1_vals'])):
        x1=X1_Color_npzfile['x1_vals'][i]
        weight_x1=X1_Color_npzfile['x1_weights'][i]
        for j in range(len(X1_Color_npzfile['c_vals'])):
            c=X1_Color_npzfile['c_vals'][j]
            weight_c=X1_Color_npzfile['c_weights'][j]
            #print x1,c,weight_x1,weight_c,weight_x1*weight_c
            ftxt.write(str(x1)+' '+str(c)+' '+str(weight_x1)+' '+str(weight_c)+' '+str(weight_x1*weight_c)+ '\n')
    ftxt.close()

    thetype='low_z '
    if 'high_z' in filename:
        thetype='high_z'

    #print X1_Color_npzfile['x1_vals'],X1_Color_npzfile['x1_weights'],np.sum(X1_Color_npzfile['x1_weights']),len(X1_Color_npzfile['x1_vals'])
    #print X1_Color_npzfile['c_vals'],X1_Color_npzfile['c_weights'],np.sum(X1_Color_npzfile['c_weights']),len(X1_Color_npzfile['c_vals'])
    
    X1_val=[]
    Color_val=[]
    
    for j in range(0,400):
        X1=np.random.choice(X1_Color_npzfile['x1_vals'],1,p=X1_Color_npzfile['x1_weights'])[0]
        indx=X1_Color_npzfile['x1_vals'] == X1
        X1_val.append((X1_Color_npzfile['x1_vals'][indx][0],X1_Color_npzfile['x1_weights'][indx][0]))
    for j in range(0,200):
        Color=np.random.choice(X1_Color_npzfile['c_vals'],1,p=X1_Color_npzfile['c_weights'])[0]
        indx=X1_Color_npzfile['c_vals'] == Color
        Color_val.append((X1_Color_npzfile['c_vals'][indx][0],X1_Color_npzfile['c_weights'][indx][0]))

    
    list_X1=list(set(X1_val))
    list_Color=list(set(Color_val))


    print 'totala',np.sum([list_X1[i][1] for i in range(len(list_X1))]),np.sum([list_Color[i][1] for i in range(len(list_Color))]),len(list_X1),len(list_Color)

    combin=[]
    for val in list_X1:
        for valb in list_Color:
            combin.append((val,valb))

    combin_tot=[]
    idxa=(X1_Color_npzfile['x1_vals']>-100.)&(X1_Color_npzfile['x1_vals']<100.)
    idxb=(X1_Color_npzfile['c_vals']>-100)&(X1_Color_npzfile['c_vals']<100.)

    r=[]
    for val,valw in zip(X1_Color_npzfile['x1_vals'],X1_Color_npzfile['x1_weights']):
        for valb,valbw in zip(X1_Color_npzfile['c_vals'],X1_Color_npzfile['c_weights']):
            r.append((val,valw,valb,valbw,thetype))

    restab=np.rec.fromrecords(r,names=['X1','X1_weight','Color','Color_weight','type'])
    #restab=np.rec.fromarrays(r,names=['X1','X1_weight','Color','Color_weight','type'])
    print 'totalb',len(restab)

    figb, axb = plt.subplots(ncols=1, nrows=2, figsize=(10,9))
    axb[0].hist(restab['X1'],weights=restab['X1_weight'],normed=True,bins=30)

    axb[1].hist(restab['Color'],weights=restab['Color_weight'],normed=False,bins=30)


    return restab


#Load x1 and c asymmetric distributions
# values from Scolnic & Kessler may 2016 arXiv:1603.01559v2

x1_mean=0.964
sig_m_x1=1.467
sig_p_x1=0.235
c_mean=-0.099
sig_m_c=0.003
sig_p_c=0.119

#Save_npz('high_z',x1_mean=0.964,sig_m_x1=1.467,sig_p_x1=0.235,c_mean=-0.099,sig_m_c=0.003,sig_p_c=0.119,min_x1=-2.0,max_x1=2.0,min_c=-0.2,max_c=0.2,x1_step=0.1,c_step=0.01)
#Save_npz('low_z',x1_mean=0.419,sig_m_x1=3.024,sig_p_x1=0.742,c_mean=-0.069,sig_m_c=0.003,sig_p_c=0.148,min_x1=-6.0,max_x1=2.0,min_c=-0.2,max_c=0.3,x1_step=0.2,c_step=0.01)

Save_npz('high_z',x1_mean=0.964,sig_m_x1=1.467,sig_p_x1=0.235,c_mean=-0.099,sig_m_c=0.003,sig_p_c=0.119,min_x1=-2.0,max_x1=2.0,min_c=-0.2,max_c=0.2,x1_step=0.2,c_step=0.02)
Save_npz('low_z',x1_mean=0.419,sig_m_x1=3.024,sig_p_x1=0.742,c_mean=-0.069,sig_m_c=0.003,sig_p_c=0.148,min_x1=-6.0,max_x1=2.0,min_c=-0.2,max_c=0.3,x1_step=0.4,c_step=0.02)

print 'low z',Gime_Numbers('Dist_X1_Color_low_z.npz')
print 'high z',Gime_Numbers('Dist_X1_Color_high_z.npz')


tab_low=Show_Results('Dist_X1_Color_low_z.npz')
tab_high=Show_Results('Dist_X1_Color_high_z.npz')

"""
tab_low=np.reshape(tab_low,(len(tab_low),1))
tab_high=np.reshape(tab_high,(len(tab_high),1))
print tab_low.shape,tab_high.shape
"""
tab_tot=np.concatenate((tab_low,tab_high))
pkl_file = open('Map_X1_C_b.pkl','wb')
pkl.dump(tab_tot, pkl_file)
pkl_file.close()

#print tab_tot

plt.show()
#print set(X1_val),set(Color_val)
