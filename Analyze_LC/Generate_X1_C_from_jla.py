import numpy as np
import matplotlib.pyplot as plt

def Get_Distrib(sel,varname,bin_width):
    r=[]
    min_val=round(np.min(sel[varname]),2)
    max_val=round(np.max(sel[varname]),2)
    """
    min_val=np.min(sel[varname])
    max_val=np.max(sel[varname])
    """
    print('hello',min_val,max_val)
    num_bins=int((max_val-min_val)/bin_width)
    hista, bin_edgesa = np.histogram(sel[varname],bins=num_bins+1,range=[min_val-bin_width/2.,max_val+bin_width/2.])
    
    bin_center = (bin_edgesa[:-1] + bin_edgesa[1:]) / 2

    print(min_val,max_val)
    for i,val in enumerate(hista):
        print(round(bin_center[i],3),float(hista[i])/float(np.sum(hista)))
        r.append((round(bin_center[i],3),float(hista[i])/float(np.sum(hista))))

    """
    plt.hist(sel[varname],histtype='step',bins=num_bins+1,range=[min_val-bin_width/2.,max_val+bin_width/2.])
    #plt.hist(sel[varname],histtype='step',bins='auto')
    plt.show()
    """
    return r

def write_in_file(dist,name):
    ftxt = open('Dist_X1_Color_jla_'+name+'.txt',"w")
    ftxt.write('# x1 c weight_x1 weight_c weight_tot\n')
    for x1 in dist['x1']:
        if x1[1] > 0.:
            for color in dist['color']:
                if color[1] > 0. :
                    print(x1[0],color[0],x1[1],color[1],x1[1]*color[1])
                    ftxt.write(str(x1[0])+' '+str(color[0])+' '+str(x1[1])+' '+str(color[1])+' '+str(x1[1]*color[1])+ '\n')
    ftxt.close()


names=['name','zcmb','zhel','dz','mb','dmb','x1','dx1','color','dcolor','3rdvar','d3rdvar','tmax','dtmax','cov_m_s','cov_m_c','cov_s_c','set','ra','dec','biascor']

jla_tab=np.loadtxt('jla_lcparams.txt',dtype={'names': tuple(names),'formats': tuple(['S15']+[np.float]*(len(names)-1))})

#print(jla_tab)

idx = jla_tab['zcmb']<0.1
#idxb = jla_tab['zcmb']>=0.1

sel=jla_tab[idx]
bin_width=dict(zip(['x1','color'],[0.1,0.02]))

dist={}
for varname in ['x1','color']:
    dist[varname]=Get_Distrib(sel,varname,bin_width[varname])
write_in_file(dist,'low_z')

 
