import numpy as np
from LC_Ana import *
#from optparse import OptionParser
from ID import *
import cPickle as pkl
import matplotlib.pyplot as plt
from scipy.spatial import distance
from Calc_NSN import *
import operator

def get_dir(stra,strb,strc):
    return stra+strb+strc

def get_operator_fn(op):
    return {
        '+' : operator.add,
        '-' : operator.sub,
        '*' : operator.mul,
        '/' : operator.div,
        '%' : operator.mod,
        '^' : operator.xor,
        '>=': operator.ge,
        '>' : operator.gt,
        '<=': operator.le,
        '<' : operator.lt,
        '==': operator.eq,
        }[op]

def Load_Selection(selname):

    test=np.loadtxt(selname,dtype={'names': ('cut','comp','val','type'),'formats': ('S150','S2','S8','S5')})

    print(test)

    sel=[]
    for cut in test:
        thecut=[]
        for val in cut['cut'].split('+'):
            thecut.append(val)
        #print(cut['type'])
        if cut['type'] != 'str':
            sel.append((thecut,get_operator_fn(cut['comp']),eval(cut['type']+'('+cut['val']+')')))
        else:
            sel.append((thecut,get_operator_fn(cut['comp']),cut['val']))
            
    return sel

def Get_Efficiencies(fieldname,fieldid,zrange,x1_c,main_dir,data_dir,selection,simu_name,preselection):
    
    for season in range(10):
        take_label=False
        effi=None
        for x1_c in x1_color:
            proc=LC_Ana(fieldname,fieldid,season,zrange,x1_c,main_dir=main_dir,data_dir=data_dir,selection=selection,simu_name=simul_name)
            proc.Effi(preselection)

def Save_Selected(fieldname,fieldid,main_dir,data_dir,selection,simu_name):
    
    for season in range(10):
        proc=LC_Ana(fieldname,fieldid,season,zrange=[],x1_color=[],main_dir=main_dir,data_dir=data_dir,selection=selection,simu_name=simul_name)
        proc.Save_Selected()

def Plot_Efficiency_Sel(axes,fieldname,fieldid,zrange,x1_c,main_dir,data_dir,selection,simu_name,preselection):

    tot_label=[]
    for season in range(10):
        take_label=False
        effi=None
        for x1_c in x1_color:
            if x1_c == (0.0,0.0):
                take_label=True
            proc=LC_Ana(fieldname,fieldid,season,zrange,x1_c,main_dir=main_dir,data_dir=data_dir,selection=selection,simu_name=simul_name)
            proc.Effi(preselection)
            proc.Plot_Effi_new(ls=myls[x1_c],color=colors[season],tot_label=tot_label,take_label=take_label,ax=axes,selection_name=selection.keys()[0])
            
            if effi is None:
                effi=proc.efficiencies[selection.keys()[0]]
            else:
                effi=np.concatenate((effi,proc.efficiencies[selection.keys()[0]]))
          
        labs = [l.get_label() for l in tot_label]
        axes.legend(tot_label, labs, ncol=1,loc='best',prop={'size':fontsize},frameon=False)
        axes.set_xlabel('z',{'fontsize': fontsize})
        axes.set_ylabel('Efficiency',{'fontsize': fontsize})
        

def Plot_Efficiencies(fieldname,fieldid,selection,main_dir,data_dir,preselection):
    #selection_names=['Sela_no_color_cut','Sela_with_color_cut','Science_Book_no_colorcut','Science_Book_with_colorcut']
    #selection_names=['Sela_no_color_cut','Sela_with_color_cut']
    #selection_names=['Science_Book_no_colorcut','Selb_10_2']
    for selection_name in selection.keys():
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(10,9))
        fig.suptitle(fieldname +' '+str(fieldid)+' - '+simul_name+ '\n'+selection_name)
        sel_dict={}
        sel_dict[selection_name]=selection[selection_name]
        Plot_Efficiency_Sel(axes,fieldname,fieldid,zrange,x1_color,main_dir=main_dir,data_dir=data_dir,selection=sel_dict,simu_name=simul_name,preselection=preselection)
        fig.savefig('Plot_Effi/effi_'+fieldname+'_'+str(fieldid)+'_'+selection_name+'.png')

def Plot_zpercentile(fieldname,fieldid,selection,main_dir,data_dir):
    
    for selection_name in selection.keys():
        sel_dict={}
        sel_dict[selection_name]=selection[selection_name]
        tab=Get_zpercentile(fieldname,fieldid,zrange,x1_color,main_dir=main_dir,data_dir=data_dir,selection=sel_dict,simu_name=simul_name)

        plt.plot(tab['season'],tab['z_90'])
        plt.plot(tab['season'],tab['z_95'])
        plt.show()


def Get_zpercentile(fieldname,fieldid,zrange,x1_c,main_dir,data_dir,selection,simu_name):
     zmax=1.225

     r=[]
     for season in range(10):
         effi=None
         for x1_c in x1_color:
             proc=LC_Ana(fieldname,fieldid,season,zrange,x1_c,main_dir=main_dir,data_dir=data_dir,selection=selection,simu_name=simul_name)
             proc.Effi()
             if effi is None:
                 effi=proc.efficiencies[selection.keys()[0]]
             else:
                 effi=np.concatenate((effi,proc.efficiencies[selection.keys()[0]]))

         idx= effi['z']<zmax
         effi=effi[idx]
         nsn=Calc_NSN(fieldname,fieldid,season,effi,simu_name,selection.keys()[0])
         z_90,z_95=nsn.zpercentile()
         r.append((fieldname,fieldid,season,selection.keys()[0],simul_name,z_90,z_95))
     return np.rec.fromrecords(r,names=['fieldname','fieldid','season','selection_name','simul_name','z_90','z_95'])
         


def Plot_N_SN(fieldname,fieldid,selection,main_dir,data_dir,preselection):
    
    for selection_name in selection.keys():
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(10,9))
        fig.suptitle(fieldname +' '+str(fieldid)+' - '+simul_name+ '\n'+selection_name)
        sel_dict={}
        sel_dict[selection_name]=selection[selection_name]
        Plot_N_SN_Sel(axes,fieldname,fieldid,zrange,x1_color,main_dir=main_dir,data_dir=data_dir,selection=sel_dict,simu_name=simul_name,preselection=preselection)
        fig.savefig('Plots_NSN/NSN_'+fieldname+'_'+str(fieldid)+'_'+selection_name+'.png')

def Plot_N_SN_Sel(axes,fieldname,fieldid,zrange,x1_c,main_dir,data_dir,selection,simu_name,preselection):
     zmax=1.225
     cumul=True
     tot_label=[]
     for season in range(10):
         effi=None
         for x1_c in x1_color:
             proc=LC_Ana(fieldname,fieldid,season,zrange,x1_c,main_dir=main_dir,data_dir=data_dir,selection=selection,simu_name=simul_name)
             proc.Effi(preselection)
             if effi is None:
                 effi=proc.efficiencies[selection.keys()[0]]
             else:
                 effi=np.concatenate((effi,proc.efficiencies[selection.keys()[0]]))

         idx= effi['z']<zmax
         effi=effi[idx]
         print(effi)
         nsn=Calc_NSN(fieldname,fieldid,season,effi,simu_name,selection.keys()[0])
         nsn.Plot_N_SN(axes,cumul=cumul,tot_label=tot_label,marker=markers[season],color=color_seas[season])
         labs = [l.get_label() for l in tot_label]
         axes.legend(tot_label, labs, ncol=1,loc='best',prop={'size':fontsize},frameon=False)
         axes.set_xlabel('z',{'fontsize': fontsize})
         axes.set_xlim([0.,nsn.zmax_eff])
         if cumul:
             axes.set_ylabel('Number of SN Ia < z',{'fontsize': fontsize})
      
def Show_Dist_Error(fieldname,fieldid,season=0):

    colors=dict(zip(['Sela_cc','Selb_10_2_cc'],['k','r']))
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(10,9))
    for selection_name in ['Sela_cc','Selb_10_2_cc']:
        
        fig.suptitle(fieldname +' '+str(fieldid)+' - '+simul_name+ '\n'+selection_name)
        proc=LC_Ana(fieldname,fieldid,season,zrange,x1_color=(0.0,0.0),main_dir=main_dir,data_dir=data_dir+'_covmb',selection_name=selection_name,selection=selection,simu_name=simul_name,add_covmb=True)
        proc.Load_Dist_Error()
        proc.Plot_Var(axes,'sigma_dist_modulus',color=colors[selection_name])
    plt.show()


min_rf_phase=-20.0
max_rf_phase=60.0

#zmax=0.7
zmax=1.5
zrange=[0.01]
step=0.025
zrange+=[val for val in np.arange(step,zmax,step)]
#z=[0.1,0.6]
print 'hello z',zrange
#z=[0.375]

add='_0_5'
add+='_'+str(min_rf_phase).replace('-','m')+'_'+str(max_rf_phase)

simul_name='sncosmo'
#selection_names=['Sela','Sela_cc','Selb_10_2','Selb_10_2_cc','Selb_10_3','Selb_10_3_cc','Selb_20_2','Selb_20_2_cc','Selb_20_3','Selb_20_3_cc']
selection_names=['Sela','Selb_10_2','Sela_cc','Selb_10_2_cc','Selb_10_3','Selb_10_3_cc']
selection_names=['Sela_cc_tot_oldcut']
#selection_names=['Selb_20_2','Selb_20_2_cc','Selb_20_3','Selb_20_3_cc']
selection={}

for selname in selection_names:
    selection[selname]=Load_Selection(selname+'.txt')
    
preselection=Load_Selection('Presel.txt')

x1_color=[(-2.0,0.2),(2.0,-0.2),(0.0,0.0)]
add='_0_5'
add='_random'
add_all='_'+str(min_rf_phase).replace('-','m')+'_'+str(max_rf_phase)
data_dir='Fitted_Light_Curves_'+simul_name
main_dir='/sps/lsst/data/dev/pgris'


myls=dict(zip(x1_color,[':','-.','-']))
colors=dict(zip(range(10),['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']))


#selection_name='Sela_no_color_cut'


fontsize=12
markers=['o','o','s','s','.','.','^','^','<','<']
color_seas=['b','r','b','r','b','r','b','r','b','r']

#to process efficiencies
"""
for fieldid in [290,744,1427,2412,2786]:
#for fieldid in [100,101,102,103]:
    Get_Efficiencies('DD',fieldid,zrange,x1_color,main_dir,get_dir(data_dir,'_0_5',add_all),selection,simul_name,preselection)
"""
#to select data from random
"""
for fieldid in [744]:
    #Save_Selected('DD',fieldid,main_dir,data_dir+'_covmb',selection,simul_name)
    data_dir_new=get_dir(data_dir,'_random',add_all)
    Save_Selected('DD',fieldid,main_dir,data_dir_new+'_covmb',selection,simul_name)
"""



#to plot efficiencies


#for fieldid in [744,1427,2412,2786]:
#for fieldid in [100,101,102,103]:
for fieldid in [744]:
    Plot_Efficiencies('DD',fieldid,selection,main_dir,get_dir(data_dir,'_0_5',add_all),preselection=None)
    Plot_N_SN('DD',fieldid,selection,main_dir,get_dir(data_dir,'_0_5',add_all),preselection=None)
    #Plot_zpercentile('DD',fieldid,selection,main_dir,get_dir(data_dir,'_0_5',add_all))
    #plt.show()



"""
for fieldid in [744]:
    Show_Dist_Error('DD',fieldid,season=0)
"""   
