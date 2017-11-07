import numpy as np
import os
import glob
import cPickle as pkl

class LC_Ana:
    def __init__(self,data_dict={},zmin=0.,zmax=1.2,bin_z=0.01,data_dir='../Fitted_Light_Curves'):

        self.data={}
        self.data_dir=data_dir
        self.Load_Data(data_dict)

        self.data_sel={}

        for key, vals in self.data.items():
            print key
            self.data_sel[key]=self.Select_Data(vals)
            print len(vals),len(self.data_sel[key])
        


    def Load_Data(self, data_dict={}):

        for key,val in data_dict.items():
            dirmeas=self.data_dir+'/'+val.thedir+'/'+str(val.fieldid)+'/Season_'+str(val.season)
            sum_file=dirmeas+'/'+val.fieldname+'_'+str(val.fieldid)+'_X1_'+str(val.X1)+'_C_'+str(val.Color)+'_all.pkl'
            
            if os.path.exists(sum_file):
                pkl_file = open(sum_file,'rb')
                print 'loading',sum_file
                loaded_file=pkl.load(pkl_file)
                #print 'done',key,tot_resu[key],tot_resu[key].dtype,sum_file
                #idx = loaded_file['z']>=self.zmin and loaded_file['z']<self.zmax
                self.data[key]=loaded_file
            else:
                files = glob.glob(dirmeas+'/'+val.fieldname+'_'+str(val.fieldid)+'*_X1_'+str(val.X1)+'_C_'+str(val.Color)+'*.pkl')
                for fi in files:
                    pkl_file = open(fi,'rb')
                    print 'loading',fi
                    if not key in self.data.keys():
                        self.data[key]=pkl.load(pkl_file)
                    else:
                        self.data[key]=np.vstack((self.data[key],pkl.load(pkl_file)))

                pkl_out = open(sum_file,'wb')
                
                pkl.dump(self.data[key], pkl_out)
                
                pkl_out.close()


    def Select_Data(self,data):

        sel=data.copy()
        logical_test={}
    
        print 'before',len(sel)
        for band in 'grizy':
            logical_test[band]=np.logical_and(sel['N_bef_'+band]>=1,sel['N_aft_'+band]>=1)

        logical_and_g_r=np.logical_and(logical_test['g'],logical_test['r'])
        logical_and_r_i=np.logical_and(logical_test['r'],logical_test['i'])
        logical_and_i_z=np.logical_and(logical_test['i'],logical_test['z'])
        logical_and_z_y=np.logical_and(logical_test['z'],logical_test['y'])
                    
        print 'logical',len(sel)
        sel=sel[np.where(np.logical_and(sel['status']=='go_fit',sel['fit_status']=='fit_ok'))]
        for val in sel:
            print val['phase_first'],(val['DayMax']-val['salt2.T0'])/(1.+val['z'])
        sel=sel[np.where(np.logical_and(sel['phase_first']<=-5,sel['phase_last']>=20))]
        
        print 'all',len(sel)
        return sel


