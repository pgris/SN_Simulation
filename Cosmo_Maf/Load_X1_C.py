import numpy as np


class Load_X1_Color:
    def __init__(self,stretch,color,percent=0.95,name_low_z='Dist_X1_Color_low_z.txt',name_high_z='Dist_X1_Color_high_z.txt'):

        #print 'alors man',stretch,color
        self.name_low_z=name_low_z
        self.name_high_z=name_high_z
        self.tab=self.Get_X1_Color_Distribution(stretch,color,percent)
        

    def Load_File(self,filename):
        res=np.loadtxt(filename, dtype={'names': ('x1', 'c', 'weight_x1','weight_c','weight_tot'),'formats': ('f8', 'f8', 'f8','f8','f8')})
        res.sort(order='weight_tot')
        res[:]=res[::-1]
    
        return res

    def Get_X1_Color_Distribution(self,stretch,color,percent):
    
        Dist_X1_Color={}

        Dist_X1_Color['low_z']=self.Load_File(self.name_low_z)
        Dist_X1_Color['high_z']=self.Load_File(self.name_high_z)

        #print 'loadad',len(Dist_X1_Color['low_z'])
        res={}
        if stretch > -90. and color > -90.:

            for val in ['low_z','high_z']:
                res[val]=self.Select(Dist_X1_Color[val],stretch,color)
        else:
            for val in ['low_z','high_z']:
                res[val]=self.Select_percent(Dist_X1_Color[val],percent)
                #print 'there',val,len(res[val])
        return res

    def Select(self,ar,stretch,color):

        #print ar
        #print 'looking for',stretch,color
        idx = (ar['x1']==stretch)&(ar['c']==color)
        sel=ar[idx]
    
        if len(sel)==0:
            print 'problem',stretch,color,'not found - Weights set to 1'
       
            return np.asarray((stretch,color,1.,1.,1.),dtype=ar.dtype)

        else:
            return sel

    def Select_percent(self,dist,per):

        sum_w=0
        res=None
        for i,val in enumerate(dist):
            sum_w+=val['weight_tot']
        #print i,val,sum_w
            if sum_w <= per:
                if res is None:
                    res=np.array(val, dtype=val.dtype)
                else:
                    res=np.vstack([res,np.array(val, dtype=val.dtype)])
            else:
                return res

        return res
