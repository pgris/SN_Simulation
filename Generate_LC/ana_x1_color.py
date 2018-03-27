import numpy as np


def Ana_File(filename):
    res=np.loadtxt(filename, dtype={'names': ('x1', 'c', 'weight_x1','weight_c','weight_tot'),'formats': ('f8', 'f8', 'f8','f8','f8')})
    #res.sort(order='weight_tot')
    #res[:]=res[::-1]
   
    print res
    """
    sum_w=0
    for i,val in enumerate(res):
        sum_w+=val['weight_tot']
        print i,val,sum_w
        if sum_w > 0.99:
            break

    print len(res)
    """

filea='Dist_X1_Color_low_z.txt'

Ana_File(filea)
