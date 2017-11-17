import numpy as np
from LC_Ana import *
from optparse import OptionParser
from ID import *
import cPickle as pkl
import matplotlib.pyplot as plt
from scipy.spatial import distance

parser = OptionParser()

parser.add_option("-f", "--fieldname", type="string", default='DD', help="filter [%default]")
parser.add_option("-i", "--fieldid", type="int", default=744, help="filter [%default]")
parser.add_option("-x", "--stretch", type="float", default=-999., help="filter [%default]")
parser.add_option("-c", "--color", type="float", default=-999., help="filter [%default]")
parser.add_option("-d", "--dirmeas", type="string", default="DD", help="filter [%default]")
parser.add_option("--dirobs", type="string", default="DD", help="filter [%default]")

opts, args = parser.parse_args()

dict_data={}

fieldid=opts.fieldid
fieldname=opts.fieldname
thedir=opts.dirmeas
thedir_obs=opts.dirobs
X1=opts.stretch
Color=opts.color

for seas in [0]:
    dict_data['DD_'+str(fieldid)+'_'+str(seas+1)]=Id(thedir=thedir,thedir_obs=thedir_obs,fieldname=fieldname,fieldid=fieldid,X1=X1,Color=Color,season=seas,colorfig='k')

LC_Ana(dict_data,zmin=0.01,zmax=1.1,bin_z=0.05)


   
