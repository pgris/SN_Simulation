import glob
import os

files = glob.glob('logs/*')
what='ImportError'

for myfile in files:
    if what in open(myfile).read():
        print('problem here',myfile)
        torepro=myfile.replace('logs','scripts').replace('.log','.sh')
        print torepro
        #os.system('sh '+torepro)
