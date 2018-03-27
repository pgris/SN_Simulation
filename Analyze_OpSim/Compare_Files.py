from Observations import *

def Cadence(obs):

    oob=obs['mjd'][1:]-obs['mjd'][:-1]

    return np.median(oob)
def sel(tab,band):

    idx = tab['band']=='LSSTPG::'+band
    return tab[idx]

filea='../../Ana_Cadence/OpSimLogs_test/DD/Observations_DD_290.txt'
fileb='OpSimLogs/DD/Observations_DD_290.txt'

obsa=Observations(fieldid=290, filename=filea)
obsb=Observations(fieldid=290, filename=fileb)

bands='g'
print len(obsa.seasons),len(obsb.seasons)
for season in range(len(obsa.seasons)):
    seasa=obsa.seasons[season]
    seasb=obsb.seasons[season]

    #print seasa
    for band in bands:
        sela=sel(seasa,band)
        selb=sel(seasb,band)
        print Cadence(sela),Cadence(selb),len(sela),len(selb)
        
        for i in range(len(sela)):
            print sela[i]['mjd'],selb[i]['mjd'],sela[i]['m5sigmadepth']-selb[i]['m5sigmadepth']

    break
