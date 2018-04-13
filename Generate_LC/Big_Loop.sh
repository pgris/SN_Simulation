#!/bin/bash

fieldname=$1
fieldid=$2
season=$3
stretch=$4
color=$5
T0step=$6
zstep=0.025


#python Loop_Prod.py --fieldname ${fieldname} --fieldid ${fieldid} --season ${season} --stretch ${stretch} --color ${color} --dirout Light_Curves_sncosmo --NT0 -1 --T0step ${T0step} --zmin 0.0 --zmax 0.1 --zstep ${zstep} --Nz 4


python Loop_Prod.py --fieldname ${fieldname} --fieldid ${fieldid} --season ${season} --stretch ${stretch} --color ${color} --dirout Light_Curves_sncosmo --NT0 -1 --T0step ${T0step} --zmin 0.8 --zmax 0.825 --zstep ${zstep} --Nz 1



