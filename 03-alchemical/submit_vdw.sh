#!/bin/bash

cd lennard-jones/

# 11 Lambda windows
categories=(0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)

# 21 Lambda windows
#categories=(0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0)

j=0
for lb in "${categories[@]}"
do
	if [ $j -lt 10 ] ; then jj=00${j} ; else jj=0${j} ; fi
	cd v$jj/
    # This is for running on an HPC cluster
    #cp ../../simulation-scripts/jq_vdw jq01
	#sed -i "s,CHANGEWINDOW,v$jj,g" jq01
	#sed -i "s,CHANGELAMBDA,$lb,g" jq01
	#qsub jq01

    python ../../simulation-scripts/simulate_lennard-jones.py $lb
	j=$(expr $j + 1)
	cd ../
done
