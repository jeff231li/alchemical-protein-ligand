#!/bin/bash

cd electrostatics-site/

# 11 Lambda windows
categories=(0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)

j=0
for lb in "${categories[@]}"
do
	if [ $j -lt 10 ] ; then jj=00${j} ; else jj=0${j} ; fi

	cd e$jj/
    # This is for running on an HPC cluster
    #cp ../../simulation-scripts/jq_elec_site jq01
	#sed -i "s,CHANGEWINDOW,e$jj,g" jq01
	#sed -i "s,CHANGELAMBDA,$lb,g" jq01
	#qsub jq01

    python ../../simulation-scripts/simulate_electrostatics-site.py $lb
	cd ../
    j=$(expr $j + 1)
done
