#!/bin/bash

cd attach/

for j in $(seq 0 1 14) 
do
	if [ $j -lt 10 ] ; then jj=00${j} ; else jj=0${j} ; fi
	cd a$jj/
    # This is for running on an HPC cluster
	#cp ../../simulation-scripts/jq_attach jq01
	#sed -i "s,CHANGEWINDOW,a$jj,g" jq01
	#qsub jq01

    python ../../simulation-scripts/simulate_attach.py
	cd ../
done

