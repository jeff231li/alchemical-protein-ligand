#!/bin/bash

ligand=$1

cd ${ligand}/

cd electrostatics-site/

categories=(0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)
#categories=(0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0)

j=0
for lb in "${categories[@]}"
do
	if [ $j -lt 10 ] ; then jj=00${j} ; else jj=0${j} ; fi

	#if [ $(expr $j % 2) -ne 0 ] ; then
	    cd e$jj/
	    echo "$(pwd) $lb"
	    cp ../../../jq_elec_site jq01
	    sed -i "s,CHANGEWINDOW,e$jj,g" jq01
	    sed -i "s,CHANGELAMBDA,$lb,g" jq01
	    qsub jq01
	    cd ../
	#fi
    j=$(expr $j + 1)
done

cd ../

cd ../
