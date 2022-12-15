#!/bin/bash

ligand=$1

cd ${ligand}/

cd electrostatics-bulk/

cp ../../jq_elec_bulk jq01
sed -i "s,CHANGEWINDOW,${ligand},g" jq01
qsub jq01

cd ../

cd ../
