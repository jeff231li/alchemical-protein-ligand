#!/bin/bash

cd electrostatics-bulk/

cp ../simulation-scripts/jq_elec_bulk jq01
sed -i "s,CHANGEWINDOW,elec-bulk,g" jq01
#qsub jq01

