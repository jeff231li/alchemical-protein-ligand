#!/bin/bash

cd electrostatics-bulk/

# This is for running on an HPC cluster
#cp ../simulation-scripts/jq_elec_bulk jq01
#sed -i "s,CHANGEWINDOW,elec-bulk,g" jq01
#qsub jq01

python ../simulation-scripts/simulate_electrostatics-bulk.py

