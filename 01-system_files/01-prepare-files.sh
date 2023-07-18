#!/bin/bash

# Prepare Ligand - Generate MOL2 with AM1-BCC charges
antechamber -fi pdb -fo mol2 -i ligand.pdb -o ligand.am1bcc.gaff2.mol2 -at gaff2 -c bcc -rn LIG -pf y

# Prepare Ligand - Generate frcmod file
parmchk2 -f mol2 -i ligand.am1bcc.gaff2.mol2 -s gaff2 -o ligand.frcmod

# Prepare protein - Strip all Hydrogen atoms
pdb4amber --nohyd protein.pdb > protein.noH.pdb

exit
