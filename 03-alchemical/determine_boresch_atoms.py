#!/usr/bin/env python
import json
import os

import numpy as np
import openmm
import openmm.app as app
import yaml
from MDAnalysis import Universe
from MDRestraintsGenerator import search
from MDRestraintsGenerator.restraints import FindBoreschRestraint
from openff.units import unit as openff_unit
from paprika.io import PaprikaEncoder


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]


def check_colvar(atomgroup, resname="LIG"):
    count = 0
    for atom in atomgroup:
        if atom.resname == resname:
            count += 1

    if len(atomgroup) == 3:
        if count == 1:
            return "theta"
        elif count == 2:
            return "beta"

    elif len(atomgroup) == 4:
        if count == 1:
            return "phi"
        elif count == 2:
            return "alpha"
        elif count == 3:
            return "gamma"

    return ""


def extract_anchor_atoms(guest_restraints, resname="LIG"):
    if resname in guest_restraints["r"]["mask1"]:
        L1 = guest_restraints["r"]["mask1"]
        P1 = guest_restraints["r"]["mask2"]
    else:
        L1 = guest_restraints["r"]["mask2"]
        P1 = guest_restraints["r"]["mask1"]

    if (
        resname in guest_restraints["alpha"]["mask1"]
        and resname in guest_restraints["alpha"]["mask2"]
    ):
        L2 = guest_restraints["alpha"]["mask1"]
        P2 = guest_restraints["alpha"]["mask4"]
    else:
        L2 = guest_restraints["alpha"]["mask4"]
        P2 = guest_restraints["alpha"]["mask1"]

    if resname in guest_restraints["gamma"]["mask1"]:
        L3 = guest_restraints["gamma"]["mask1"]
    else:
        L3 = guest_restraints["gamma"]["mask4"]

    if resname in guest_restraints["phi"]["mask1"]:
        P3 = guest_restraints["phi"]["mask4"]
    else:
        P3 = guest_restraints["phi"]["mask1"]

    anchor_atoms = {
        "P1": guest_restraints["phi"]["mask2"],
        "P2": guest_restraints["phi"]["mask3"],
        "P3": guest_restraints["phi"]["mask4"],
        "L1": guest_restraints["gamma"]["mask3"],
        "L2": guest_restraints["gamma"]["mask2"],
        "L3": guest_restraints["gamma"]["mask1"],
    }
    return anchor_atoms


# Selection
ligand_resname = "LIG"
ligand_sel = f"resname {ligand_resname} and not name H*"
protein_sel = "protein and name CA"

# Load PDB and XML files from equilibration step
pdbfile = app.PDBFile("../02-equilibration/extra_equilibrated.pdb")
with open("../02-equilibration/protein_ligand.xml", "r") as f:
    system = openmm.XmlSerializer.deserialize(f.read())

# Create MDAnalysis Universe
universe = Universe(
    f"../01-system_files/protein-ligand-ff19SB.prmtop",
    f"../02-equilibration/extra_trajectory.dcd",
)

# Find ligand atoms
ligand_atoms = search.find_ligand_atoms(
    universe,
    l_selection=ligand_sel,
    p_align=protein_sel,
)

# Find protein atoms
atom_set = []
for l_atoms in ligand_atoms:
    psearch = search.FindHostAtoms(
        universe,
        l_atoms[0],
        p_selection=protein_sel,
    )
    psearch.run()
    atom_set.extend([(l_atoms, p) for p in psearch.host_atoms])

# Create the boresch finder analysis object
boresch = FindBoreschRestraint(universe, atom_set)

# Run the restraint analysis
boresch.run()

# Write the plots out the statistics
os.makedirs("boresch", exist_ok=True)
boresch.restraint.plot(path="./boresch")

# ---------------------------- Extract Guest Restraint --------------------------- #
guest_restraints = {
    "r": {"target": None, "mask1": None, "mask2": None},
    "theta": {"target": None, "mask1": None, "mask2": None, "mask3": None},
    "phi": {
        "target": None,
        "mask1": None,
        "mask2": None,
        "mask3": None,
        "mask4": None,
    },
    "alpha": {
        "target": None,
        "mask1": None,
        "mask2": None,
        "mask3": None,
        "mask4": None,
    },
    "beta": {"target": None, "mask1": None, "mask2": None, "mask3": None},
    "gamma": {
        "target": None,
        "mask1": None,
        "mask2": None,
        "mask3": None,
        "mask4": None,
    },
}

# ---------------------------- Extract Bond Restraint ---------------------------- #
mean_value = boresch.restraint.bond.mean
picked_value = boresch.restraint.bond.values[boresch.restraint.min_frame]
print(f"r -> (mean) {mean_value:.2f} A | (picked) {picked_value:.2f} A")

guest_restraints["r"]["target"] = picked_value * openff_unit.angstrom
for atom in boresch.restraint.bond.atomgroup:
    if atom.resname == ligand_resname:
        guest_restraints["r"]["mask1"] = f":{atom.resname}@{atom.name}"
    else:
        guest_restraints["r"]["mask2"] = f":{atom.resid}@{atom.name}"

# ---------------------------- Extract Angle Restraint --------------------------- #
for angle in boresch.restraint.angles:
    mean_value = angle.mean
    picked_value = angle.values[boresch.restraint.min_frame]

    colvar = check_colvar(angle.atomgroup)
    guest_restraints[colvar]["target"] = picked_value * openff_unit.degree

    print(f"{colvar} -> (mean) {mean_value:.2f} deg | (picked) {picked_value:.2f} deg")

    for i, atom in enumerate(angle.atomgroup):
        if atom.resname == ligand_resname:
            guest_restraints[colvar][f"mask{i+1}"] = f":{atom.resname}@{atom.name}"
        else:
            guest_restraints[colvar][f"mask{i+1}"] = f":{atom.resid}@{atom.name}"

# ---------------------------- Extract Dihedral Restraint ------------------------ #
for dihedral in boresch.restraint.dihedrals:
    mean_value = dihedral.mean
    picked_value = dihedral.values[boresch.restraint.min_frame]

    colvar = check_colvar(dihedral.atomgroup)
    guest_restraints[colvar]["target"] = picked_value * openff_unit.degree

    print(f"{colvar} -> (mean) {mean_value:.2f} deg | (picked) {picked_value:.2f} deg")

    for i, atom in enumerate(dihedral.atomgroup):
        if atom.resname == ligand_resname:
            guest_restraints[colvar][f"mask{i+1}"] = f":{atom.resname}@{atom.name}"
        else:
            guest_restraints[colvar][f"mask{i+1}"] = f":{atom.resid}@{atom.name}"

with open("boresch/guest_colvars.json", "w") as file:
    dumped = json.dumps(guest_restraints, cls=PaprikaEncoder)
    file.write(dumped)

# ------------------------------------ Anchor Atoms ------------------------------- #
anchor_atoms = extract_anchor_atoms(guest_restraints)
print(anchor_atoms)

with open("boresch/anchor_atoms.yaml", "w") as file:
    documents = yaml.dump(anchor_atoms, file)
