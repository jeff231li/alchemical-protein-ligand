import json
import os
import shutil

import numpy as np
import openmm.unit as openmm_unit
import openmm as openmm
import openmm.app as app
import parmed as pmd
import yaml
from openff.units import unit as openff_unit
from openmmtools.alchemy import (
    AbsoluteAlchemicalFactory,
    AlchemicalRegion,
    AlchemicalState,
)
from paprika.io import PaprikaDecoder, save_restraints
from paprika.restraints import DAT_restraint, create_window_list
from paprika.restraints.openmm import apply_dat_restraint
from paprika.restraints.utils import parse_window
from tqdm import tqdm


def get_guest_restraints(
    structure,
    anchor_atoms,
    guest_colvars,
    lambda_values,
    kdist,
    kangle,
    return_wall=False,
):
    guest_restraints = []

    # Translational Restraints
    # colvar - r
    r = DAT_restraint()
    r.mask1 = anchor_atoms["P1"]
    r.mask2 = anchor_atoms["L1"]
    r.topology = structure
    r.auto_apr = False
    r.continuous_apr = False
    r.amber_index = False

    r.attach["target"] = guest_colvars["r"]["target"]
    r.attach["fraction_list"] = lambda_values
    r.attach["fc_final"] = kdist

    r.initialize()
    guest_restraints.append(r)

    # colvar - theta
    r = DAT_restraint()
    r.mask1 = anchor_atoms["P2"]
    r.mask2 = anchor_atoms["P1"]
    r.mask3 = anchor_atoms["L1"]
    r.topology = structure
    r.auto_apr = False
    r.continuous_apr = False
    r.amber_index = False

    r.attach["target"] = guest_colvars["theta"]["target"]
    r.attach["fraction_list"] = lambda_values
    r.attach["fc_final"] = kangle

    r.initialize()
    guest_restraints.append(r)

    # colvar - phi
    r = DAT_restraint()
    r.mask1 = anchor_atoms["P3"]
    r.mask2 = anchor_atoms["P2"]
    r.mask3 = anchor_atoms["P1"]
    r.mask4 = anchor_atoms["L1"]
    r.topology = structure
    r.auto_apr = False
    r.continuous_apr = False
    r.amber_index = False

    r.attach["target"] = guest_colvars["phi"]["target"]
    r.attach["fraction_list"] = lambda_values
    r.attach["fc_final"] = kangle

    r.initialize()
    guest_restraints.append(r)

    # Rotational restraints
    # colvar - alpha
    r = DAT_restraint()
    r.mask1 = anchor_atoms["P2"]
    r.mask2 = anchor_atoms["P1"]
    r.mask3 = anchor_atoms["L1"]
    r.mask4 = anchor_atoms["L2"]
    r.topology = structure
    r.auto_apr = False
    r.continuous_apr = False
    r.amber_index = False

    r.attach["target"] = guest_colvars["alpha"]["target"]
    r.attach["fraction_list"] = lambda_values
    r.attach["fc_final"] = kangle

    r.initialize()
    guest_restraints.append(r)

    # colvar - beta
    r = DAT_restraint()
    r.mask1 = anchor_atoms["P1"]
    r.mask2 = anchor_atoms["L1"]
    r.mask3 = anchor_atoms["L2"]
    r.topology = structure
    r.auto_apr = False
    r.continuous_apr = False
    r.amber_index = False

    r.attach["target"] = guest_colvars["beta"]["target"]
    r.attach["fraction_list"] = lambda_values
    r.attach["fc_final"] = kangle

    r.initialize()
    guest_restraints.append(r)

    # colvar - gamma
    r = DAT_restraint()
    r.mask1 = anchor_atoms["P1"]
    r.mask2 = anchor_atoms["L1"]
    r.mask3 = anchor_atoms["L2"]
    r.mask4 = anchor_atoms["L3"]
    r.topology = structure
    r.auto_apr = False
    r.continuous_apr = False
    r.amber_index = False

    r.attach["target"] = guest_colvars["gamma"]["target"]
    r.attach["fraction_list"] = lambda_values
    r.attach["fc_final"] = kangle

    r.initialize()
    guest_restraints.append(r)

    # Wall Restraint
    wall_restraints = []

    r = DAT_restraint()
    r.mask1 = anchor_atoms["P1"]
    r.mask2 = anchor_atoms["L1"]
    r.topology = structure
    r.auto_apr = False
    r.continuous_apr = False
    r.amber_index = False

    r.attach["target"] = guest_colvars["r"]["target"] + 8.0 * openff_unit.angstrom
    r.attach["fraction_list"] = [1.0] * len(lambda_values)
    r.attach["fc_final"] = kdist

    r.custom_restraint_values["r1"] = 0.0
    r.custom_restraint_values["r2"] = 0.0
    r.custom_restraint_values["rk2"] = 0.0

    r.initialize()
    wall_restraints.append(r)

    if return_wall is True:
        return guest_restraints, wall_restraints

    return guest_restraints


# Restraints and fractions lists
kdist = 5.0 * openff_unit.kcal / openff_unit.mole / openff_unit.angstrom**2
kangle = 100.0 * openff_unit.kcal / openff_unit.mole / openff_unit.radian**2

# 15 non-linearly spaced windows
attach_string = (
    "0.00 0.40 0.80 1.60 2.40 4.00 5.50 8.65 11.80 18.10 24.40 37.00 49.60 74.80 100.00"
)
attach_fractions = [float(i) / 100 for i in attach_string.split()]

# 11 electrostatics and 21 lennard-jones
lambda_elec_values = np.linspace(0, 1, 11)[::-1]
lambda_vdw_values = np.linspace(0, 1, 21)[::-1]

# Load guest colvars
with open("boresch/guest_colvars.json", "r") as file:
    guest_colvars = json.load(file, cls=PaprikaDecoder)

# Load Anchor Atoms
with open("boresch/anchor_atoms.yaml", "r") as file:
    anchor_atoms = yaml.safe_load(file)

pdbfile = app.PDBFile("../02-equilibration/extra_equilibrated.pdb")
with open("../02-equilibration/protein_ligand.xml", "r") as f:
    system = openmm.XmlSerializer.deserialize(f.read())

structure = pmd.openmm.load_topology(pdbfile.topology, system, xyz=pdbfile.positions)

# Alchemical stuff
ligand_resname = "LIG"
alchemical_atoms = [
    atom.index
    for atom in pdbfile.topology.atoms()
    if atom.residue.name == ligand_resname
]
alchemical_factory = AbsoluteAlchemicalFactory(
    consistent_exceptions=False,
    alchemical_pme_treatment="exact",
    disable_alchemical_dispersion_correction=False,
)
alchemical_region = AlchemicalRegion(
    alchemical_atoms=alchemical_atoms,
    annihilate_electrostatics=True,
    annihilate_sterics=False,
)

# ---------------------------- Create Attach Windows --------------------------------- #
# Generate Guest restraints
guest_restraints, wall_restraints = get_guest_restraints(
    structure,
    anchor_atoms,
    guest_colvars,
    attach_fractions,
    kdist,
    kangle,
    return_wall=True,
)
# Save restraints to file
save_restraints(
    guest_restraints + wall_restraints, filepath="boresch/boresch_restraints.json"
)

os.makedirs("attach", exist_ok=True)
window_list = create_window_list(guest_restraints)

# Write Attach Windows
for window in tqdm(window_list):
    folder = f"attach/{window}"
    os.makedirs(folder, exist_ok=True)

    window_number, phase = parse_window(window)

    #print(f"Creating XML for window {window} ...")

    # Copy PDBFile
    shutil.copy("../02-equilibration/extra_equilibrated.pdb", f"{folder}/system.pdb")

    # Load OpenMM XML
    with open("../02-equilibration/protein_ligand.xml", "r") as file:
        system = openmm.XmlSerializer.deserialize(file.read())

    # Add restraints to system
    for restraint in guest_restraints:
        apply_dat_restraint(system, restraint, phase, window_number, force_group=10)
    for restraint in wall_restraints:
        apply_dat_restraint(system, restraint, phase, window_number, force_group=11)

    # Update XML object
    system_xml = openmm.XmlSerializer.serialize(system)
    with open(f"{folder}/system.xml", "w") as file:
        file.write(system_xml)

# ---------------------------- Create Elec-Bulk Files --------------------------------- #
os.makedirs("electrostatics-bulk", exist_ok=True)

prmtop = app.AmberPrmtopFile("../01-system_files/ligand-ff19SB.prmtop")
inpcrd = app.AmberInpcrdFile("../01-system_files/ligand-ff19SB.rst7")
system = prmtop.createSystem(
    nonbondedMethod=app.NoCutoff,
    constraints=app.HBonds,
    implicitSolvent=app.OBC2,
    removeCMMotion=False,
)
with open("electrostatics-bulk/ligand.pdb", "w") as f:
    app.PDBFile.writeFile(
        prmtop.topology,
        inpcrd.positions,
        f,
    )
with open("electrostatics-bulk/ligand.xml", "w") as f:
    f.write(openmm.XmlSerializer.serialize(system))

# ---------------------------- Create Elec-Site Windows --------------------------------- #
os.makedirs("electrostatics-site", exist_ok=True)

with open("../02-equilibration/protein_ligand.xml", "r") as file:
    system = openmm.XmlSerializer.deserialize(file.read())

# Add restraints
for restraint in guest_restraints:
    apply_dat_restraint(
        system, restraint, "attach", len(attach_fractions) - 1, force_group=10
    )

# Convert to Alchemical System
system_alch = alchemical_factory.create_alchemical_system(system, alchemical_region)
integrator = openmm.VerletIntegrator(0.001*openmm_unit.picoseconds)
context = openmm.Context(system_alch, integrator)

# Write individual windows
for i, lb in enumerate(tqdm(lambda_elec_values)):
    window = f"e{i:03}"
    os.makedirs(f"electrostatics-site/{window}", exist_ok=True)
    shutil.copy(
        "../02-equilibration/extra_equilibrated.pdb",
        f"electrostatics-site/{window}/system.pdb",
    )
    with open(f"electrostatics-site/{window}/system.xml", "w") as file:
        file.write(openmm.XmlSerializer.serialize(system_alch))

# ---------------------------- Create vdw Windows --------------------------------- #
os.makedirs("lennard-jones", exist_ok=True)

with open("../02-equilibration/protein_ligand.xml", "r") as file:
    system = openmm.XmlSerializer.deserialize(file.read())

# Add restraints
for restraint in guest_restraints:
    apply_dat_restraint(
        system, restraint, "attach", len(attach_fractions) - 1, force_group=10
    )

# Convert to Alchemical System
system_alch = alchemical_factory.create_alchemical_system(system, alchemical_region)

for i, lb in enumerate(tqdm(lambda_vdw_values)):
    window = f"v{i:03}"
    os.makedirs(f"lennard-jones/{window}", exist_ok=True)
    shutil.copy(
        "../02-equilibration/extra_equilibrated.pdb",
        f"lennard-jones/{window}/system.pdb",
    )
    with open(f"lennard-jones/{window}/system.xml", "w") as file:
        file.write(openmm.XmlSerializer.serialize(system_alch))
