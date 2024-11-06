import logging
from importlib import reload

import numpy as np
import openmm.unit as unit
from openmm import LangevinIntegrator, Platform, System, XmlSerializer
from openmm.app import (
    OBC2,
    AmberInpcrdFile,
    AmberPrmtopFile,
    CheckpointReporter,
    DCDReporter,
    HBonds,
    NoCutoff,
    PDBFile,
    Simulation,
    StateDataReporter,
)


def run_molecular_dynamics(
    system: System,
    prmtop: AmberPrmtopFile,
    inpcrd: AmberInpcrdFile,
    time_steps: dict,
    phase: str,
    restart_file: str = None,
    minimize_energy: bool = False,
    write_energy: bool = True,
    write_trajectory: bool = True,
    write_checkpoint: bool = True,
    write_pdb: bool = True,
    properties: dict = {"Precision": "mixed"},
):
    """
    Module to perform Molecular Dynamics simulations.
    """
    logger.info(f"Running simulation for {phase} phase ...")

    # Thermostat
    integrator = LangevinIntegrator(
        time_steps["integrator"]["temperature"],
        time_steps["integrator"]["collision"],
        time_steps["integration_time"][phase],
    )

    # Simulation Object
    simulation = Simulation(
        prmtop.topology,
        system,
        integrator,
        Platform.getPlatformByName("CUDA"),
        properties,
    )
    # Set positions from PDB file
    simulation.context.setPositions(inpcrd.positions)

    # Restart from previous run
    if restart_file:
        logger.info(f"\t\tRestarting from {restart_file} ...")
        simulation.loadState(restart_file)
        simulation.currentStep = 0

    # Minimize Energy
    if minimize_energy:
        logger.info("\t\tMinimizing energy ...")
        simulation.minimizeEnergy()

    # Reporters
    if write_trajectory:
        simulation.reporters.append(
            DCDReporter(f"{phase}.dcd", time_steps["output_steps"][phase])
        )
    if write_checkpoint:
        simulation.reporters.append(
            CheckpointReporter(
                f"{phase}_restart.xml",
                time_steps["output_steps"][phase] * 10,
                writeState=True,
            )
        )
    if write_energy:
        simulation.reporters.append(
            StateDataReporter(
                f"{phase}.log",
                time_steps["output_steps"][phase],
                step=True,
                time=True,
                potentialEnergy=True,
                kineticEnergy=True,
                totalEnergy=True,
                temperature=True,
                speed=True,
                separator=",",
            )
        )

    # MD Step
    simulation.step(time_steps["simulation_steps"][phase])
    logger.info(f"\t\tCompleted running {phase} phase ...")

    # Save files
    simulation.saveState(f"{phase}.xml")

    # Save final frame to file
    if write_pdb:
        positions = simulation.context.getState(getPositions=True).getPositions()
        with open(f"{phase}.pdb", "w") as f:
            PDBFile.writeFile(
                prmtop.topology,
                positions,
                f,
            )


# ------------------------------------------------------- #
# GPU option
# ------------------------------------------------------- #
platform_properties = {"Precision": "mixed"}
# Use the line below if you want to specify a specific GPU id on the Server
# platform_properties = {"Precision": "mixed", "DeviceIndex": "2"} #<-- change "0" to "2" for GPU no 2

# ------------------------------------------------------- #
# File Options
# ------------------------------------------------------- #
# Load AMBER Files
prmtop = AmberPrmtopFile("../01-system_files/protein-ligand-ff19SB.prmtop")
inpcrd = AmberInpcrdFile("../01-system_files/protein-ligand-ff19SB.rst7")

# Create System - in vacuum with OBC2 implicit solvent model
system = prmtop.createSystem(
    nonbondedMethod=NoCutoff,
    constraints=HBonds,
    implicitSolvent=OBC2,
    removeCMMotion=False,
)

# Save system to XML file
with open("protein_ligand.xml", "w") as f:
    f.write(XmlSerializer.serialize(system))

# Initialize logger
reload(logging)

logger = logging.getLogger()
logging.basicConfig(
    filename="paprika_equilibration.log",
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %I:%M:%S %p",
    level=logging.INFO,
)

# ------------------------------------------------------- #
# Setup Simulation Time and Output
# ------------------------------------------------------- #
time_steps = {
    "integrator": {
        "temperature": 298.15 * unit.kelvin,
        "collision": 1.0 / unit.picoseconds,
    },
    "integration_time": {
        "production-extra": 2.0 * unit.femtoseconds,
    },
    # ------------------------------------------------------- #
    # Change the time below if you need to run less MD steps
    "simulation_time": {
        "production-extra": 3.0 * unit.nanoseconds,
    },
    # -------------------------------------------------------#
    "simulation_steps": {
        "production-extra": None,
    },
    "output_frequency": {
        "production-extra": 5.0 * unit.picoseconds,
    },
    "output_steps": {
        "production-extra": None,
    },
}

# Determine number of steps
for phase in ["production-extra"]:
    time_steps["simulation_steps"][phase] = int(
        np.ceil(
            time_steps["simulation_time"][phase] / time_steps["integration_time"][phase]
        )
    )
    time_steps["output_steps"][phase] = int(
        np.ceil(
            time_steps["output_frequency"][phase]
            / time_steps["integration_time"][phase]
        )
    )
    logger.info(f"{phase} intergation time -- {time_steps['integration_time'][phase]}")
    logger.info(f"{phase} simulation time  -- {time_steps['simulation_time'][phase]}")
    logger.info(f"{phase} output steps     -- {time_steps['output_steps'][phase]}")
    logger.info(f"{phase} total steps      -- {time_steps['simulation_steps'][phase]}")

# ------------------------------------------------------- #
# Molecule Dynamics Simulation
# ------------------------------------------------------- #
# Production
run_molecular_dynamics(
    system=system,
    prmtop=prmtop,
    inpcrd=inpcrd,
    time_steps=time_steps,
    phase="production-extra",
    restart_file="production.xml",
    minimize_energy=False,
    write_energy=True,
    write_trajectory=True,
    write_checkpoint=True,
    write_pdb=True,
    properties=platform_properties,
)
