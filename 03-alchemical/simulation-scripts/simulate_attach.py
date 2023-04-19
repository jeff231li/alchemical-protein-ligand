import logging
from importlib import reload

import numpy as np
import openmm.unit as unit
from openmm import LangevinIntegrator, Platform, XmlSerializer
from openmm.app import (
    CheckpointReporter,
    DCDReporter,
    PDBFile,
    Simulation,
    StateDataReporter,
)


def run_molecular_dynamics(
    system_xml: str,
    pdbfile: str,
    time_steps: dict,
    phase: str,
    restart_file: str = None,
    minimize_energy: bool = False,
    write_energy: bool = True,
    write_trajectory: bool = True,
    write_checkpoint: bool = True,
    properties: dict = {"Precision": "mixed"},
):
    """
    Module to perform Molecular Dynamics simulations.
    """
    logger.info(f"Running simulation for {phase} phase ...")

    # Open system
    with open(system_xml, "r") as file:
        system = XmlSerializer.deserialize(file.read())
    coords = PDBFile(pdbfile)

    # Thermostat
    integrator = LangevinIntegrator(
        time_steps["integrator"]["temperature"],
        time_steps["integrator"]["collision"],
        time_steps["integration_time"][phase],
    )

    # Simulation Object
    simulation = Simulation(
        coords.topology,
        system,
        integrator,
        Platform.getPlatformByName("CUDA"),
        properties,
    )
    # Set positions from PDB file
    simulation.context.setPositions(coords.positions)

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
    simulation.saveState(f"{phase}.xml")
    logger.info(f"\t\tCompleted running {phase} phase ...")


# ------------------------------------------------------- #
# GPU option
# ------------------------------------------------------- #
platform_properties = {"Precision": "mixed"}
# Use the line below if you want to specify a specific GPU id on the Server
# platform_properties = {"Precision": "mixed", "DeviceIndex": "2"} #<-- change "0" to "2" for GPU no 2

# ------------------------------------------------------- #
# File Options
# ------------------------------------------------------- #
system_xml = "system.xml"
pdbfile = "system.pdb"

# Initialize logger
reload(logging)

logger = logging.getLogger()
logging.basicConfig(
    filename="paprika_attach.log",
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
        "thermalization": 1.0 * unit.femtoseconds,
        "equilibration": 2.0 * unit.femtoseconds,
        "production": 4.0 * unit.femtoseconds,
    },
    # ------------------------------------------------------- #
    # Change the time below if you need to run less MD steps
    "simulation_time": {
        "thermalization": 25.0 * unit.picoseconds,
        "equilibration": 100.0 * unit.picoseconds,
        "production": 1.0 * unit.nanoseconds,
    },
    # -------------------------------------------------------#
    "simulation_steps": {
        "thermalization": None,
        "equilibration": None,
        "production": None,
    },
    "output_frequency": {
        "thermalization": 25.0 * unit.picoseconds,
        "equilibration": 20.0 * unit.picoseconds,
        "production": 5.0 * unit.picoseconds,
    },
    "output_steps": {
        "thermalization": None,
        "equilibration": None,
        "production": None,
    },
}

# Determine number of steps
for phase in ["thermalization", "equilibration", "production"]:
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
# Minimization and Thermalization
run_molecular_dynamics(
    system_xml=system_xml,
    pdbfile=pdbfile,
    time_steps=time_steps,
    phase="thermalization",
    restart_file=None,
    minimize_energy=True,
    write_energy=True,
    write_trajectory=False,
    write_checkpoint=False,
    properties=platform_properties,
)

# Equilibration
run_molecular_dynamics(
    system_xml=system_xml,
    pdbfile=pdbfile,
    time_steps=time_steps,
    phase="equilibration",
    restart_file="thermalization.xml",
    minimize_energy=False,
    write_energy=True,
    write_trajectory=False,
    write_checkpoint=True,
    properties=platform_properties,
)

# Production
run_molecular_dynamics(
    system_xml=system_xml,
    pdbfile=pdbfile,
    time_steps=time_steps,
    phase="production",
    restart_file="equilibration.xml",
    minimize_energy=False,
    write_energy=True,
    write_trajectory=True,
    write_checkpoint=True,
    properties=platform_properties,
)
