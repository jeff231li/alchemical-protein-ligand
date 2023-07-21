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
from openmmtools.alchemy import (
    AbsoluteAlchemicalFactory,
    AlchemicalRegion,
    AlchemicalState,
)


def run_molecular_dynamics_alchemical(
    system_xml: str,
    pdbfile: str,
    time_steps: dict,
    phase: str,
    ligand_resname: str,
    window_number: int = 1,
    lambda_electrostatics: float = 1.0,
    lambda_sterics: float = 1.0,
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
    logger.info(
        f"Running simulation for {phase} phase for l_elec = {lambda_electrostatics}  l_vdw = {lambda_sterics}..."
    )

    # Open system
    with open(system_xml, "r") as file:
        system = XmlSerializer.deserialize(file.read())
    coords = PDBFile(pdbfile)

    # Alchemical stuff
    alchemical_atoms = [
        atom.index
        for atom in coords.topology.atoms()
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
    alchemical_system = alchemical_factory.create_alchemical_system(
        system, alchemical_region
    )
    alchemical_state = AlchemicalState.from_system(alchemical_system)

    # Thermostat
    integrator = LangevinIntegrator(
        time_steps["integrator"]["temperature"],
        time_steps["integrator"]["collision"],
        time_steps["integration_time"][phase],
    )

    # Simulation Object
    simulation = Simulation(
        coords.topology,
        alchemical_system,
        integrator,
        Platform.getPlatformByName("CUDA"),
        properties,
    )
    # Set positions from PDB file
    simulation.context.setPositions(coords.positions)

    # Set alchemical state
    alchemical_state.lambda_electrostatics = lambda_electrostatics
    alchemical_state.lambda_sterics = lambda_sterics
    alchemical_state.apply_to_context(simulation.context)

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
            DCDReporter(
                f"{phase}-{window_number}.dcd", time_steps["output_steps"][phase]
            )
        )
    if write_checkpoint:
        simulation.reporters.append(
            CheckpointReporter(
                f"{phase}-{window_number}_restart.xml",
                time_steps["output_steps"][phase] * 10,
                writeState=True,
            )
        )
    if write_energy:
        simulation.reporters.append(
            StateDataReporter(
                f"{phase}-{window_number}.log",
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
    if phase == "production":
        simulation.saveState(f"{phase}-{window_number}.xml")
    else:
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
system_xml = "ligand.xml"
pdbfile = "ligand.pdb"

# Initialize logger
reload(logging)

logger = logging.getLogger()
logging.basicConfig(
    filename="alchemical_elec.log",
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
        "equilibration": 1.0 * unit.femtoseconds,
        "production": 2.0 * unit.femtoseconds,
    },
    # ------------------------------------------------------- #
    # Change the time below if you need to run less MD steps
    "simulation_time": {
        "equilibration": 50.0 * unit.picoseconds,
        "production": 500.0 * unit.picoseconds,
    },
    # -------------------------------------------------------#
    "simulation_steps": {
        "equilibration": None,
        "production": None,
    },
    "output_frequency": {
        "equilibration": 20.0 * unit.picoseconds,
        "production": 1.0 * unit.picoseconds,
    },
    "output_steps": {
        "equilibration": None,
        "production": None,
    },
}

# Determine number of steps
for phase in ["equilibration", "production"]:
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

time_steps["total_iterations"] = int(
    time_steps["simulation_time"]["production"]
    / time_steps["output_frequency"]["production"]
)
# ------------------------------------------------------- #
# Alchemical Options
# ------------------------------------------------------- #
# Number of Lambda windows
lambda_electrostatics = np.linspace(0, 1, 11)[::-1]
n_lambda = len(lambda_electrostatics)

lambda_sterics = 1.0

# Residue name of Ligand
resname = "LIG"

# ------------------------------------------------------- #
# Molecule Dynamics Simulation
# ------------------------------------------------------- #
# Minimization and Equilibration
run_molecular_dynamics_alchemical(
    system_xml=system_xml,
    pdbfile=pdbfile,
    time_steps=time_steps,
    phase="equilibration",
    ligand_resname=resname,
    lambda_electrostatics=lambda_electrostatics[0],
    lambda_sterics=lambda_sterics,
    restart_file=None,
    minimize_energy=True,
    write_energy=True,
    write_trajectory=False,
    write_checkpoint=True,
    properties=platform_properties,
)

# Production
for k, lb1 in enumerate(lambda_electrostatics):
    run_molecular_dynamics_alchemical(
        system_xml=system_xml,
        pdbfile=pdbfile,
        time_steps=time_steps,
        phase="production",
        ligand_resname=resname,
        window_number=k + 1,
        lambda_electrostatics=lb1,
        lambda_sterics=lambda_sterics,
        restart_file="equilibration.xml" if k == 0 else f"production-{k}.xml",
        minimize_energy=False,
        write_energy=True,
        write_trajectory=True,
        write_checkpoint=True,
        properties=platform_properties,
    )
