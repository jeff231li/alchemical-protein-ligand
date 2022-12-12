import logging
from importlib import reload

import numpy as np
import simtk.unit as openmm_unit
from openmmtools.alchemy import (
    AbsoluteAlchemicalFactory,
    AlchemicalRegion,
    AlchemicalState,
)
from simtk.openmm import LangevinIntegrator, Platform, XmlSerializer
from simtk.openmm.app import (
    CheckpointReporter,
    DCDReporter,
    PDBFile,
    Simulation,
    StateDataReporter,
)

reload(logging)

logger = logging.getLogger()
logging.basicConfig(
    filename="alchemical_elec.log",
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %I:%M:%S %p",
    level=logging.INFO,
)

# Init
properties = {"CudaPrecision": "mixed"}
temperature = 298.15 * openmm_unit.kelvin
pressure = 1.01325 * openmm_unit.bar
kT = temperature * openmm_unit.BOLTZMANN_CONSTANT_kB * openmm_unit.AVOGADRO_CONSTANT_NA
dt_therm = 1.0 * openmm_unit.femtoseconds
dt_equil = 2.0 * openmm_unit.femtoseconds
dt_prod = 2.0 * openmm_unit.femtoseconds
equil_steps = 50000
prod_steps = 500000
out_freq = 500

# MBAR stuff
resname = "LIG"
total_iterations = int(prod_steps / out_freq)
lambda_values = np.linspace(0, 1, 11)[::-1]
n_lambda = len(lambda_values)

# Open system
with open("ligand.xml", "r") as file:
    system = XmlSerializer.deserialize(file.read())
coords = PDBFile("ligand.pdb")

# Alchemical stuff
alchemical_atoms = [
    atom.index for atom in coords.topology.atoms() if atom.residue.name == resname
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
system = alchemical_factory.create_alchemical_system(system, alchemical_region)
alchemical_state = AlchemicalState.from_system(system)

# Minimization and Thermalisation
# --------------------------------------------------------------------#
logger.info("Minimizing and Thermalisation ...")

# Thermostat
integrator = LangevinIntegrator(temperature, 1.0 / openmm_unit.picoseconds, dt_therm)

# Simulation Object
simulation = Simulation(
    coords.topology,
    system,
    integrator,
    Platform.getPlatformByName("CUDA"),
    properties,
)
simulation.context.setPositions(coords.positions)

# Set alchemical state
alchemical_state.lambda_sterics = 1.0
alchemical_state.lambda_electrostatics = lambda_values[0]
alchemical_state.apply_to_context(simulation.context)

# Minimize Energy
simulation.minimizeEnergy()

# Reporters
check_reporter = CheckpointReporter("equilibration.chk", out_freq * 10)
dcd_reporter = DCDReporter("equilibration.dcd", out_freq * 10)
state_reporter = StateDataReporter(
    "equilibration.log",
    out_freq * 10,
    step=True,
    kineticEnergy=True,
    potentialEnergy=True,
    totalEnergy=True,
    temperature=True,
    volume=True,
    speed=True,
    separator=",",
)
simulation.reporters.append(check_reporter)
simulation.reporters.append(dcd_reporter)
simulation.reporters.append(state_reporter)

# MD Step
simulation.step(equil_steps)
simulation.saveState("equilibration.xml")

# Production
# --------------------------------------------------------------------#
logging.info("Running production ...")

# Thermostat
integrator = LangevinIntegrator(temperature, 1.0 / openmm_unit.picoseconds, dt_prod)

# Create simulation object
simulation = Simulation(
    coords.topology,
    system,
    integrator,
    Platform.getPlatformByName("CUDA"),
    properties,
)

# Run over lambda
for k, lb1 in enumerate(lambda_values):
    logging.info(f"Running short equilibration for lambda: {lb1}")
    alchemical_state.lambda_sterics = 1.0
    alchemical_state.lambda_electrostatics = lb1
    alchemical_state.apply_to_context(simulation.context)

    # Names
    if k == 0:
        previous = "equilibration"
    else:
        previous = f"production-{k}"
    current = f"production-{k+1}"

    simulation.loadState(f"{previous}.xml")
    simulation.currentStep = 0

    # Short equilibration step before production
    simulation.step(equil_steps)

    # Reporters
    check_reporter = CheckpointReporter(f"{current}.xml", out_freq * 10, writeState=True)
    dcd_reporter = DCDReporter(f"{current}.dcd", out_freq)
    state_reporter = StateDataReporter(
        f"{current}.log",
        out_freq * 10,
        step=True,
        kineticEnergy=True,
        potentialEnergy=True,
        totalEnergy=True,
        temperature=True,
        volume=True,
        speed=True,
        separator=",",
    )
    simulation.reporters.append(check_reporter)
    simulation.reporters.append(dcd_reporter)
    simulation.reporters.append(state_reporter)

    # MD step
    logging.info(f"Running production for lambda: {lb1}")
    for j in range(total_iterations):
        alchemical_state.lambda_sterics = 1.0
        alchemical_state.lambda_electrostatics = lb1
        alchemical_state.apply_to_context(simulation.context)
        simulation.step(out_freq)

    simulation.saveState(f"{current}.xml")

    simulation.reporters.pop(-1)
    simulation.reporters.pop(-1)
    simulation.reporters.pop(-1)
