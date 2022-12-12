import logging
from importlib import reload

import simtk.unit as openmm_unit
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
    filename="paprika_attach.log",
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
therm_steps = 50000
equil_steps = 500000
prod_steps = 5000000
out_freq = 500

# Open system
with open("system.xml", "r") as file:
    system = XmlSerializer.deserialize(file.read())
coords = PDBFile("system.pdb")

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

# Minimize Energy
simulation.minimizeEnergy()

# Reporters
dcd_reporter = DCDReporter("thermalisation.dcd", out_freq * 10)
state_reporter = StateDataReporter(
    "thermalisation.log",
    out_freq * 10,
    step=True,
    kineticEnergy=True,
    potentialEnergy=True,
    totalEnergy=True,
    temperature=True,
    speed=True,
    separator=",",
)
simulation.reporters.append(dcd_reporter)
simulation.reporters.append(state_reporter)

# MD Step
simulation.step(therm_steps)
simulation.saveState("thermalisation.xml")

# Equilibration
# --------------------------------------------------------------------#
logging.info("Running equilibration ...")

# Thermostat
integrator = LangevinIntegrator(temperature, 1.0 / openmm_unit.picoseconds, dt_equil)

# Simulation object
simulation = Simulation(
    coords.topology,
    system,
    integrator,
    Platform.getPlatformByName("CUDA"),
    properties,
)
simulation.loadState("thermalisation.xml")

# Reporters
dcd_reporter = DCDReporter("equilibration.dcd", out_freq * 10)
state_reporter = StateDataReporter(
    "equilibration.log",
    out_freq * 10,
    step=True,
    kineticEnergy=True,
    potentialEnergy=True,
    totalEnergy=True,
    temperature=True,
    speed=True,
    separator=",",
)
simulation.reporters.append(dcd_reporter)
simulation.reporters.append(state_reporter)

# MD step
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
simulation.loadState("equilibration.xml")

# Reset current step
simulation.currentStep = 0

# Reporters
dcd_reporter = DCDReporter("production.dcd", out_freq, append=False)
check_reporter = CheckpointReporter("production.xml", out_freq * 100, writeState=True)
state_reporter = StateDataReporter(
    "production.log",
    out_freq * 10,
    step=True,
    kineticEnergy=True,
    potentialEnergy=True,
    totalEnergy=True,
    temperature=True,
    speed=True,
    separator=",",
    append=False,
)
simulation.reporters.append(dcd_reporter)
simulation.reporters.append(state_reporter)
simulation.reporters.append(check_reporter)

# Run production
simulation.step(prod_steps)
simulation.saveState("production.xml")
