import openmm.unit as openmm_unit
from openmm import LangevinIntegrator, Platform, XmlSerializer
from openmm.app import (
    OBC2,
    AmberInpcrdFile,
    AmberPrmtopFile,
    DCDReporter,
    HBonds,
    NoCutoff,
    Simulation,
    StateDataReporter,
)

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

temperature = 298.15 * openmm_unit.kelvin
collision = 1.0 / openmm_unit.picoseconds
time_step = 2.0 * openmm_unit.femtoseconds

# Langevin Thermostat
integrator = LangevinIntegrator(temperature, collision, time_step)

# Simulation object
simulation = Simulation(
    prmtop.topology,
    system,
    integrator,
    Platform.getPlatformByName("CUDA"),
    # {"CudaPrecision": "mixed", "CudaDeviceIndex": "0"},
    {"CudaPrecision": "mixed"},
)
simulation.context.setPositions(inpcrd.positions)

# Add Reporters
state_reporter = StateDataReporter(
    "openmm_statistics.csv",
    5000,
    step=True,
    kineticEnergy=True,
    potentialEnergy=True,
    totalEnergy=True,
    temperature=True,
    speed=True,
    separator=",",
)
dcd_reporter = DCDReporter(
    "trajectory.dcd",
    5000,
)
simulation.reporters.append(state_reporter)
simulation.reporters.append(dcd_reporter)

# Energy minimization
simulation.minimizeEnergy()

# Run MD simulation (10 ns)
simulation.step(10000)
