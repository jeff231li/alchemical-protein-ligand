import openmm.unit as openmm_unit
from openmm import LangevinIntegrator, Platform, XmlSerializer
from openmm.app import DCDReporter, PDBFile, Simulation, StateDataReporter

# Load OpenMM Files
pdbfile = PDBFile("equilibrated.pdb")
with open("protein_ligand.xml", "r") as f:
    system = XmlSerializer.deserialize(f.read())

temperature = 298.15 * openmm_unit.kelvin
collision = 1.0 / openmm_unit.picoseconds
time_step = 2.0 * openmm_unit.femtoseconds

# Langevin Thermostat
integrator = LangevinIntegrator(temperature, collision, time_step)

# Simulation object
simulation = Simulation(
    pdbfile.topology,
    system,
    integrator,
    Platform.getPlatformByName("CUDA"),
    # {"CudaPrecision": "mixed", "CudaDeviceIndex": "0"},
    {"CudaPrecision": "mixed"},
)
simulation.context.setPositions(pdbfile.positions)

# Add Reporters
state_reporter = StateDataReporter(
    "extra_statistics.csv",
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
    "extra_trajectory.dcd",
    5000,
)
simulation.reporters.append(state_reporter)
simulation.reporters.append(dcd_reporter)

# Run MD simulation (3 ns)
simulation.step(1500000)

# Write final frame to file
with open("extra_equilibrated.pdb", "w") as f:
    PDBFile.writeFile()
