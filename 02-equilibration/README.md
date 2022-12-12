In this step we will run molecular dynamics simulation on the protein-ligand system to see if the ligand is stable in the binding pocket. If the ligand is not stably bound then it would not be beneficial to run the binding free energy calculations.

The snippets below is from the file `equilibration.py` and I will explain the code. In the first part of the code, we load the AMBER files we previously generated and create an OpenMM system.
```python
from openmm.app import AmberPrmtopFile, AmberInpcrdFile, HBonds, NoCutoff, OBC2
from openmm import XmlSerializer

# Load AMBER Files
prmtop = AmberPrmtopFile("../01-system_files/protein-ligand-ff19SB.prmtop")
inpcrd = AmberInpcrdFile("../01-system_files/protein-ligand-ff19SB.rst7")

# Create System - with OBC2 implicit solvent model
system = prmtop.createSystem(
    nonbondedMethod=NoCutoff,
    constraints=HBonds,
    implicitSolvent=OBC2,
    removeCMMotion=False,
)

# Save system to XML file
with open("protein_ligand.xml", "w") as f:
    f.write(XmlSerializer.serialize(system))
```
The OpenMM system will have `NoCutoff` (i.e. include all atoms in the nonbonded calculations), hydrogen bond lengths will be constrained `HBonds`, and uses the OBC2 implicit solvent model. We will also write the OpenMM system to file so than we can skip this step in the future. Once the file is written (`protein_ligand.xml`) we can check what forces have been registered by typing 
```bash
grep forceGroup protein_ligand.xml
```
you should see the following `forces`:
```bash
<Force forceGroup="0" name="HarmonicBondForce" type="HarmonicBondForce" usesPeriodic="0" version="2">
<Force forceGroup="0" name="HarmonicAngleForce" type="HarmonicAngleForce" usesPeriodic="0" version="2">
<Force forceGroup="0" name="PeriodicTorsionForce" type="PeriodicTorsionForce" usesPeriodic="0" version="2">
<Force forceGroup="0" name="CMAPTorsionForce" type="CMAPTorsionForce" usesPeriodic="0" version="2">
<Force alpha="0" cutoff="1" dispersionCorrection="1" ewaldTolerance=".0005" exceptionsUsePeriodic="0" forceGroup="0" includeDirectSpace="1" ljAlpha="0" ljnx="0" ljny="0" ljnz="0" method="0" name="NonbondedForce" nx="0" ny="0" nz="0" recipForceGroup="-1" rfDielectric="1" switchingDistance="-1" type="NonbondedForce" useSwitchingFunction="0" version="4">
<Force cutoff="1" forceGroup="0" method="0" name="GBSAOBCForce" soluteDielectric="1" solventDielectric="78.5" surfaceAreaEnergy="2.25936" type="GBSAOBCForce" version="2">
```

Next, we configure the thermostat to maintain the system temperature at a constant level and create the `simulation` object. Depending on the platform available, you will need to play around with `Platform.getPlatformByName` and the option afterwards (`"CudaDeviceIndex": "0"` tells OpenMM to run on GPU number 0 if there are multiple on your machine)
```python
import openmm.unit as openmm_unit
from openmm import LangevinIntegrator, Platform
from openmm.app import Simulation

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
    #{"CudaPrecision": "mixed", "CudaDeviceIndex": "0"},
    {"CudaPrecision": "mixed"}, 
)
simulation.context.setPositions(inpcrd.positions)
```
To help with monitoring the simulations, we will add a `StateDataReporter` that tracks the current time and energies, and a `DCDReporter` to write the trajectory.
```python
from openmm.app import DCDReporter, StateDataReporter

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
```
Finally, we will perform energy minimization followed by 10 ns of MD simulations. The state and coordinates of the final frame will also be written to file.
```python
# Energy minimization
simulation.minimizeEnergy()

# Run MD simulation (10 ns)
simulation.step(5000000)

# Save the simulation state
simulation.saveState("equilibrated.xml")

# Save final frame to file
positions = simulation.context.getState(getPositions=True).getPositions()
with open("equilibrated.pdb", "w") as f:
    PDBFile.writeFile(
        prmtop.topology,
        positions,
        f,
    )
```
View the trajectory in VMD and **make sure** that the ligand is stable in the protein. 

If stable, let's run another 3 ns of MD simulation, which we will use in the next step to determine the anchor atoms. Take a look at the file `extra_equilibration.py` to see how to perform a simulation from a previous run. Run the following to get the extra 3 ns of trajectory.
```bash
python extra_equilibration.py
```