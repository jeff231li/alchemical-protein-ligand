# alchemical-protein-ligand
This repository contains python scripts to perform [absolute binding free enegy](http://www.alchemistry.org/wiki/Absolute_Binding_Free_Energy_-_Gromacs_2016) (ABFE) for a protein-ligand system. The binding free energy will be calculated using the [double-decoupling method](http://www.alchemistry.org/wiki/Decoupling_and_annihilation) (DDM) and the free energy is obtained with [thermodynamic integration](http://www.alchemistry.org/wiki/Thermodynamic_Integration) (TI).

# Dependencies
The scripts here will require some libraries, which I will list and go through the installation process. 

> **_NOTE:_** The libraries we will use will only work on Linux and MacOS. For Windows, you will need to install [WSL](https://learn.microsoft.com/en-us/windows/wsl/).

If you don't already have Anaconda installed in your machine, then go to https://anaconda.org and install version 3. Most of the dependencies will be included when we install pAPRika. I recommend installing pAPRika in a fresh `conda` environment. 

These are the steps to install pAPRika:
1. Download [the latest release](https://github.com/GilsonLabUCSD/pAPRika/releases), extract it, and change to the `paprika` directory:
2. Change the `name` field in `devtools/conda-envs/test_env.yaml` to be `paprika`.
3. Create the environment: `conda env create -f devtools/conda-envs/test_env.yaml`.
4. Activate the environment: `conda activate paprika`
5. Install `paprika` in the environment: `pip install .`

Next, we will install [MDRestraintGenerator](https://github.com/IAlibay/MDRestraintsGenerator), which we will use to determine the anchor atoms for the free energy calculations.

```bash
pip install MDRestraintsGenerator
```

Finally, we will need to install [OpenMMTools](https://github.com/choderalab/openmmtools)

```bash
conda install -c conda-forge openmmtools
```

# Tutorial
I broke this up into three sections:
1. [01-system_files](01-system_files) System preparation of protein-ligand complex
2. [02-equilibration](02-equilibration) Initial molecular dynamics (MD) simulation.
3. [03-alchemical](03-alchemical) Free energy calculations

The first section will create the topology and coordinates files using AmberTools and pAPRika modules. This assumes that we have the ligand bound to the protein (either crystal structure or docking). The second part we run molecular dynamics (MD) simulations to make sure that the ligand is stable in the binding pocket. If not, then it would not benefit performing the ABFE calculations in the third step.
