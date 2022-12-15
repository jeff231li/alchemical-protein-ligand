Here, we will generate the files for running the double-decoupling method (DDM). We will calculate the contributions from the electrostatics and lennard-jones potential separately. For the electrostatics we will annihilate the charges but for the lennard-jones we will use decoupling. The protocol follows that from this [paper](https://www.nature.com/articles/s42004-022-00721-4). The calculation for the absolute binding free energy (ABFE) is broken up into five stages:
1. (Binding site) apply restraints on the ligand relative to the protein
2. (Binding site) annihilate partial charges on the ligand
3. (Binding site) decouple lennard-jones interactions between protein and ligand
4. (Bulk) scale partial charges on the ligand
5. (Bulk) release restraints on ligand and add correction for 1 M standard-state.

The last part is estimated analytically so no simulation will be needed. I've configured the scripts to make signs  match so that the binding free energy is just the sum of all the components: 


DG_b = DG(attach) + DG(elec_site) + DG(vdw) + DG(elec_bulk) + DG(release)
