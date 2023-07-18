from openmm.app import PDBFile, Modeller

protein = PDBFile("protein.noH.pdb")
ligand = PDBFile("ligand.pdb")

protein_ligand = Modeller(protein.topology, protein.positions)
protein_ligand.add(ligand.topology, ligand.positions)

with open("protein_ligand.pdb", "w") as f:
    PDBFile.writeFile(
        protein_ligand.topology,
        protein_ligand.positions,
        f,
    )

