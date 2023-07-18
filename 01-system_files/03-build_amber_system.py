from paprika.build.system import TLeap

# Build AMBER Topologies for Ligand
system = TLeap()
system.output_path = "./"
system.output_prefix = "ligand-ff14SB"
system.neutralize = False
system.pbc_type = None
system.template_lines = [
    "set default PBRadii mbondi2",
    "source leaprc.gaff2",
    "loadamberparams ligand.frcmod",
    "LIG = loadmol2 ligand.am1bcc.gaff2.mol2",
    "model = loadpdb ligand.pdb",
]
system.build(clean_files=False)
# Do HMR so we can increase integration time step
system.repartition_hydrogen_mass()

# Build AMBER Topologies for Protein-Ligand system
system = TLeap()
system.output_path = "./"
system.output_prefix = "protein-ligand-ff14SB"
system.neutralize = False
system.pbc_type = None
system.template_lines = [
    "set default PBRadii mbondi2",
    "source leaprc.protein.ff14SB",
    #"source leaprc.water.opc",
    "source leaprc.gaff2",
    "loadamberparams ligand.frcmod",
    "LIG = loadmol2 ligand.am1bcc.gaff2.mol2",
    "model = loadpdb protein_ligand.pdb",
]
system.build(clean_files=False)
# Do HMR so we can increase integration time step
system.repartition_hydrogen_mass()

