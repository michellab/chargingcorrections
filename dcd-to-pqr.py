import MDAnalysis
from MDAnalysis.core.topologyattrs import Radii

PRM7 = "input/SYSTEM.top"
DCD = "input/traj000000001.dcd"

u = MDAnalysis.Universe(PRM7, DCD)

# select solute atoms

selection = u.select_atoms('resname LIG')


# construct restricted ligand Universe

ligand_u = MDAnalysis.core.universe.Merge(selection)


# study how to provide customised radii. 

radO = 2.27315
radC = 2.37934
radH = 1.67673
radO3 = 2.03297
radO4 = 1.98610 

rad = [radO,radC,radO,radC,radC,radC,radC,radC,radC,radO3,radC,radO4,radC,radH,radH,radH,radH,radH,radH,radH]

ligand_u.add_TopologyAttr(Radii(rad))            


# write pqr file

ligand_u.atoms.write('output/ASA.pqr')




