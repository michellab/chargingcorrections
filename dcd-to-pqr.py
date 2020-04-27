import MDAnalysis
PRM7 = "input/SYSTEM.top"
DCD = "input/traj000000001.dcd"
u = MDAnalysis.Universe(PRM7, DCD)

# select solute atoms
# write pqr file
# study how to provide customised radii. 
