import MDAnalysis
from MDAnalysis.core.topologyattrs import Radii
from MDAnalysis.core.topologyattrs import Charges
import MDAnalysis.topology.TOPParser
import os

#assigning lambda value

print('enter lambda value :')

L = float(input())


# towards the automation of the code 

#path = __file__
#dirpath = os.path.dirname(path)
#dirname = os.path.basename(dirpath)

#L = float(dirname[7:])



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

# how to customize the charges

topo = MDAnalysis.topology.TOPParser.TOPParser("input/SYSTEM.top")

topology_parsed = topo.parse()

tab_charges = topology_parsed.charges.values

first_charges = []

for i in range(len(selection)):
    first_charges.append(tab_charges[i])
  
lambda_charges = []

for charge in first_charges:
    lambda_charges.append(charge*(1-L))

ligand_u.add_TopologyAttr(Charges(lambda_charges))
    
# write pqr file

ligand_u.atoms.write('pqr_snapshots/LIG_free_1st_L_'+str(L)+'.pqr', frames=[0])




