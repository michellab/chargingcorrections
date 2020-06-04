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



#PRM7 = "input/SYSTEM.top"
#DCD = "input/traj000000001.dcd"
#CRD = "input/SYSTEM.crd"

PRM7 = "input/bound/SYSTEM.top"
DCD = "input/bound/traj000000001.dcd"
CRD = "input/bound/SYSTEM.crd"


# Let's load a BSS topology so we can retrieve LJ params per atoms
# and use the sigma values to guess atomic radii
rads = []
charges = []
import BioSimSpace as BSS
system = BSS.IO.readMolecules([PRM7,CRD])
molecules = system.getMolecules()
for mol in molecules:
    # skip solvent. This will include ions !!!
    if mol.getResidues()[0].name() == 'WAT':
        continue
    atoms = mol.getAtoms()
    for atom in atoms:
        sigma = atom._sire_object.property("LJ").sigma().value()
        # Heuristic is to set rad to rmin which is 2^(1/6.)*sigma/2.0
        # set minimum rad to 0.5 for polarH
        rad = max(sigma*2**(1/6.)/2., 0.5000)
        rads.append(rad)
        charge = atom._sire_object.property("charge").value()
        if mol.getResidues()[0].name() == 'LIG':
            charge = (1-L)*charge
        charges.append(charge)
#sys.exit(-1)
print (rads)
print (charges)

u = MDAnalysis.Universe(PRM7, DCD)

# select solute atoms

# This will include ions !!
selection = u.select_atoms('not resname WAT')


# construct restricted ligand Universe

system_u = MDAnalysis.core.universe.Merge(selection)
system_u.add_TopologyAttr(Radii(rads))            
system_u.add_TopologyAttr(Charges(charges))

# how to customize the charges
#topo = MDAnalysis.topology.TOPParser.TOPParser("input/SYSTEM.top")
#topology_parsed = topo.parse()
# write pqr file

system_u.atoms.write('output_'+str(L)+'.pqr')




