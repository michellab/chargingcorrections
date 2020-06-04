import MDAnalysis
from MDAnalysis.core.topologyattrs import Radii
from MDAnalysis.core.topologyattrs import Charges
import MDAnalysis.topology.TOPParser
import os

lambda_values = [0.000, 0.050, 0.100, 0.200, 0.300, 0.400, 0.500, 0.600, 0.700, 0.800, 0.900, 1.000]

PRM7 = "../input/SYSTEM.top"
CRD="../input/SYSTEM.crd"

rads=[]

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

for L in lambda_values : # lambda_values is a list of lambda values used during the discharge step that will have to be defined

                
    # go to lambda folder
        
    os.chdir('lambda-'+format(L, '.3f'))

    DCD = "traj000000001.dcd"
    PRM7 = "../../input/SYSTEM.top"
    
    # Let's load a BSS topology so we can retrieve LJ params per atoms
    # and use the sigma values to guess atomic radii
    charges = []
    
    for mol in molecules:
        # skip solvent. This will include ions !!!
        if mol.getResidues()[0].name() == 'WAT':
            continue
        atoms = mol.getAtoms()
        for atom in atoms:
            charge = atom._sire_object.property("charge").value()
            if mol.getResidues()[0].name() == 'LIG':
                charge = (1-L)*charge
            charges.append(charge)
    #sys.exit(-1)
    
    #creating the universe
      
    u = MDAnalysis.Universe(PRM7, DCD)
    

    # select solute atoms

    #selection = u.select_atoms('resname LIG')
    
    #to include ions we should use 
    selection = u.select_atoms('not resname WAT')


    # construct restricted ligand or complex Universe

    
    system_u = MDAnalysis.core.universe.Merge(selection)
    system_u.add_TopologyAttr(Radii(rads))            
    system_u.add_TopologyAttr(Charges(charges))


    # write pqr file

    system_u.atoms.write('../../pqr_snapshots/rad_bss/LIG_opt_free_L_'+str(L)+'.pqr')

    os.chdir("../")