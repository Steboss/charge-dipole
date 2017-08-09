from parmed.amber import * 
import sys


top = sys.argv[1]
crd = sys.argv[2]
base = AmberParm(top,crd)

base.residues[0].atoms[0].charge = -0.834
base.residues[0].atoms[1].charge =  0.834

base.residues[0].atoms[0].bonds[0].type.req = 0.586

base.remake_parm()
base.write_parm("dipole.prmtop")
base.write_rst7("dipole.rst7")
print("Warning: the final coordinates file must be fixed")
