#May 2017 Stefano Bosisio
#Script to merge sodium topology with water one. in this way we will have a
#sodium + 1 tip3p water molecule topology and crd
#Usage: python merge.py  na.prmtop na.rst7 water.prmtop water.rst7
import parmed
from parmed.amber import *
import os,sys,math,time,subprocess,random
import multiprocessing
from multiprocessing import Pool, Process
from functools import partial


###MAIN SCRIPT###

natop = sys.argv[1]
nacrd = sys.argv[2]

wattop = sys.argv[3]
watcrd = sys.argv[4]
print("Loading topologies...")
base = AmberParm(wattop,watcrd)
ion = AmberParm(natop,nacrd)
res_numb = len(base.residues) + 1
print("Creation of sodium ion...")
ion_atom = ion.residues[0].atoms[0]
ion_atom.xx = random.uniform(1,3)
ion_atom.xy = random.uniform(1,3)
ion_atom.xz = random.uniform(1,3)
print("Adding sodium to water topology...")
base.add_atom(ion_atom,ion_atom.name,res_numb)
#now  recompute vdW
base.remake_parm()
print("Fixing vdW terms...")
radius = ion_atom.rmin
epsilon = ion_atom.epsilon
parmed.tools.addLJType(base,":Na+",radius=radius,epsilon=epsilon).execute()
print("Saving files...")
base.write_parm("system.prmtop")
base.write_rst7("system.rst7")
print("All done!")
