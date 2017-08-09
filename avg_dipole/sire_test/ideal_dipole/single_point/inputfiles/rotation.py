#July 2017 Stefano Bosisio
#Compute single theta values for charge-dipol interactions

import os,re, sys, shutil
import math, random
import numpy as np
#benchmark the code
import time
from multiprocessing import Pool, Process
from functools import partial
import numba
from numba import jit, float64,int64
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D

from Sire.IO import *
from Sire.Mol import *
from Sire.CAS import *
from Sire.System import *
from Sire.Move import *
from Sire.MM import *
from Sire.FF import *
from Sire.Units import *
from Sire.Vol import *
from Sire.Maths import *
from Sire.Base import *
from Sire.Qt import *
from Sire.ID import *
from Sire.Config import *
from Sire.Analysis import *
from Sire.Tools.DCDFile import *
from Sire.Tools import Parameter, resolveParameters
import Sire.Stream

rf_dielectric = Parameter("reaction field dielectric", 82.0,
                          """Dielectric constant to use if the reaction field cutoff method is used.""")

cutoff_type = Parameter("cutoff type", "cutoffperiodic", """The cutoff method to use during the simulation.""")

cutoff_dist = Parameter("cutoff distance", 30 * angstrom,
                        """The cutoff distance to use for the non-bonded interactions.""")
combining_rules = Parameter("combining rules", "arithmetic",
                            """Combining rules to use for the non-bonded interactions.""")

#global variables
global pi
pi = math.pi
cos = math.cos
sin = math.sin #usage sin(x)
#############

def createSystem(molecules):

    print("Creating the system...")

    moleculeNumbers = molecules.molNums()
    moleculeList = []

    for moleculeNumber in moleculeNumbers:
        molecule = molecules.molecule(moleculeNumber).molecule()
        moleculeList.append(molecule)

    molecules = MoleculeGroup("molecules")
    ions = MoleculeGroup("ions")

    for molecule in moleculeList:
        natoms = molecule.nAtoms()
        if natoms == 1:
            ions.add(molecule)
        else:
            molecules.add(molecule)

    all = MoleculeGroup("all")
    all.add(molecules)
    all.add(ions)

    # Add these groups to the System
    system = System()

    system.add(all)
    system.add(molecules)
    system.add(ions)

    return system



def setupForcefields(system, space, cutoff_type,cutoff_dist,rf_dielectric):

    print("Creating force fields... ")

    all = system[MGName("all")]
    molecules = system[MGName("molecules")]
    ions = system[MGName("ions")]
    energiesDict = {}

    inter_ions_molecules_nonbondedff = InterGroupCLJFF("ions:molecules")
    energiesDict["E_{ions:molecules}^{coulomb}"] = [ 0.0, 0]

    if (cutoff_type != "nocutoff"):
        inter_ions_molecules_nonbondedff.setUseReactionField(True)
        inter_ions_molecules_nonbondedff.setReactionFieldDielectric(rf_dielectric)

    inter_ions_molecules_nonbondedff.add(ions, MGIdx(0))
    inter_ions_molecules_nonbondedff.add(molecules, MGIdx(1))


    # Here is the list of all forcefields
    forcefields = [inter_ions_molecules_nonbondedff]

    for forcefield in forcefields:
        system.add(forcefield)

    system.setProperty("space", space)
    system.setProperty("switchingFunction", CHARMMSwitchingFunction(cutoff_dist.val))
    system.setProperty("combiningRules", VariantProperty(combining_rules.val))

    total_nrg = inter_ions_molecules_nonbondedff.components().total()
    e_total = system.totalComponent()

    system.setComponent(e_total, total_nrg)

    return system,energiesDict

def ion_set_x(system,r):
    r"""This function set the sodium at a distance r from the origin (dipole middle point)

    Parameters
    ----------
    system:     Sire System
                System created from topology and coordinate files

    r:          float64
                distance between the ion and dipole middle point
    Returns
    ----------
    """
    #a bit tricky to change coordinates
    molnums = system.molNums()
    #create a changedion group, to store all the ion coords we are changing
    changedion = MoleculeGroup("changedion")
    for molnum in molnums:
        #gather all the info about the molecule
        molecule = system.molecule(molnum).molecule()
        molname = str(molecule.residue().name().value())
        molatoms = molecule.atoms()
        molnatoms = molecule.nAtoms()
        #if the molecule is the sodium ion set the coords
        if molname == "Na+":
            #print("Setting Na+ coordinates...")
            #create an AtomCoords object
            new_coords = AtomCoords(molecule.property("coordinates"))
            for x in range(0,molnatoms):
                #give the new coords for sodium ion
                atom_coords = Vector(r,0,0)
                na_atom = molatoms[x] #take the sodium ion
                ion_charge = na_atom.property("charge").value()
                cgnaidx = na_atom.cgAtomIdx() #detect is cgIdx
                new_coords.set(cgnaidx,atom_coords) #set cgIdx and new coords into the AtomCoords object
            #outside the if-cycle edit sodium coordinates
            molecule = molecule.edit().setProperty("coordinates",new_coords).commit()
            #put sodium ion in the changed molecules group
            changedion.add(molecule)
    #update these changes to the system
    system.update(changedion)

    return atom_coords,ion_charge

def ion_check(system):
    r"""Sanity check function on ion coordinates, to see if the ion is in the
    desired position

    Parameters
    ----------
    system:     Sire System
                System created from topology and coordinate files

    Returns
    ----------
    """
    molnums = system.molNums()

    for molnum in molnums:
        molecule = system.molecule(molnum).molecule()
        molname = str(molecule.residue().name().value())
        molnatoms =molecule.nAtoms()
        if molname=="Na+":
            coords=molecule.atom(AtomName("Na+")).property("coordinates")
            print("Actual coordinates for Na+ ion...")
            print(coords)

def dipole_set_x(system,scale=None):
    r"""Function to fix the dipole onto the x axis

    Parameters
    ----------
    system:     Sire System
                System created from topology and coordinate files
    scale:		float
    			Bond length we want


    Returns
    ----------
    """

    molnums = system.molNums()
    changedorigin = MoleculeGroup("changedorigin")
    space = system.property("space")
    for molnum in molnums:
        molecule = system.molecule(molnum).molecule()
        molname = str(molecule.residue().name().value())
        molatoms = molecule.atoms()
        molnatoms = molecule.nAtoms()
        if molname=="LIG":
            new_coords = AtomCoords(molecule.property("coordinates"))
            print("Setting dipole middle point in the origin...")
            #bond lenght  = 0.586
            for x in range(0,molnatoms):
                atom = molatoms[x]
                if atom.name().value()=="H1":
                    h1_atom = atom
                    h1_coords = h1_atom.property("coordinates")
                    h1charge = h1_atom.property("charge").value()
                    cgh1 = h1_atom.cgAtomIdx()
                elif atom.name().value()=="H2":
                    h2_atom = atom
                    h2_coords = h2_atom.property("coordinates")
                    h2charge = h2_atom.property("charge").value()
                    cgh2 = h2_atom.cgAtomIdx()
                else:
                    pass
            #now compute the distance
            bond = space.calcDist(h1_coords,h2_coords)
            #from the initial bond start the scaling:

            if scale is not None :
            	bond = scale
            	h1_newcoords=Vector(bond/2.0,0,0)
            	h2_newcoords = Vector(-bond/2.0,0,0)



            else:
				#set the new coords
                h1_newcoords = Vector(bond/2.0,0,0)
                h2_newcoords = Vector(-bond/2.0,0,0)

            #set the new coords
            new_coords.set(cgh1,h1_newcoords)
            new_coords.set(cgh2,h2_newcoords)
            molecule = molecule.edit().setProperty("coordinates",new_coords).commit()
            changedorigin.add(molecule)

    system.update(changedorigin)
    #return the bond length
    return bond,h1_newcoords,h2_newcoords,h1charge,h2charge

def dipole_check(system):
    r"""Sanity check function on dipole coordinates

    Parameters
    ----------
    system:     Sire System
                System created from topology and coordinate files

    Returns
    ----------
    """
    molnums = system.molNums()

    for molnum in molnums:
        molecule = system.molecule(molnum).molecule()
        molname = str(molecule.residue().name().value())
        molnatoms =molecule.nAtoms()
        if molname=="LIG":
            h1_coords=molecule.atom(AtomName("H1")).property("coordinates")
            h2_coords=molecule.atom(AtomName("H2")).property("coordinates")
            print("Actual coordinates for dipole atoms...")
            print(h1_coords,h2_coords)





def pythag(coords1,coords2,angle,bond,neg=False):
    r"""Function to compute the distance between two set of coordinates

    Parameters
    ----------
    coord1: 1D array
        First set of coordinates (x,y,z)

    coord2: 1D array
        Second set of coordinates

    Returns
    ----------
    dist:   flaot64
        Distance between First-Second set of coordinates
    """
    if neg:
        r = math.sqrt((coords1[0]-coords2[0])**2 + (coords1[1]-coords2[1])**2 +(coords1[2]-coords2[2])**2 )
        a_fact = 0.5*bond*cos(angle)
        dist = r-a_fact
    else:
        r = math.sqrt((coords1[0]-coords2[0])**2 + (coords1[1]-coords2[1])**2 +(coords1[2]-coords2[2])**2 )
        a_fact = 0.5*bond*cos(angle)
        dist = r+a_fact

    return dist
#@jit(float64(float64,float64,float64,float64,float64,float64,float64),nopython=True,nogil=True)
def BarkerWatts(qion,h1,h2,r1,r2,cutoff,dielectric):
    r"""Function to compute BW reaction field Coulombic interaction potential

    Parameters
    ----------
    qion:   float64
            ion charge ( +1.0 e)
    h1:     float64
            first dipole atom charge
    h2:     float64
            second dipole atom charge
    r1:     float64
            distance h1-ion
    r2:     float64
            distance h2-ion
    cutoff: float64
            BW radius of cutoff in A
    dielectric: float64
            BW dielectric constant
    Returns
    ----------
    pot:   float64
        Coulombic potential between ion and h1 and h2
    """
    #function to compute the coulombic energy with BW without shift
    #BW = sum_i (one_over_four_pi_eps0)* (q_ion*q_i)* (1/r_ion-i + krf*r_ion-i**2)

    krf = (1/cutoff**3)*((dielectric-1)/(2*dielectric + 1))
    #crf = (1/cutoff)*(3*dielectric)/(2*dielectric + 1)  #crf to compare with Sire results
    pot_ion1 = one_over_four_pi_eps0*(qion*h1)*( (1/r1) + krf*r1**2)# -crf)
    pot_ion2 = one_over_four_pi_eps0*(qion*h2)*( (1/r2) + krf*r2**2)# -crf)
    pot = pot_ion1 + pot_ion2

    return pot

def derived_RF(qion,dipolecharge,dipolebond,theta,r,cutoff,dielectric):
	#in this function for each distance compute the derive charge-dipole Reaction field for a specific theta

	mu = abs(dipolecharge)*dipolebond

	charge_part = -one_over_four_pi_eps0*qion*mu*cos(theta)
	rf_part = (1/r**2 - 2*(r/cutoff**3)*((dielectric-1)/(1+2*dielectric)) + (1/r**2)*(3/(1+2*dielectric)))

	nrg = charge_part*rf_part
	return nrg


########################################MAIN#################################################################

top_file = sys.argv[1]  #topology
crd_file = sys.argv[2]  #coordiantes
#create the Sire System
amber = Amber()
molecules, space = amber.readCrdTop(crd_file, top_file)
system = createSystem(molecules)
#create a force field to have a reference for Coulombic interactions, it may be useful
system, energiesDict = setupForcefields(system, space,cutoff_type.val,cutoff_dist,rf_dielectric.val)
#algorithm:
#1 Fix the dipol onto the x axis (b/2,0,0) (-b/2,0,0)
#2 fix the dipoole at a distnace from the dipole, start -bond?
#3 compute the interaction with the BWRF with distance charge-atom1 and charge-atom2 as r -= 0.5*l*cos(theta)
#all above must be done in a cycle and modifying the  ion from start to 30 Angstrom
#4 save the val on a file open as append
#5 create the new directory by rotating the dipol around y-axis by 15/30 deg
print("First single point calculation")
print("Fixing dipole onto x-axis")
bond,h1coords,h2coords,h1charge,h2charge=dipole_set_x(system) #bond is the dipole bond length
#dipole_line = (1,0,0) #initially te dipole isonto the x axis
middle_point = (0,0,0) #Check if the middle point is in the origin always
dipole_check(system) #sanity check if the coords are all right
#set the ion onto the x axis
ioncoord,qion=ion_set_x(system,3.0) #fix the dipole at an initial distance = 3.0
#ion_line = (1,0,0) # the ion will always be onto the x line
#define all the distance for the ion
r_vector = np.linspace(3,30,100) #
rot_rad = 3*pi/180.0
n_rotations = 120 #cover all the curves rom 0 to 360 step 3 deg
#thetas = np.linspace(0,2*pi,24) #rotation of 15 degrees

#dipole bond lengths scaled up and down by a factor:
#5,4,3,2.5 ,2, 1.95,1.5,1.3,1.2,1,1.1,1.2,1.3,1.5,1.95,2,2.5,3,4,5
#scales = [0.1172,0.1465,0.1953,0.234,0.293,0.300,0.390,0.450,0.488,0.586,0.644,0.703,0.761,0.879,1.142,1.172,1.465,1.758,2.344,2.929]
scales = np.linspace(0.05,3.0,100)
#cycle through the bonds
#cycle through each theta angle, defining the position of the dipole w.r.t the dipole
#cycle through the r position and save the file in the theta folder
for scale in scales:

    print("Dipole bond will be set at:%.4f A" %scale)
    #scale the dipole bond length
    bond,h1coords,h2coords,h1charge,h2charge=dipole_set_x(system,scale)
    #sanity check on the dipole
    dipole_check(system)
    bonddir = "l_%.3f" % bond
    if not os.path.exists(bonddir):
        os.makedirs(bonddir)
    print("Created directory %s" % bonddir)

    for i in range(0,n_rotations+1):

    	theta = i*rot_rad
    	theta_deg = theta*(180/pi)
    	#create a new theta directory
    	outputfold = "%s/%d_deg" % (bonddir,math.ceil(theta_deg))
    	if not os.path.exists(outputfold):
    		os.makedirs(outputfold)

    	computefile = outputfold + "/computed.dat"
    	computed = open(computefile,"w")
    	theoryfile = outputfold + "/theory.dat"
    	theory = open(theoryfile,"w")

    	print("Angle between ion and dipole %.2f" % theta_deg)
    	#initially we don't need to rotate the dipole, then rotate the initial h1coords/h2coords
    	for r in r_vector:
    		#print("Setting ion coordinates for distance %.4f" % r)
    		ioncoord,qion = ion_set_x(system,r)
    		#the angle between ion and dipole is given by theta, which is the rotation around the y axis fro the dipole
    		r1 = pythag(ioncoord,middle_point,theta,bond,neg=True) #negative charge
    		r2 = pythag(ioncoord,middle_point,theta,bond,neg=False) #positive charge
    		#Barker Watts Coulombic energy
    		nrg = BarkerWatts(qion,h1charge,h2charge,r1,r2,cutoff_dist.val.value(),rf_dielectric.val)

    		rf_nrg = derived_RF(qion,h1charge,bond,theta,r,cutoff_dist.val.value(),rf_dielectric.val)
    		#print("r,ioncoord,h1coords,h2coords,qion,h1charge,h2charge,r1,r2")
    		#print(r,theta,ioncoord,h1coords,h2coords,qion,h1charge,h2charge,r1,r2)
    		#time.sleep(2)
    		#write on the file
    		computed.write("%.4f,%.8f\n" % (r,nrg))
    		theory.write("%.4f,%.8f\n" %(r,rf_nrg))


    	#close the file
    	computed.close()
    	theory.close()
    	#once all the distances are computed and save modify the dipole in theta
    	#hence, rotatino of the dipole coords around the y axis with negative angles
    	#create the matrix
    	r00 = cos(-rot_rad)
    	r01 = -sin(-rot_rad)
    	r10 = sin(-rot_rad)
    	r11 = cos(-rot_rad)
    	r22 = 1
    	#compute the newh1 and h2 position
    	h1coords = (h1coords[0]*r00 + h1coords[1]*r01,\
    				h1coords[0]*r10 + h1coords[1]*r11,\
    				h1coords[2]*r22)
    	h2coords = (h2coords[0]*r00 + h2coords[1]*r01,\
    				h2coords[0]*r10 + h2coords[1]*r11,\
    				h2coords[2]*r22)
