#July 2017 Stefano Bosisio
#Compute single point interaction energy between dipole and charge


import os,re, sys, shutil
import math, random
import numpy as np
#benchmark the code
import time
from multiprocessing import Pool, Process
from functools import partial
import matplotlib.pyplot as plt
import seaborn as sbn
sbn.set_style("whitegrid")


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


def sodium_settings(system,r):
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
                cgnaidx = na_atom.cgAtomIdx() #detect is cgIdx
                na_charge = na_atom.property("charge").value()

                new_coords.set(cgnaidx,atom_coords) #set cgIdx and new coords into the AtomCoords object
            #outside the if-cycle edit sodium coordinates
            molecule = molecule.edit().setProperty("coordinates",new_coords).commit()
            #put sodium ion in the changed molecules group
            changedion.add(molecule)
    #update these changes to the system
    system.update(changedion)

    return atom_coords,na_charge

def sodium_check(system):
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


def dipole_charges(system):
    r"""This function returns the value of H1-H2 atoms charges

    Parameters
    ----------
    system:     Sire System
                System created from topology and coordinate files

    Returns
    ----------


    h1_charge:  float64
                Current charge for H1 atom

    h2_charge:  float64
                Current charge for H2 atom

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
            dipole = molecule
            for x in range(0,molnatoms):
                atom = molatoms[x]
                if atom.name().value()=="H1":
                    h1_atom = atom
                    h1_charge = h1_atom.property("charge").value()
                elif atom.name().value()=="H2":
                    h2_atom = atom
                    h2_charge = h2_atom.property("charge").value()
                else:
                    pass

    return h1_charge,h2_charge

def dipole_bond(system):
    r"""This function return the dipole bond length

    Parameters
    ----------
    system:     Sire System
                System created from topology and coordinate files



    Returns
    ----------

    bond:       float64
                Dipole bond length

    """

    molnums = system.molNums()
    space = system.property("space")
    for molnum in molnums:
        molecule = system.molecule(molnum).molecule()
        molname = str(molecule.residue().name().value())
        molatoms = molecule.atoms()
        molnatoms = molecule.nAtoms()
        if molname=="LIG":
            new_coords = AtomCoords(molecule.property("coordinates"))
            dipole = molecule
            for x in range(0,molnatoms):
                atom = molatoms[x]
                if atom.name().value()=="H1":
                    h1_atom = atom
                    h1_coords = h1_atom.property("coordinates")            

                elif atom.name().value()=="H2":
                    h2_atom = atom
                    h2_coords = h2_atom.property("coordinates")                

                else:
                    pass
            #Compute the original bond length
            bond = space.calcDist(h1_coords,h2_coords)
            

    return bond


def dipole_scaling(system,scaling_factor,bondlength=None,h1_orig=None,h2_orig=None):
    r"""This function set the  dipole position onto the x-axis. The bond is scaled
    by bond_scale and the angle w.r.t the charge is given by theta

    Parameters
    ----------
    system:     Sire System
                System created from topology and coordinate files

    scaling_factor: float64
                scaling factor for bond and charges, to find the best agreement
                with Coulombic interactions

    bondlength: float64
                Original bond length to scale


    h1_orig:    float64
                Original H1 charge, None if scaling_charges = False

    h2_orig:    float64
                Original H2 charge, None if scaling_charges = False


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
            dipole = molecule
            for x in range(0,molnatoms):
                atom = molatoms[x]
                if atom.name().value()=="H1":
                    h1_atom = atom
                    h1_coords = h1_atom.property("coordinates")
                    h1_idx = x
                    cgh1 = h1_atom.cgAtomIdx()

                elif atom.name().value()=="H2":
                    h2_atom = atom
                    h2_coords = h2_atom.property("coordinates")
                    h2_idx = x
                    cgh2 = h2_atom.cgAtomIdx()

                else:
                    pass

            #set the new bond length and charges

            print("Scaling bond length")
            bond = bondlength/scaling_factor
            

            print("Scaling charges")
            h1_charge = h1_orig*scaling_factor
            h2_charge = h2_orig*scaling_factor
            #commit the change in charges
            dipole = dipole.edit().atom(AtomIdx(h1_idx)).setProperty("charge",h1_charge*mod_electron).molecule().commit()
            dipole = dipole.edit().atom(AtomIdx(h2_idx)).setProperty("charge",h2_charge*mod_electron).molecule().commit()


            #now re locate the dipole onto the y-axis with the new distance
            h1_newcoords = Vector(bond/2.0,0.0,0.0)
            h2_newcoords = Vector(-bond/2.0,0.0,0.0)
            print("New Dipole Coords")
            print(h1_newcoords)
            print(h2_newcoords)
            #uncomment below to visualize the output
            #time.sleep(2)
            new_coords.set(cgh1,h1_newcoords)
            new_coords.set(cgh2,h2_newcoords)
            dipole = dipole.edit().setProperty("coordinates",new_coords).commit()

            #add to the modification group
            changedorigin.add(dipole)
    #update everything
    system.update(changedorigin)


def dipole_rotate(system,theta):
    r"""This function rotate the dipole around the z-axis by an angle theta.
    The rotation computed through a rotation matrix

    Parameters
    ----------
    system:     Sire System
                System created from topology and coordinate files

    theta:      float64
                Angle of rotation


    Returns
    ----------

    """
    molnums = system.molNums()

    for molnum in molnums:
        molecule = system.molecule(molnum).molecule()
        molname = str(molecule.residue().name().value())
        molatoms = molecule.atoms()
        molnatoms = molecule.nAtoms()
        if molname=="LIG":
            new_coords = AtomCoords(molecule.property("coordinates"))
            dipole = molecule
            for x in range(0,molnatoms):
                atom = molatoms[x]
                if atom.name().value()=="H1":
                    h1_atom = atom
                    h1_coords = h1_atom.property("coordinates")
                    h1_idx = x
                    cgh1 = h1_atom.cgAtomIdx()

                elif atom.name().value()=="H2":
                    h2_atom = atom
                    h2_coords = h2_atom.property("coordinates")
                    h2_idx = x
                    cgh2 = h2_atom.cgAtomIdx()

                else:
                    pass


    #now rotate the dipole accordign to theta
    #compute the matrix elements
    rot00 = cos(theta)
    rot01 = -sin(theta)
    rot02 = 0
    rot10 = sin(theta)
    rot11 = cos(theta)
    rot12 = 0
    rot20 = 0
    rot21 = 0
    rot22 = 1

    rotmat = Matrix(rot00,rot01,rot02,\
                    rot10,rot11,rot12,\
                    rot20,rot21,rot22)
    newh1 = rotmat*h1_coords
    newh2 = rotmat*h2_coords



    return newh1,newh2


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
            h1_charge = molecule.atom(AtomName("H1")).property("charge")
            h2_charge = molecule.atom(AtomName("H2")).property("charge")
            print("Actual charges for H1 and H2 atoms...")
            print(h1_charge,h2_charge)


def pythag(coords1,coords2):
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
    dist =math.sqrt((coords1[0]-coords2[0])**2 + (coords1[1]-coords2[1])**2 + (coords1[2]-coords2[2])**2)
    return dist

def Coulomb(qion,h1,h2,r1,r2):

    pot_ion1 = one_over_four_pi_eps0*(qion*h1)*( (1/r1))*(1/82.0)
    pot_ion2 = one_over_four_pi_eps0*(qion*h2)*( (1/r2))*(1/82.0)
    pot = pot_ion1 + pot_ion2

    return pot

def theory_rotational_average(distances,outputfile):
    r"""Function to compute the rotationally average interaction
    in vacuum

    Parameters
    ----------
    distances:   list
                 list of distances between ion and dipole
    theta:       float64
                 Dipole's angle we are computing the interaction
    outputfile:  file
                 file to write out distances and energies

    Returns
    ----------
    nrgs:        list
                 list of computed energies at theta for all the distances
    """

    q1, mu, epsw = 1, (2.35*0.20819434), 82.0 #  charge q1 = 1e , dipole mu = 0.47 eA, dielectric TIP3P water epsw = 82.0
    kT = -(1/(3*0.592186))
    nrgs = []
    for r in distances:
        qdip =  (one_over_four_pi_eps0)*((q1*mu)**2)*(r**-4)*(1/epsw**2)*(kT)
        nrgs.append(qdip)
        outputfile.write("%.4f,%.4f\n" % (r,qdip))

    return nrgs

########################################MAIN#################################################################


top_file = sys.argv[1]  #topology
crd_file = sys.argv[2]  #coordiantes
#create the Sire System
amber = Amber()
molecules, space = amber.readCrdTop(crd_file, top_file)
system = createSystem(molecules)
#create a force field to have a reference for Coulombic interactions, it may be useful
system, energiesDict = setupForcefields(system, space,cutoff_type.val,cutoff_dist,rf_dielectric.val)
#create a vector of theta angles we want to test and bond distances to fix the dipole-coulombic agreement
thetas = np.linspace(0,2*pi,13) #go 30deg by 30deg
rot_angle = 30*(pi/180.0)#thetas[1] #rotation angle
n_rotations=13 #total number of rotations
r_vector = np.linspace(5,200,50) #vector of distances for the ion
bond_scale = 2# max factor we are scaling bonds and charges (20)
#Wjat we need to do:
#1) Take down the dipole charges - original --> dipole_charges
#2) The bond length --> dipole_bond
#3) Set the dipole on the x-axis in order to have a theta =0.0  --> dipole_scaling
#4) scale if necessary boh charges and bond --> dipole_scale
#5) rotate the dipole --> dipole_rotate
h1_orig, h2_orig = dipole_charges(system) #charges
h1h2_bond= dipole_bond(system) #bond
#set the dipole onto the y axis
dipole_scaling(system,1.0,h1h2_bond,h1_orig,h2_orig) #if the scaling factor =1.0 
#the dipole will remain the same

#output parameters
#outputresults directory
output = "results/"
if not os.path.exists(output):
  os.makedirs(output)
#####

#store all the results in two dictionaries - they may be useful 
theory = {} #theoretical behavior
computed = {} #computed one

#cycle 
print("Computing charge-dipole rotationally average interaction")

outputfile = "results/theory.csv" 
ofile = open(outputfile,"w")

nrgs = theory_rotational_average(r_vector,ofile)
ofile.close()



print("Computing charge-dipole interactions with Coulombic scheme...")


avg = 0.0

for r in r_vector:
    print("distance %.4f" % r)
    counter = 0.0
    for theta in thetas:
        #rotate the dipole and comput ethe energy at each theta
        
        #print("Theta %.4f" % theta)
        outputfile =  "results/computed.csv" 
        ofile = open(outputfile,"w")
        h1_coords,h2_coords = dipole_rotate(system,theta)
        ion_coords,ion_charge = sodium_settings(system,r)
        h1_charge,h2_charge = dipole_charges(system)
        r1 = pythag(ion_coords,h1_coords)
        r2 = pythag(ion_coords,h2_coords)
        nrg = Coulomb(ion_charge,h1_charge,h2_charge,r1,r2)
        if counter==0:
            avg = nrg
            counter+=1
        else:
            avg = avg + (nrg - avg)/counter
            counter+=1
    print(avg)
    ofile.write("%.4f,%.16f\n" %(r,avg))


