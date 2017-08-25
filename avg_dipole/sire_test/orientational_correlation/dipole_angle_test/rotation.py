#July 2017 Stefano Bosisio
#Compute single theta values for charge-water interactions

import os,re, sys, shutil
import math, random
import numpy as np
#benchmark the code
import time

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

def water_info(system):
    r"""Function to retrieve the current position (initial) of the
    water molecule, the oxygen MUST stay in the origin

    Parameters
    ----------
    system:     Sire System
                System created from topology and coordinate files


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
        if molname=="WAT":
            water = molecule
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
                elif atom.name().value()=="O":
                    o_atom = atom
                    o_coords = o_atom.property("coordinates")
                    o_charge = o_atom.property("charge").value()

                else:
                    pass

    #compute the angle between h1-o-h2
    angle = space.calcAngle(h1_coords,o_coords,h2_coords)

    rad_angle = angle.value()
    angle_half = rad_angle/2.0

    #set h1 and h2 in order to have the dipole at 0 degrees
    h1_x = h1_coords[0]*cos(-angle_half) + h1_coords[1]*(-sin(-angle_half)) + h1_coords[2]*0.0
    h1_y = h1_coords[0]*sin(angle_half) + h1_coords[1]*( cos(-angle_half)) + h1_coords[2]*0.0
    h1_z = h1_coords[2]*1.0

    h2_x = h2_coords[0]*cos(-angle_half) + h2_coords[1]*(-sin(-angle_half)) + h2_coords[2]*0.0
    h2_y = h2_coords[0]*sin(-angle_half) + h2_coords[1]*( cos(-angle_half)) + h2_coords[2]*0.0
    h2_z = h2_coords[2]*1.0

    h1_new = Vector(h1_x,h1_y,h1_z)
    h2_new = Vector(h2_x,-h2_y,h2_z)
    new_coords = AtomCoords(water.property("coordinates"))
    new_coords.set(cgh1,h1_new)
    new_coords.set(cgh2,h2_new)
    water.edit().setProperty("coordinates",new_coords).commit()
    changedorigin.add(water)
    system.update(changedorigin)

    return h1_new,h2_new,h1charge,h2charge,o_charge,rad_angle


def water_check(system):
    r"""Sanity check function on water coordinates

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
    dist:   float64
        Distance between First-Second set of coordinates
    """
    r = math.sqrt((coords1[0]-coords2[0])**2 + (coords1[1]-coords2[1])**2 +(coords1[2]-coords2[2])**2 )

    return r

#@jit(float64(float64,float64,float64,float64,float64,float64,float64),nopython=True,nogil=True)
def BarkerWatts(qion,oxy,h1,h2,rO,r1,r2,cutoff,dielectric):
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
    pot_ionH1 = one_over_four_pi_eps0*(qion*h1)*( (1/r1) + krf*r1**2)# -crf)
    pot_ionH2 = one_over_four_pi_eps0*(qion*h2)*( (1/r2) + krf*r2**2)# -crf)
    pot_ionO  = one_over_four_pi_eps0*(qion*oxy)*((1/rO) + krf*r2**2)
    pot = pot_ionH1 + pot_ionH2 + pot_ionO

    return pot

def derived_RF(qion,r,theta,cutoff,dielectric):
	#in this function for each distance compute the derive charge-dipole Reaction field for a specific theta

    mu = (0.834)*(0.586) #hardcoding
    
    charge_part = -one_over_four_pi_eps0*qion*mu*cos(theta)
    rf_part = (1/r**2 - 2*(r/cutoff**3)*((dielectric-1)/(1+2*dielectric)) + (1/r**2)*(3/(1+2*dielectric)))

    nrg = charge_part*rf_part
    return nrg


def generateCoordFiles(system,idx):
    r"""Function to save mdcrd coordinate files, in order to create a trajectory
    file to visualize rotations

    Parameters
    ----------
    system:     Sire System
                System created from topology and coordinate files

    idx:        int
                A numerical index to save the i-th structure
    Returns
    ----------

    """
    #here we are genereting coordiantes file to visualize the situation
    tot_atoms = 4
    wat_coords = system[MGName("all")].moleculeAt(0).molecule().property("coordinates")
    na_coords = system[MGName("all")].moleculeAt(1).molecule().property("coordinates")
    coordinates = []
    for coord in wat_coords.toVector():
        coordinates.append(coord)
    for coord in na_coords.toVector():
        coordinates.append(coord)
    #print(coordinates)


    string_coord = str(coordinates)
    string_rev   = re.sub("[^0-9.e-]", " ", string_coord)
    string_split = string_rev.split()
    directory = "coordinates"

    if not os.path.exists(directory):
        os.makedirs(directory)

    file_name = directory + "/structure_%03d.mdcrd" %(idx)
    singlecoord = open(file_name, "w")
    singlecoord.write("TP3\n")
    singlecoord.write("    %s\n" %tot_atoms)
    index=0
    #a bit of rules to correctly write coordinates
    for f in string_split:
     if float(f)<0:
         singlecoord.write("  %.7f" % float(f))
     elif float(f)>=10:
         singlecoord.write("  %.7f" % float(f))
     else:
         singlecoord.write("   %.7f" % float(f))
     index+=1
     if not index%6:
         singlecoord.write("\n")




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
print("Retrieving water atoms coordinates...")
#First retrieve the coords of the atoms of water
h1coords,h2coords,h1charge,h2charge,o_charge,wat_angle=water_info(system)
print(h1coords,h2coords)
print("Water molecule set to have dipole degree = 0.0")
angle_shift = wat_angle/2.0
print("Set the ion onto the x-axis...")
#set the ion onto the x axis
ioncoord,qion=ion_set_x(system,3.0) #fix the dipole at an initial distance = 3.0
#define all the distance for the ion
r_vector = np.linspace(3,30,100) #
rot_rad = 15*pi/180.0
n_rotations = 2*pi/rot_rad #cover all the curves rom 0 to 360 step 3 deg
n_rotations= int(n_rotations)
#cycle through each theta angle, angle of rotation of the water w.r.t the z-axis
#cycle through the r position and save the file in the theta folder

counter = 0 
for i in range(0,n_rotations+1):

    
#rotate 180 for convection
    theta = i*rot_rad
    theta_deg = theta*(180/pi)

    #stating from th ecoordinates, compute COM
    oxy_coords = [0,0,0] # always in the origin
    num_comx = (oxy_coords[0]*16.0  + h1coords[0]*1.008 + h2coords[0]*1.008)
    num_comy = (oxy_coords[1]*16.0  + h1coords[1]*1.008 + h2coords[1]*1.008)
    num_comz = (oxy_coords[2]*16.0  + h1coords[2]*1.008 + h2coords[2]*1.008)
    den_com = 16.0+1.008+1.008
    COM =[ num_comx/den_com, num_comy/den_com, num_comz/den_com]
    
    #compute the angle between  ion - COM- oxygen, this should be equal to theta
    ion_com = pythag(ioncoord,COM)
    ion_oxy = pythag(ioncoord,oxy_coords)
    oxy_com = pythag(oxy_coords,COM)
    numerator =oxy_com**2 + ion_com**2 - ion_oxy**2
    denomin   = 2*oxy_com*ion_oxy
    ratio = numerator/denomin
    if ratio > 1.00:
        #print(ratio)
        ratio = 1.0
    
    angle = math.acos(ratio) 
    #fa = math.sqrt(ioncoord[0]**2 + ioncoord[1]**2 + ioncoord[2]**2)
    #fb = math.sqrt(COM[0]**2 + COM[1]**2 + COM[2]**2)
    #fdot = ioncoord[0]*COM[0] + ioncoord[1]*COM[1] + ioncoord[2]*COM[2]
    #dotangle = math.acos(fdot/(fa*fb))
    fcrossX = ioncoord[1]*COM[2] - ioncoord[2]*COM[1]
    fcrossY = ioncoord[2]*COM[0] - ioncoord[0]*COM[2]
    fcrossZ = ioncoord[0]*COM[1] - ioncoord[1]*COM[0]
    fcross  =math.sqrt(fcrossX *fcrossX + fcrossY*fcrossY + fcrossZ*fcrossZ)
    fdot = ioncoord[0]*COM[0] + ioncoord[1]*COM[1] + ioncoord[2]*COM[2]
    dotangle = math.atan2(fcross,fdot)
    print("Carnot angle : %.4f vs  dotangle %.4f Theta %.4f" % (angle*57.29,dotangle*57.29,theta*57.29) )


    #rotate the hydrogens...
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

'''
    molnums = system.molNums()
    changedorigin = MoleculeGroup("changedorigin")
    space = system.property("space")
    for molnum in molnums:
        molecule = system.molecule(molnum).molecule()
        molname = str(molecule.residue().name().value())
        molatoms = molecule.atoms()
        molnatoms = molecule.nAtoms()
        if molname=="WAT":
            new_coords = AtomCoords(molecule.property("coordinates"))
            #bond lenght  = 0.586
            for x in range(0,molnatoms):
                atom = molatoms[x]
                if atom.name().value()=="H1":
                    h1_atom = atom
                    cgh1 = h1_atom.cgAtomIdx()
                elif atom.name().value()=="H2":
                    h2_atom = atom
                    cgh2 = h2_atom.cgAtomIdx()
                else:
                    pass
            
            #set the new coords
            new_coords.set(cgh1,Vector(h1coords))
            new_coords.set(cgh2,Vector(h2coords))
            molecule = molecule.edit().setProperty("coordinates",new_coords).commit()
            changedorigin.add(molecule)

    system.update(changedorigin)
    generateCoordFiles(system,counter)

'''
