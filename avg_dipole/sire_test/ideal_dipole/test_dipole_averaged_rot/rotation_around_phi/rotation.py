#June 2017 Stefano Bosisio
#Script to compute Coulombic interaction betweenan ideal dipole and an ion
#the dipole is rotated along all the possible rotaitonal degree of freedom (theta, phi, psi)
#Usage: ~/sire.app/bin/python rotation.py system.prmtop system.rst7

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

def sodium_coordinates(system,r):
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
            print("Setting Na+ coordinates...")
            #create an AtomCoords object
            new_coords = AtomCoords(molecule.property("coordinates"))
            for x in range(0,molnatoms):
                #give the new coords for sodium ion
                atom_coords = Vector(0,0,r)
                na_atom = molatoms[x] #take the sodium ion
                cgnaidx = na_atom.cgAtomIdx() #detect is cgIdx
                new_coords.set(cgnaidx,atom_coords) #set cgIdx and new coords into the AtomCoords object
            #outside the if-cycle edit sodium coordinates
            molecule = molecule.edit().setProperty("coordinates",new_coords).commit()
            #put sodium ion in the changed molecules group
            changedion.add(molecule)
    #update these changes to the system
    system.update(changedion)

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

def dipole_set_x(system):
    r"""This function is used only once and set the dipole bond middle point in the
    origin of a Cartesian system of reference, while the two atoms will lie on
    +b/2, 0, 0 and -b/2, 0, 0 respectively, where b is the dipole bond length

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
        if molname=="LIG":
            new_coords = AtomCoords(molecule.property("coordinates"))
            print("Setting dipole middle point in the origin...")
            #bond lenght  = 0.586
            for x in range(0,molnatoms):
                atom = molatoms[x]
                if atom.name().value()=="H1":
                    h1_atom = atom
                    h1_coords = h1_atom.property("coordinates")
                    cgh1 = h1_atom.cgAtomIdx()
                elif atom.name().value()=="H2":
                    h2_atom = atom
                    h2_coords = h2_atom.property("coordinates")
                    cgh2 = h2_atom.cgAtomIdx()
                else:
                    pass
            #now compute the distance
            bond = space.calcDist(h1_coords,h2_coords)
            #set the new coords
            h1_newcoords = Vector(0.0,0.0,bond/2.0)
            h2_newcoords = Vector(0.0,0.0,-bond/2.0)
            #set the new coords
            new_coords.set(cgh1,h1_newcoords)
            new_coords.set(cgh2,h2_newcoords)
            molecule = molecule.edit().setProperty("coordinates",new_coords).commit()
            changedorigin.add(molecule)

    system.update(changedorigin)

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




def rotation_matrix(system,crd_counter,dcd_creator=None):
    r"""Rotation of the ion around the z axis

    Parameters
    ----------
    system:     Sire System
                System created from topology and coordinate files
    dcd_creator:File
                file to write all the coordinates file, if saved, to create a
                final trajectory file (netcdf) to visualize the rotation


    Returns
    ----------

    nrgs:       1D array
                It stores all the Coulombic energies computed for each rotation
                from here the mean and std deviation are computed
    """

    molnums = system.molNums()
    nrgs = [] #list of energy values
    #preparation  - this routine takes just 0.0001 s
    natoms = 0
    for molnum in molnums:
        #gather all the info about the molecule
        molecule = system.molecule(molnum).molecule()
        molname = str(molecule.residue().name().value())
        molatoms = molecule.atoms()
        molnatoms = molecule.nAtoms()
        natoms+=molnatoms
        if molname=="LIG":
            dipole = molecule
            for at in range(0,molnatoms):
                atomx = molatoms[at]
                if atomx.name().value()=="H1":
                    h1 = atomx
                    cgh1 = atomx.cgAtomIdx()
                    h1idx = at  #create an index for hydrogen 1
                    h1charge = h1.property("charge").value()
                    init_h1coords = Vector(h1.property("coordinates"))
                elif atomx.name().value()=="H2":
                    h2 = atomx
                    cgh2 = atomx.cgAtomIdx()
                    h2idx = at
                    h2charge = h2.property("charge").value()
                    init_h2coords = Vector(h2.property("coordinates"))
                else:
                    continue
        else:
            #store ion coordiantes and charge
            ioncoord = molecule.atom(AtomName("Na+")).property("coordinates")
            qion = molecule.atom(AtomName("Na+")).property("charge").value()


    TWOPI = 2*pi
    #1) set the coord, if i==0 start on the z axis otherwise shift the coords as
    #h1_new = init_h1[0]8sin(rot),0,h1[0]*cos(rot)
    #rotate around z axis
    thetas = np.linspace(0,pi,50)
    #phis = np.linspace(0,2*pi,100)

    start =time.time()

    for theta in thetas:
        #print(theta*(180/pi))
        #init_h1coords = half of the dipole bond length
        h1_current = (init_h1coords[2]*sin(theta),0.0,init_h1coords[2]*cos(theta))
        h2_current = (init_h2coords[2]*sin(theta),0.0,init_h2coords[2]*cos(theta))

        #now rotate around the z-axis
        for phi in phis:

            r00 = cos(phi)
            r01 = -sin(phi)

            r10 = sin(phi)
            r11 = cos(phi)

            r22 = 1

            h1_new = (h1_current[0]*r00 + h1_current[1]*r01,\
            		  h1_current[0]*r10 + h1_current[1]*r11,\
            		  h1_current[2]*r22)

            h2_new = (h2_current[0]*r00 + h2_current[1]*r01,\
            		  h2_current[0]*r10 + h2_current[1]*r11,\
            		  h2_current[2]*r22)

            #print(h1_new,h2_new)
            #time.sleep(2)
            #compute the energy
            r1 = pythag(ioncoord,h1_new)
            r2 = pythag(ioncoord,h2_new)
            #Barker Watts attack!
            nrg = BarkerWatts(qion,h1charge,h2charge,r1,r2,cutoff_dist.val.value(),rf_dielectric.val)
            nrgs.append(nrg)



    end= time.time()
    #print(orientations)
    print("Total time %.4f" % (end-start))

    return nrgs


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

def  plot_3d(coords_list):
    r"""This function takes as input all the corods we have generated throught
    the rotations and plot them onto a 3D plane. To test the goodness of the rotation
    we should obtain a sphere

    Parameters
    ----------
    coords_list:	list
    				list with all teh new coords for both dipole atoms after each
    				rotation

    Returns
    ---------
    """
    outputf = open("coords_list.dat","w")

    for coord in coords_list:
        outputf.write("%.8f,%.8f,%.8f\n" % coord)

'''
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_aspect("equal")

    ax.scatter(*zip(*coords_list),marker="o",color="r")
    ax.set_xlabel("X axis",fontsize=20)
    ax.set_ylabel("Y axis",fontsize=20)
    ax.set_zlabel("Z axis",fontsize=20)

    fig.savefig("Full_rotation.png",dpi=300)
    outpf = open("coords_list.dat","w")
    for coords in coords_list:
        outpf.write(coords)
'''
def write_pdb(system,natoms,idx):
    r"""Function to save a pdb of dipole-ion system for each rotation

    Parameters
    ----------
    system:     Sire System
                System created from topology and coordinate files

    natoms:     int
                Total number of atoms in the Sire System

    idx:        int
                A numerical index to save the i-th structure

    Returns
    ----------

    """
    directory = "pdb"

    if not os.path.exists(directory):
        os.makedirs(directory)

    pdbfilename = "pdb/structure_%d.pdb" % idx
    pdbfile = open(pdbfilename,"w")

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
    #create a group of coords
    counter=0
    array = []
    for f in string_split:
        array.append(float(f))

    outline=""

    for f in range(0,natoms):
        if array[3*f]<0.0 and array[3*f+1]<0.0 and array[3*f+2]<0.0:
            coords = "      %.3f  %.3f  %.3f" %(array[3*f],array[3*f+1],array[3*f+2])
        elif array[3*f+1]<0.0 and array[3*f+2] <0.0:
            coords = "       %.3f  %.3f  %.3f" %(array[3*f],array[3*f+1],array[3*f+2])
        elif array[3*f+1]<0.0 and array[3*f+2] <0.0:
            coords = "       %.3f  %.3f  %.3f" %(array[3*f],array[3*f+1],array[3*f+2])
        elif array[3*f]<0.0 and array[3*f+2] <0.0:
            coords = "      %.3f   %.3f  %.3f" %(array[3*f],array[3*f+1],array[3*f+2])
        elif array[3*f]<0.0 and array[3*f+1] <0.0:
            coords =  "      %.3f  %.3f   %.3f" %(array[3*f],array[3*f+1],array[3*f+2])
        elif array[3*f+2] < 0.0:
            coords = "       %.3f   %.3f  %.3f" %(array[3*f],array[3*f+1],array[3*f+2])
        elif array[3*f+1] <0.0:
            coords = "       %.3f  %.3f   %.3f" %(array[3*f],array[3*f+1],array[3*f+2])
        elif array[3*f]<0.0:
            coords = "      %.3f   %.3f   %.3f" %(array[3*f],array[3*f+1],array[3*f+2])
        else:
            coords = "       %.3f   %.3f   %.3f"%(array[3*f],array[3*f+1],array[3*f+2])
        outline+="ATOM      %d  H1  LIG     1%s  1.00  0.00           O\n" % (f, coords)
    outline+="END"
    pdbfile.write(outline)
    pdbfile.close()
    #sys.exit(-1)

def generateCoordFiles(system,idx,axis):
    r"""Function to save mdcrd coordinate files, in order to create a trajectory
    file to visualize rotations

    Parameters
    ----------
    system:     Sire System
                System created from topology and coordinate files

    idx:        int
                A numerical index to save the i-th structure

    axis:       string
                Set to "all", but it can be "x" "y" "z" in case of single axis rotation

    Returns
    ----------

    """
    #here we are genereting coordiantes file to visualize the situation
    tot_atoms = 3
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

    file_name = directory + "/structure_%s_%03d.mdcrd" %(axis,idx)
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
    #generate a pdb structure as well?

def histogram(nrgs,r,counter,maxlim=None):
    r"""Function to create a histogram with energies distribution for each
    ion-dipole distance

    Parameters
    ----------
    nrgs:       list
                list with all the energies at a distance r

    r:          float64
                distance ion-dipole

    Returns
    ----------

    """

    #the histogram is frequency vs energies
    directory = "histograms"
    if not os.path.exists(directory):
        os.makedirs(directory)
    name = "histograms/pot_%.4f.png" % (r)
    #50 bins = 10000 point for each bin
    hist,bins = np.histogram(nrgs,bins=10)#50)
    #design the maximum limit
    if counter==0:
        maxlim = 0.0
        for val in hist:
            if val > maxlim:
                maxlim = val
        maxlim = math.ceil(maxlim)
    else:
        pass

    width = 0.7*(bins[1] - bins[0])
    center = (bins[:-1] + bins[1:])/2
    fig,ax = plt.subplots()
    ax.bar(center,hist,align="center",width=width)
    ax.set_ylim(0,maxlim)
    x_axis_name = r"V(r), r=%.2f $\AA$ /kcal$\cdot$mol$^{-1}$" % r
    plt.xlabel(x_axis_name)
    plt.ylabel("Occurrence")
    plt.tight_layout()
    plt.savefig(name,transparent=True,dpi=300)
    plt.clf()
    plt.close("all")
    return maxlim

def bootstrap(values,idx):
    r"""Resampling energy list values for bootstraping

    Parameters
    ----------
    values:     list
                list of energies at a distance r

    idx:        int
                index of the process, necessary to use multiprocessing

    Returns
    ----------

    """
    nvals = len(values)
    new_values = []
    for x in range(0,nvals):
        i = random.randint(0,nvals-1)
        new_values.append(values[i])

    mean = np.mean(new_values)
    return mean

########################################MAIN#################################################################

top_file = sys.argv[1]  #topology
crd_file = sys.argv[2]  #coordiantes
#create the Sire System
amber = Amber()
molecules, space = amber.readCrdTop(crd_file, top_file)
system = createSystem(molecules)
#create a force field to have a reference for Coulombic interactions, it may be useful
system, energiesDict = setupForcefields(system, space,cutoff_type.val,cutoff_dist,rf_dielectric.val)
#Fix the dipole in space, so te middle point of the bond will be in the origin
#while the atoms will be in  b/2, 0, 0 and -b/2, 0, 0  where b is the bond lenth
#1) Fix the dipole onto the x axis (b/2,0,0) and (-b/2,0,0) --> dipole_set_x
#2) rotate the dipole. For each rotation at r comute nrg
print("Fixing dipole in space...")
dipole_set_x(system)
dipole_check(system)
r_vector = np.linspace(5,80,20) #vector of distances for the ion
outputf = open("r_coul.dat","w") #save distance, coulombic avg, coulombic std
dcdfile = open("dcd_creator.traj","w") #add an option for this file. It's useful to see
# the final structures - but not necessary if we see them correct
#now cycle for each ion distance we have
max_counter = 0
for r in r_vector:
	 print("Setting Sodium ion coordiantes: %.4f" % r)
	 sodium_coordinates(system,r) #fix the sodium coordinates
	 sodium_check(system) #sanity check
	 #now rotate teh dipole along all the possible axis, in order to have
	 #the best rotaitionally averaged approximation as possible
	 print("Rotation dipole for ion distance %.4f" % r)
	 #for each distance we have we are going to compute a histogram, to see the distribution
	 #of energies. What we want is a uniform distribution (so - and + correctly explore the
	 #same space)
	 if max_counter == 0 :
	     nrgs  = rotation_matrix(system,max_counter,dcdfile)
	     #maxlim = histogram(nrgs,r,max_counter)
	     dcdfile.write("trajout first.nc\n")
	     dcdfile.write("go\n")
	     dcdfile.close()
	     max_counter+=1
	 else:
	     nrgs = rotation_matrix(system,max_counter)
	     #histogram(nrgs,r,max_counter,maxlim)
	     max_counter+=1
	 #should we add the command: cpptraj  -p system.prmtop -i dcd_creatore.traj ?
	 #bootstraping routine
	 start = time.time()
	 print("Bootstrap ...")
	 nbootstrap = 100 #how many time we iterate
	 pool = Pool(processes=5) #create 5 processes to run in parallel, so  bootstrap will be
	 #really fast
	 idxs = np.linspace(0,nbootstrap,nbootstrap-1) #a series of indexes as an input for pool.map
	 partialboot = partial(bootstrap,nrgs) #partialise the bootstrap function
	 #in this way we can use pool.map and we are giving nrgs as a whole
	 means = pool.map(partialboot,idxs) #this will give as output means which is
	 #a list of all the mean computed by bootstrap function
	 pool.close() #close and join the processes
	 pool.join()

	 avg = np.mean(nrgs)
	 stdev = np.std(means)
	 end = time.time()
	 print(avg,stdev)
	 print("Bootstraping time %.4f" % (end-start))
	 #print("Bootstraping time: %.4f" % (end-start))
	 outputf.write("%.4f,%.18f,%.18f\n" % (r,avg,stdev))
