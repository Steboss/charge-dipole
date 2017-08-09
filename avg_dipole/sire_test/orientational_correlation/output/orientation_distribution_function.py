	#July 2017 Stefano Bosisio
#Test to compute the orientation distribution function for ion and water

import sys,os
import math
import numpy as np
import time
import mdtraj as md
from parmed.amber import *
import multiprocessing
from multiprocessing import Pool, Process
from functools import partial

one_over_four_pi_eps0 = 332.0637090025476
dielectric = 82.0
mu = 0.48 #eA
qion = 1.0
epsilon1 = (dielectric -1)/(1+2*dielectric)
epsilon2 = (3/(1+2*dielectric))
constqmu = -one_over_four_pi_eps0*(qion*mu)**2


def BarkerWatts(na_coords,res_coords,r_cutoff):
    r"""Function to compute BW reaction field Coulombic interaction potential

    Parameters
    ----------
    na_coords:   list
                 list of  coordinates (x,y,z) of sodium ion
    res_coords:  list
                 list of coordinates for the three atoms in the water molecule
    r_cutoff:    float64
                 shell length
    Returns
    ----------
    pot:   		float64
        		Coulombic potential between ion and h1 and h2 and o
    """
	#r_cutoff = shell length
    nacharge = 1.0

    krf = (1/r_cutoff**3)*((dielectric-1)/(2*dielectric+1))
    pot = 0.0

    for key in res_coords.keys():
    	if "H" in key:
            atcharge = 0.417
     	elif "O" in key:
    	    atcharge = -0.834
     	else:
     	    continue
        #compute distance na and res atoms
        dist=compute_distance(na_coords,res_coords[key])
        pot += one_over_four_pi_eps0*(nacharge*atcharge)*((1/dist) + krf*dist**2)

    return pot

def theory(natoms,cutoff,na_coords,shell_coords):
    r"""Function to compute RF rot-avg charge dipole interaction

    Parameters
    ----------
    natoms:	int
    		number of atoms in a shell
    cutoff: float64
    		radius of cutoff, which is the shell length
    na_coords: list
    		list of coordinates of the sodium ion at a particular frame
    shell_coords: list
    		list of coordinates of the water-COM within a shell
    ----------
    pot_avg:   float64
    		average RF rot-avg potential for that cell for natoms

    """
    pot = 0.0
    counter = 0

    for coords in shell_coords:
    	dist = compute_distance(na_coords,coords)
    	term1 = (4*(epsilon1 + epsilon2)**2)*((-1/(3*dist**3)) + (1/(3*cutoff**3)) )
        term2 = (4*(epsilon1**2)/cutoff**3) *((dist**4)/4 - (cutoff**4)/4)
        term3 = -2*epsilon1*(math.log(dist) - math.log(cutoff))
        term4 = -2*epsilon2*(math.log(dist) - math.log(cutoff))
        pot += (constqmu*(term1 + term2 + term3 + term4))
        counter +=1

    pot_avg = pot/counter

    return pot_avg


def compute_distance(na_coords,com):
    r"""Function to compute the distance between sodium ion and coms

    Parameters
    ----------
    na_coords: list
    		list of coordinates of the sodium ion at a particular frame
    com: 	list
    		list of coordinates of the water-COM within a shell
    ----------
    dist:   float64
    		distance betwene sodium and COM

    """
    dist = math.sqrt((na_coords[0] - com[0])**2 + (na_coords[1]-com[1])**2 +(na_coords[2]-com[2])**2)
    return dist



def readdcd(frame,top):
    r"""Main function to read each frame of the trajectory

    Parameters
    ----------
    frame:	mdtraj_array
    		array of coordinates from frame nth
    top:	mdtraj_topology
    		topology information of the system
    ----------

    nrg_dict: dictionary
    		dictionary with list of energies for that shell, for example:
    		nrg_dict[6] = [1,2,3,45,6] where 6 is in Angstrom and energies are
    		stored in a list
    shell_coords: dictionary
    		dictionary with list of coordinate for all the water COM  in a shell
    		for example:
    	    shell_coords[6] = [X,Y,Z], [X,Y,Z] ...
    	    where 6 is in Anstrom, is the shell length
    	    X,Y, Z the COM coordinates
    na_coords: list
    		list of coordinates for the sodium ion at a particular frame

    """
    nrg_dict =  {}
    shell_coords = {}
    coordinates = frame[0]
    natoms = len(coordinates)

    na_coords = (coordinates[0][0],coordinates[0][1],coordinates[0][2])

    counter = 0
    print("Extracting energies for each shell...")
    for res in top.residues:
		if "HOH" in res.name:
			den_com = 0.0 #sum ovr all the masses
			num_comx = 0.0 #sum of masses*coords
			num_comy = 0.0
			num_comz = 0.0
			res_coords = {}
			for at in res.atoms:
				#compute the com
				atname = at.name
				atidx = at.index
				res_coords[atname] = coordinates[atidx]
				if "H" in atname:
					atmass = 1.008 #masses and charges are not implemented in mdtraj
					atcharge = 0.417
				elif "O" in atname:
					atmass = 16.0
					atcharge = -0.834
				else:
					continue #skip if it's an ion - shouldn't

				den_com +=atmass
				num_comx +=(atmass*coordinates[atidx][0])
				num_comy +=(atmass*coordinates[atidx][1])
				num_comz +=(atmass*coordinates[atidx][2])

			#com calculation
			X_com = num_comx/den_com
			Y_com = num_comy/den_com
			Z_com = num_comz/den_com
			COM = [X_com,Y_com,Z_com]
			#now we have two alternatives
			#1) Compute now the BW for any shell we are in
			#2) Store in a dictionary the shell number and the res index
			#then compute back in another function for at in residues[index]
			#define in what shell we are by taking the ceil
			dist = compute_distance(na_coords,COM)
			#define the shell where we are
			shell = math.ceil(dist)
			#SANITY CHECK: save the COM coords in each "shell" and then create
			#a corodiantes file to visualize it
			if shell in shell_coords:
				shell_coords[shell].append(COM)
			else:
				shell_coords[shell] = [COM]
				#now compute the energy BW
			nrg = BarkerWatts(na_coords,res_coords,shell)
			#for each shell store the energy values in a dictionary
			if shell in nrg_dict:
				nrg_dict[shell].append(nrg)
			else:
				nrg_dict[shell] = [nrg]
			
			counter+=1

		else:
			continue

    return nrg_dict,shell_coords,na_coords


def write_pdb(coords,nacoords,idx):

    directory = "pdb"

    if not os.path.exists(directory):
        os.makedirs(directory)

    pdbfilename = "pdb/structure_%d.pdb" %idx
    pdbfile = open(pdbfilename,"w")

    keys = []

    for key in coords.keys():
    	keys.append(key)

    keys.sort()
    atcount = 0

    for key in keys:
    	for val in coords[key]:
    	    atcount+=1
    	    x = val[0]
    	    y = val[1]
    	    z = val[2]

    	    #outline="%8.3f%8.3f%8.3f" %(x,y,z)

    	    #final ="ATOM      %d  H1  LIG     1%s  1.00  0.00           O\n" % (atcount, outline)
    	    final="{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format("ATOM",atcount,"H1"," ","LIG","A",1," ",x,y,z,1.00,1.00,"H","0")
            pdbfile.write(final)

    final="{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format("ATOM",atcount+1,"Na"," ","NAA","A",1," ",nacoords[0],nacoords[1],nacoords[2],1.00,1.00,"Na","0")
    pdbfile.write(final)
    pdbfile.write("END")
    pdbfile.close()
    #sys.exit(-1)



############
####MAIN####
############

#devid the trajectory in shells
#for each shell compute the Barker Watts coulombic interactions  with each water molecule
#compute with the theoretical derivation N_water in the sshell * the integral
#what if a water is between two shells?
#for i in shells :
#at fshells 1 find out the watermolecules within that wshell for each frame
#for each frame compute the distance betwen ion and the centre of mass of th emolecule
traj = sys.argv[1]
top = sys.argv[2]
crd = sys.argv[3]

print("Trajectory loading...")
print("This may take a while...")

mdtrajtop = md.load_prmtop(top)
mdtrajdcd = md.open(traj,"r")
nframes = len(mdtrajdcd)

mean_dict = {}
theory_dict = {}
counter=0
#now cycle throught each frame
for framenumber in range(1,nframes):
	print("Processing frame %d..." % framenumber)
	current, cell_lengths, angles = mdtrajdcd.read(n_frames=1)
	nrg_dict,coord_dict,na_coords= readdcd(current,mdtrajtop)
	#take the average energy for each shell according to the theory
	#for shell in nrg_dict.keys():
	#	th_pot_avg = theory(len(nrg_dict[shell]),shell, na_coords,coord_dict[shell])
	#	if shell in theory_dict:
	#		theory_dict[shell].append(th_pot_avg)
	#	else:
	#		theory_dict[shell] = [th_pot_avg]

	write_pdb(coord_dict,na_coords,counter)
	counter+=1
	#if counter ==100:
	#	break
	#now compute the average for each shell
	for key in nrg_dict.keys():

		if key in mean_dict:
			avg = np.mean(nrg_dict[key])
			#std = np.std(nrg_dict[key])
			mean_dict[key].append(avg)
		else:
			avg = np.mean(nrg_dict[key])
			#std = np.std(nrg_dict[key])
			mean_dict[key] = [avg]

#now compute the average energy for each shell

for key in mean_dict.keys():
	mean_dict[key] = [np.mean(mean_dict[key]),np.std(mean_dict[key])]
#compute the theoretical mean and std
#for key in theory_dict.keys():
#	theory_dict[key] = [np.mean(theory_dict[key]),np.std(theory_dict[key])]
for key in mean_dict.keys():
	print("COMPUTED Shell %d  nrg %.4f +/- %.4f" % (key, mean_dict[key][0],mean_dict[key][1]))
	#print("THEORY   Shell %d  nrg %.4f +/- %.4f" % (key,theory_dict[key][0],theory_dict[key][1]))
