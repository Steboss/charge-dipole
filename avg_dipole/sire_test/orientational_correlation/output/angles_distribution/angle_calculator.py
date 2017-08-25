#July 2017 Stefano Bosisio
#Test to compute the orientation distribution function for ion and water
#Usage: python orientation_distribution_function.py  traj TOP WIDTH
#WIDTH is the width for the shell, it can be 1 A (5-6) or 0.5 A ( 5-5.5-6) or 0.2 (5-5.2-5.4...) and so on...

import sys,os
import math
import numpy as np
import time
import mdtraj as md
from parmed.amber import *
import multiprocessing
from multiprocessing import Pool, Process
from functools import partial
import matplotlib.pyplot as plt 

def compute_distance(na_coords,COM):
    r"""Function to compute the distance between sodium ion and coms

    Parameters
    ----------
    na_coords: list
    		list of coordinates of the sodium ion at a particular frame
	COM: 	list
    		list of coordinates of the water-COM within a shell
    ----------
    dist:   float64
    		distance betwene sodium and COM

    """
    dist = math.sqrt((na_coords[0] - COM[0])**2 + (na_coords[1]-COM[1])**2 +(na_coords[2]-COM[2])**2)
    return dist



def compute_angle(na_coords,COM,oxy_coords):
	r"""Function to compute the angle between ion and the COM, which is 
	the middle point of the dipole

	Parameters
	----------
	na_coords: list
				list of coordinates of the sodium ion at a particular frame

	COM:		list
				list of coordinates of water COM
	
	oxy_coords:	list
				list of coordinates of water oxygen

	Returns
	---------

	angle: 		float64
				angle ( in rad) between ion - COM - oxy, namely the angle between 
				the ion and the dipole

	"""

	#try cosine rule
	v1 = [na_coords[0] - COM[0], na_coords[1] - COM[1], na_coords[2] - COM[2] ]
	v2 = [oxy_coords[0] - COM[0], oxy_coords[1] - COM[1], oxy_coords[2] - COM[2]]
	#magnitude
	v1mag = math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
	v1norm = [ v1[0]/v1mag, v1[1]/v1mag, v1[2]/v1mag]

	v2mag = math.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
	v2norm = [ v2[0]/v2mag, v2[1]/v2mag, v2[2]/v2mag]

	res = v1norm[0]*v2norm[0] + v1norm[1]*v2norm[1] + v1norm[2]*v2norm[2]

	angle = math.acos(res)
	return angle #angle in rads

def readdcd(frame,framenumb,top,rcutoff,width=1.0,pdb=False,pdbcomoxy=False):
	r"""Main Function: this function reads each snapshot and computes the Coulombic energies:
	1) Barker Watts Coulombic interaction between the ion and the water charges-->BWRF-->BarkerWatts()
	2) Barker Watts Coulombic interaction between the ion and the dipole charges-->BWRF_charge_charge-->compute_ion_dipole
	3) Reaction Field Coulombic interaction between the ion and the dipole mu (as derived from Onsager)-->RF_charge_dipole-->compute_ion_dipole
	4) Rotationally averaged Reaction Field Coulombic interaction between ion and the dipole, using the 
	actual position of water molecules-->AVG_RF_charge_dipole-->compute_ion_dipole_rotationally_avg

	Parameters
	----------
	frame:	mdtraj_array
			array of coordinates of the n-th frame
	
	framenumb:	int
			frame number

	top:	mdtraj_topology
			topology information of the system

	pdb:	bool
			Whether or not save the whole system in pdb files
	
	rcutoff: float64	
			Radius of cutoff copmuted as the smallest box size
		
	width:  float64
			width of each shell, default 1 Angstrom
	
	pdbcomoxy:	bool
			Wheter or not save pdbs of COM and oxygens, to verify the ion-dipole angle

	Returns
	----------
	BWRF:	dictionary
			dictionary with all the energies for a specific shell
			e.g.  BWRF[3] = {1.00, 2.00,3.00...} where 3, in Angstrom, is the shell length
			and 1.00,2.00,3.00 the energies computed from all the ion-water interactions
	
	RF_charge_dipole: dictionary
			dictionary with all the energies for the RF charge-dipole interactions

	AVG_RF_charge_dipole: dictionary
			dictionary with all the energies for the RF AVG ROT charge-dipole interaction

	angles:	dictionary
			dictionary with all the ion-COM-oxygen angles, namely ion - dipole angle theta, for 
			each shell
	"""
	#energy dictionaries:
	#UNCOMMENT if you want to check the angles ion-COM-oxy = ion-dipole
	com_oxy_coords = {} # this saves the coords of com and oxygen for angle check
	#take the coordinates
	coordinates = frame[0]
	#number of atoms
	natoms = len(coordinates)
	#sodium coodinates
	na_coords = (coordinates[0][0],coordinates[0][1],coordinates[0][2])
	#UNCOMMENT below if you want to check angles
	#anglefile = open("angles.dat","w")
	angles_selection = {} # for comoxypdb only
	angles = {}
	counter = 0
	print("Extracting energies for each shell...")
	for res in top.residues:
		if "HOH" in res.name: #if the residue is water
			#sum over all the masses
			den_com = 0.0 
			#numerator = sum of masses*coords_ith
			num_comx = 0.0
			num_comy = 0.0
			num_comz = 0.0
			#
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
					oxy_coords = coordinates[atidx].tolist()
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
			#compute ion-COM distance
			dist = compute_distance(na_coords,COM)
			#define the shell where we are in 
			#first define the bounds
			if width == 1: 
				shell = math.ceil(dist) #take always the upper limit to enclose the distance within the hsell
				#print("Distance %.4f in shell %.4f" % (dist,shell))
			else:
				lower_shell = math.floor(dist)
				upper_shell = math.ceil(dist)
				#define all the shells
				all_shells = np.arange(lower_shell,upper_shell,width)#create as many shell as required by the width
				for limits in all_shells:
					#print("Comparing %.4f with %.4f" %(dist,limits))
					if dist > limits:
						shell = upper_shell #if the COM distance is greater than any other shell, the upper limits will be the shell
					else: 
						shell = limits #else the shell will be found and break the cycle
						break
				
			#BWRF charge*charge and charge*dipole
			angle_ion_com = compute_angle(na_coords,COM,oxy_coords)

			if shell in angles.keys(): 

				angles[shell].append(angle_ion_com)

			else:
				angles[shell] = [angle_ion_com]
			
			#UNCOMMNET below if you want to check angles
			#if shell < 10.0: #take only the  first shells, to simplify the visualization
			#	if shell in com_oxy_coords:
			#		com_oxy_coords[shell].append([COM,oxy_coords]) 
			#		angles[shell].append(angle_ion_com)
			#	else:
			#		com_oxy_coords[shell]=[]
			#		com_oxy_coords[shell].append([COM,oxy_coords])
			#		angles[shell] = [angle_ion_com]
			counter+=1

		else:
			continue
	
	#check for PDB creation
	if pdb:
		write_pdb(shell_coords,na_coords,framenumb)
	#if pdb_com_oxy
	if pdbcomoxy:
		#write the angle file
		keys = []
		for key in angles.keys():
			keys.append(key)
		keys.sort()
		for key in keys:
			for angle in angles[key]:
				#shell and angle
				anglefile.write("%.4f,%.4f\n" % (key,angle*(180/math.pi)))
		anglefile.close()
		write_pdbcomoxy(com_oxy_coords,na_coords,framenumb)

	return angles


def write_pdbcomoxy(coords,nacoords,idx):
	r"""Function to write pdb for ion-COM-oxygen system, to visualize the angle between
	them, which is the ion-dipole angle used in the calculation

	Parameters
	----------
	coords:		dictionary
				dictionary with the COM and oxygen coordinates:
				coords[shell] = [[COM],[OXY]] [[...]]
	nacoords:	list
				list of coordiantes of the sodium ion at frame idx-th
	idx:		int
				frame number

	Returns
	---------

	"""
	directory = "pdb_com_oxy"

	if not os.path.exists(directory):
		os.makedirs(directory)

	#sort the keys
	keys = []

	for key in coords.keys():
		keys.append(key)
	
	keys.sort()

	pdbfilename = "%s/structure_%d.pdb" %(directory,idx)
	pdbfile = open(pdbfilename,"w")



	atcount = 0
	chaincount = 0 

	for key in keys:
		for val in coords[key]:
			atcount+=1#atom number
			chaincount +=1 #necessary to divide different chain/oxygen-COM
			#
			#COM coords
			comx = val[0][0]
			comy = val[0][1]
			comz = val[0][2]
			#oxygen coords
			oxyx = val[1][0]
			oxyy = val[1][1]
			oxyz = val[1][2]

			#write COM coords
			final="{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format("ATOM",atcount,"H1"," ","LIG","A",chaincount," ",comx,comy,comz,1.00,1.00,"H","H")
			pdbfile.write(final)
			#write OXY
			atcount+=1 #update the atom counter
			final="{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format("ATOM",atcount,"O"," ","LIG","A",chaincount," ",oxyx,oxyy,oxyz,1.00,1.00,"O","0")
			pdbfile.write(final)


	final="{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format("ATOM",atcount+1,"Na"," ","NAA","A",chaincount+1," ",nacoords[0],nacoords[1],nacoords[2],1.00,1.00,"Na","0")
	pdbfile.write(final)
	pdbfile.write("END")
	pdbfile.close()



def write_pdb(coords,nacoords,idx):
	r"""Function to write pdb for all the system for the frame idx-th

	Parameters
	----------
	coords:		dictionary
				Coordinates of COMS and oxygens

	nacoords:	list
				list of coordiantes of the sodium ion at frame idx-th
	idx:		int
				frame number

	Returns
	---------

	"""
	directory = "pdb"

	if not os.path.exists(directory):
		os.makedirs(directory)

	#sort the keys
	keys = []

	for key in coords.keys():
		keys.append(key)

	keys.sort()

	pdbfilename = "%s/structure_%d.pdb" %(directory,idx)
	pdbfile = open(pdbfilename,"w")

	atcount = 0

	for key in keys:
		for val in coords[key]:
			atcount+=1
			x = val[0]
			y = val[1]
			z = val[2]

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

#Load traj, top and coord
traj = sys.argv[1]
top = sys.argv[2]
width = float(sys.argv[3]) # e.g. 1, 0.5, 0.2
print("Trajectory loading...")
print("This may take a while...")

mdtrajtop = md.load_prmtop(top) #this function allows the usage of SYSTEM.top
mdtrajdcd = md.open(traj,"r") #this function does not require huge memory allocation


try :
	#pass this info via a master script, to run this script at the same time
	#processing  different frames and taking less amount of time
	start = int(sys.argv[4]) # initial frame
	end = int(sys.argv[5]) #last frame to analyse
	#we can skip this part and analyse the whole trajectory with this script
except:
	start = 0 
	end = len(mdtrajdcd)
nframes = len(mdtrajdcd) # this is necessary to compute the radius of cutoff
#namely the smallest box size


angles_dict = {}


#comptue the average radius of cutoff, and skip the first 1000 frames
minimum = 100000 #set a minumim (higher than the expected box length)
print("Computing the minimum box length...")
for framenumber in range(0,nframes):#
	current, cell_lengths, angles = mdtrajdcd.read(n_frames=1) #for each frame extract the cell_lengths
	for length in cell_lengths[0]: #comparison
		if length< minimum:
			minimum = length

#set the cutoff
rcutoff =minimum
print("Computed average radius of cutoff %.4f" % rcutoff)
#reverse the trajectory and let's the main function commence
mdtrajdcd.seek(0)
#run through the shell calculation
for framenumber in range(start,end):
	print("Processing frame %d..." % framenumber)
	current, cell_lengths, angles = mdtrajdcd.read(n_frames=1)
	#main function:
	angles_shell= readdcd(current,framenumber,mdtrajtop,rcutoff,width)

	
	#add all the  results to each dictionary
	for key in angles_shell.keys():
		if key in angles_dict: 
			for single_angle in angles_shell[key]:
				angles_dict[key].append(single_angle)

		else:
			counter = 0 
			for single_angle in angles_shell[key]:
				if counter ==0:
					angles_dict[key] = [single_angle]
					counter+=1
				else:
					angles_dict[key].append(single_angle)
			counter = 0 

#once we have stored all the average energies for each shell, compute the 
#total average
#create an output directory

#first plot single histograms
max_val = 0 
outputfold = "output_single_%s_%s_%.2f" % (start,end,width) 
if  not os.path.exists(outputfold):
    os.makedirs(outputfold)

for key in angles_dict.keys():
	filename = ("%s/hist_%.4f.png" % (outputfold,key))

	hist,bins = np.histogram(angles_dict[key],bins=20)
	histwidth = 0.7*(bins[1] - bins[0])
	center = (bins[:-1] + bins[1:])/2
	fig,ax = plt.subplots()
	ax.bar(center,hist,align="center",width=histwidth)
	plt.xlabel(r"$\theta$ / rad", fontsize=20)
	plt.ylabel("Occurrence",fontsize=20)
	plt.savefig(filename)
	plt.clf()
	plt.close("all")
	tmp_max = max(hist)
	if tmp_max > max_val:
		max_val = tmp_max
