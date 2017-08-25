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

one_over_four_pi_eps0 = 332.0637090025476
dielectric = 82.0

def Coulombic(na_coords,res_coords,rcutoff,distionCOM):
    r"""Function to compute BW reaction field Coulombic interaction potential

    Parameters
    ----------
    na_coords:   list
                 list of  coordinates (x,y,z) of sodium ion
    res_coords:  list
                 list of coordinates for the three atoms in the water molecule
    rcutoff:     float64
				 average minimal box length, through all the snapshots
                 
    Returns
    ----------
    pot:   		float64
        		Coulombic potential between ion and water's atoms
    """

    nacharge = 1.0  #give the ion charge in e
    pot = 0.0 #initialise potential 

    for key in res_coords.keys():
		#hard code the charges, since mdtraj cannot retrieve these info

		if "H" in key:
			atcharge = 0.417
		elif "O" in key:
			atcharge = -0.834
		else:
			continue

		
        #compute distance na and res atoms
		dist=compute_distance(na_coords,res_coords[key]) 
		pot += one_over_four_pi_eps0*(nacharge*atcharge)*((1/dist))


    return pot


def compute_ion_dipole(na_coords,COM,rcutoff,angle):
	r"""Function to compute ion-water potential with water as dipole:
	1) The potential is computed as ion*charge, with charge from dipole 
	charges, +/-0.834 
	2) The potential is computed as ion*dipole, from the dervied reaction field
	of Onsager

    Parameters
    ----------
    na_coords:   list
                 list of  coordinates (x,y,z) of sodium ion

    COM:  		list
    			list of coordinates of the water-COM within a shell

    rcutoff:     float64
				 average minimal box length, through all the snapshots
	angle:		float64
				angle in rad  betwene  the ion and the COM of the water
                 
    Returns
    ----------
    pot_charge_charge:   		float64
        		Coulombic potential computed as BWRF * dipole charges
	pot_charge_dipole:	float64
				Coulombic potential computed as RF for ion-dipole
	"""

	#in this case we have the water as a dipole
	#
	#Euclidean distance
	r = compute_distance(na_coords,COM)
	bondlength = 0.586 #A
	h1 = -0.834 #e
	h2 = +0.834 #e
	qion = +1 #e
	mu = h2*bondlength # dipole moment in eA
	#the dipole bond length is 0.586 A
	#dist1 = r - 0.5*bondlength*math.cos(angle) # for the negative charge -0.834 e
	#dist2 = r + 0.5*bondlength*math.cos(angle) # positive charge  +0.834 e
	#compute the energy in two ways:
	#1) charge * charge BWRF
	#2) charge * dipole BWRF
	#krf= (1/rcutoff**3)*((dielectric-1)/(2*dielectric + 1))
	#pot_ion1 = one_over_four_pi_eps0*(qion*h1)*( (1/dist1) + krf*dist1**2)
	#pot_ion2 = one_over_four_pi_eps0*(qion*h2)*( (1/dist2) + krf*dist2**2)
	#pot_charge_charge = pot_ion1 + pot_ion2

	#2)charge*dipole, derived from Onsager
	pot_charge_dipole = -one_over_four_pi_eps0*qion*mu*(math.cos(angle))*(1/r**2)  
	
	return pot_charge_dipole #they should be the same

def compute_ion_dipole_rotationally_avg(na_coords,COM,rcutoff):
	r"""Function to compute the rotationally averaged charge-dipole reaction field
	interaction

    Parameters
    ----------
    na_coords:   list
                 list of  coordinates (x,y,z) of sodium ion
    COM:  		list
    			list of coordinates of the water-COM within a shell
    rcutoff:     float64
				 average minimal box length, through all the snapshots
                 
    Returns
    ----------
    pot:   		float64
        		Coulombic potential computed as  RotAvgRF
	"""
	r = compute_distance(na_coords,COM)
	qion = 1.0
	bondlength = 0.586 # A
	mu = 0.834*bondlength
	inv_KT = -1/(3*0.592)

	pot = inv_KT*((qion*mu)**2)*(one_over_four_pi_eps0**2)*(1/r)**4

	return pot

def compute_ion_dipole_rotationally_avg_integrated(nwaters,position,rcutoff):
	r"""Function to compute the rotationally averaged charge-dipole reaction field as
	an integral, thus = Nwaters*potential
	interaction

    Parameters
    ----------
    nwaters:   int
               number of water moelcule in the shell
	position: float64
				position of the water molecule, given by the averaged of the shell bounds
	rcutoff:	float64
				radius of cutoff, minimum box length
                 
    Returns
    ----------
    pot:   		float64
        		Coulombic potential computed as  RotAvgRF

	"""
	r = position
	ratio = (r**3)/(rcutoff**3)
	qion = 1.0
	bondlength = 0.586 # A
	mu = 0.834*bondlength#eA
	inv_KT = -1/(3*0.592)

	npot = nwaters*inv_KT*((qion*mu)**2)*(one_over_four_pi_eps0**2)*(1/r)**4
	return npot
	


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
	Coulombic:	dictionary
			dictionary with all the energies for a specific shell with Coulombic potential
			e.g.  Coulombic[3] = {1.00, 2.00,3.00...} where 3, in Angstrom, is the shell length
			and 1.00,2.00,3.00 the energies computed from all the ion-water interactions
	
			print(dist)
			sys.exit(-1)
	charge_dipole: dictionary
			dictionary with all the energies for the charge-dipole interactions

	avg_rot_charge_dipole: dictionary
			dictionary with all the energies for the rotationally averaged charge-dipole interaction

	angles:	dictionary
			dictionary with all the ion-COM-oxygen angles, namely ion - dipole angle theta, for 
			each shell
	"""
	#energy dictionaries:
	coulomb =  {} #-->Coulombic 
	#BWRF_charge_charge = {} #-->compute_ion_dipole
	charge_dipole = {} #-->compute_ion_dipole
	avg_rot_charge_dipole = {} #-->compute_ion_dipole_rotationally_avg
	shell_coords = {} # this is necessary for checking all the system via pdbs
	angles = {} # dictionary to store ion-dipole angle values
	#UNCOMMENT if you want to check the angles ion-COM-oxy = ion-dipole
	#com_oxy_coords = {} # this saves the coords of com and oxygen for angle check
	#take the coordinates
	coordinates = frame[0]
	#number of atoms
	natoms = len(coordinates)
	#sodium coodinates
	na_coords = (coordinates[0][0],coordinates[0][1],coordinates[0][2])
	#UNCOMMENT below if you want to check angles
	#anglefile = open("angles.dat","w")
	#angles = {}
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
				
				
			#BWRF energy
			nrg = Coulombic(na_coords,res_coords,rcutoff,dist)
			#BWRF charge*charge and charge*dipole
			angle_ion_com = compute_angle(na_coords,COM,oxy_coords)
			ion_dipole_nrg = compute_ion_dipole(na_coords,COM,rcutoff,angle_ion_com)
			#Rot-avg RF. In this case compute the potential as a sum of contributions
			#between ion and water-dipole using as a distance the COM
			ion_dipole_rot_nrg = compute_ion_dipole_rotationally_avg(na_coords,COM,rcutoff)
			#now store every values:
			if shell in shell_coords: 
				#shellcoords is useful to generate pdb of the whole system
				#here if the shell is in shell_coords it means that
				#it's present also in all the other dictionary
				shell_coords[shell].append(COM) #pdb purposes
				coulomb[shell]+=nrg #BWRF
				#BWRF_charge_charge[shell]+=ion_charge_nrg # BW ion charge charge
				charge_dipole[shell]+=ion_dipole_nrg # RF ion dipole interaction
				avg_rot_charge_dipole[shell]+=ion_dipole_rot_nrg#Rot AVG RF ion dipole interaction, with
				angles[shell].append(angle_ion_com)
				#position (r) given by the COM

			else:
				shell_coords[shell] = [COM]
				coulomb[shell] = nrg
				#BWRF_charge_charge[shell]=ion_charge_nrg
				charge_dipole[shell]=ion_dipole_nrg
				avg_rot_charge_dipole[shell]=ion_dipole_rot_nrg
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
		try:
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
		except:
			print("Trying to write the pdb for ion-COM-oxygen system but some errors occur")
			print("\n Did you try to uncomment the relevant part for the angle script to work?")


	return coulomb,charge_dipole,avg_rot_charge_dipole,shell_coords,angles


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


#intialize the dictionary
coulomb = {} # BWRF charge-charge-charge --> BarkerWatts()
#BWRF charge vs dipole's charges, in this case the distance ion-charges are computed
#as a function of the angle between the ion and the middle point of the dipoleg
#the whole dipole moment mu and multiplying it with cos(theta)
charge_dipole = {} #-->compute_ion_dipole
#rotationally averaged charge dipole reaction field interaction
#this is compute by placing as position the actual position of the waters'COM
avg_rot_charge_dipole = {} #-->compute_ion_dipol_rotationally_avg
#rotationally averaged charge dipole reaction field interaction
#computed as N*potential where N is the number of water molecule per shell
#and waters'positions are considered at half of the shell length
avg_rot_charge_dipole_int = {} #-->compute_ion_dipole_rotataionlly_avg_integrated
#create a dictionary to store the average number of water molecule in each shell, this 
#could be useful
average_water_number = {}
#store the values of angles-dipole for each shell
angles_dict = {}

#define a counter
counter=0
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
	nrg_coulomb,nrg_charge_dipole,nrg_avg_charge_dipole,shells,angles_shell= readdcd(current,framenumber,mdtrajtop,rcutoff,width)

	#compute the integrated rotationally averaged charge-dipole RF interaction
	for key in shells.keys(): #take the nrg_BWRF keys, which are the same for all the dictionaries
		position = key - (width/2.0) #this is the average position of the molecule - 1 A width shell
		nwaters = len(shells[key]) #this is the number of water molecule within the shell
		nrg_avg_integral = compute_ion_dipole_rotationally_avg_integrated(nwaters,position,rcutoff)

		if key in avg_rot_charge_dipole_int:
			avg_rot_charge_dipole_int[key].append(nrg_avg_integral)
			average_water_number[key].append(len(shells[key]))
		else:
			avg_rot_charge_dipole_int[key] = [nrg_avg_integral]#REMEMBER TO AVERAGE IT!
			average_water_number[key] = [len(shells[key])]

	#add all the  results to each dictionary
	for key in shells.keys():
		if key in coulomb: #this means that the key is in all the dictionaries
			coulomb[key].append(nrg_coulomb[key])
			charge_dipole[key].append(nrg_charge_dipole[key])
			avg_rot_charge_dipole[key].append(nrg_avg_charge_dipole[key])
			#compute the average angle for each shell
			angles_dict[key].append(np.mean(angles_shell[key]))

		else:
			coulomb[key]=[nrg_coulomb[key]]
			charge_dipole[key]=[nrg_charge_dipole[key]]
			avg_rot_charge_dipole[key]=[nrg_avg_charge_dipole[key]]
			angles_dict[key] = [np.mean(angles_shell[key])]

	counter+=1
	#end = time.time()
	#print("Time %.4f" %(end-start))
	#if counter ==100:
	#	break

#once we have stored all the average energies for each shell, compute the 
#total average
#create an output directory


outputfold = "output_%s_%s_%.2f" % (start,end,width) 
if  not os.path.exists(outputfold):
    os.makedirs(outputfold)

#compute the average number of water mols for each shell
water_file = open("%s/n_waters.csv" % outputfold,"w")
for key in average_water_number.keys():
	average_water_number[key]=[np.mean(average_water_number[key]),np.std(average_water_number[key])]
	#write the average water molecules in a file
	#mean,std
	water_file.write("%.4f,%.4f,%.4f\n" % (key,average_water_number[key][0],average_water_number[key][1]))
#write the file with all the angles
angle_file = open("%s/angles.csv" % outputfold,"w")
for key in angles_dict.keys():
	angles_dict[key]= [np.mean(angles_dict[key]),np.std(angles_dict[key])]
	angle_file.write("%.4f,%.4f,%.4f\n" % (key,angles_dict[key][0],angles_dict[key][1]))

#outputfile should containt
#shell,,<energy>,std.dev(energy)
#where shell is the shell length in A, <water mols> the average number of water molecules within that shell
#std.dev(water mols) the standard deviation of the nwater mols numb,er, <energy> the average energy of interaction
#between ion and water in the shell, std.dev(energy) the std.dev of the energy

coulomb_file = open("%s/Coulomb.csv"%outputfold,"w")
for key in coulomb.keys():

	coulomb[key] = [np.mean(coulomb[key]),np.std(coulomb[key])]

	coulomb_file.write("%.4f,%.8f,%.8f\n" % (key,coulomb[key][0],coulomb[key][1]) )

charge_dipole_file = open("%s/charge_dipole.csv"%outputfold,"w")
for key in charge_dipole.keys():
	charge_dipole[key] = [np.mean(charge_dipole[key]),np.std(charge_dipole[key])]
	charge_dipole_file.write("%.4f,%.8f,%.8f\n" % (key,charge_dipole[key][0],charge_dipole[key][1]))

avg_rot_charge_dipole_file = open("%s/avg_rot_charge_dipole.csv"%outputfold,"w")
for key in avg_rot_charge_dipole.keys():
	avg_rot_charge_dipole[key] = [np.mean(avg_rot_charge_dipole[key]),np.std(avg_rot_charge_dipole[key])]
	avg_rot_charge_dipole_file.write("%.4f,%.8f,%.8f\n" % (key,avg_rot_charge_dipole[key][0],avg_rot_charge_dipole[key][1]))


avg_rot_charge_dipole_int_file = open("%s/avg_rot_charge_dipole_int.csv"%outputfold,"w")
for key in avg_rot_charge_dipole_int.keys():
	avg_rot_charge_dipole_int[key] = [np.mean(avg_rot_charge_dipole_int[key]),np.std(avg_rot_charge_dipole_int[key])]
	avg_rot_charge_dipole_int_file.write("%.4f,%.8f,%.8f\n" % (key,avg_rot_charge_dipole_int[key][0],avg_rot_charge_dipole_int[key][1]))

#add a readme to each folder
readme_file = open("%s/README" % outputfold,"w")
readme_file.write("""Coulomb.csv:  Coulombic interaction 
charge_dipole.csv : Charge-Dipole interaction (dipole-ion)
avg_rot_charge_dipole.csv: Rotationally Averaged charge-dipole interaction, setting position r with COM positions
avg_rot_charge_dipole_int.csv: Rotationally Averaged charge-dipole  interaction, computed as N_waters*potential, setting molecules in the middle of each shell 
angles.csv: average angle and std.dev for each shell, in rads
n_waters.csv: average numbe rof water molecule for each shell and std.dev""")
water_file.close()
angle_file.close()
coulomb_file.close()
charge_dipole_file.close()
avg_rot_charge_dipole_file.close()
avg_rot_charge_dipole_int_file.close()
readme_file.close()