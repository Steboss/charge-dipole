#AUG 2017 Stefano Bosisio
#Script to extract the average distacne between all teh ssimulated cases
#USAGE:
#python average_distance.py "l_0.117/*_deg/computed.dat" "l_0.117/*_deg/theory.dat" 0.117
#as input: computed files, theory files, dipole bond length
import sys,os
import math
import numpy as np
import re
import glob
import glob

def extract_potential(folders):
	"""Function to extract the potential energy from each
	file.
	For example, in folder l_0.117 we are simulating a dipole
	length of 0.117 A. In this folder there are *_deg folders
	with results for the theory (theory.dat) and computed (computed.dat)
	ion-dipole interactions"""

	results = {}

	for ifolder in folders:
		distance = []
		pots = []
		#we are giving l_*/*_deg/computed
		deg = ifolder.split("/")[1].split("_")[0]

		reader = open(ifolder,"r").readlines()
		for read in reader:
			r = float(read.split(",")[0])
			pot = float(read.split(",")[1])
			distance.append(r)
			pots.append(pot)

		results[deg] = [distance, pots]

		#so results[deg][0] will be dist
		#results[deg][1] potential

	return results

def natural_sort(list_to_sort):
	"""Function to sort a list in a natural order"""
	convert = lambda text: int(text) if text.isdigit() else text.lower()
	alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)',key)]

	return sorted(list_to_sort,key=alphanum_key)


#MAIN

computedfiles = sys.argv[1]
theoryfiles = sys.argv[2]
bondlength = sys.argv[3]

print("Extracting  potential from input files...")
computed = extract_potential(glob.glob(computedfiles))
theory = extract_potential(glob.glob(theoryfiles))

#now compute when the computed and theory differences
#are less than <0.001 kcal/mol

match = {}
print("Matching theory with computational results...")


counter = 0 
for key in computed.keys():
	#cycle through the potentials
	for i,pots in enumerate(computed[key][1],0):
		diff = pots - theory[key][1][i]
		if abs(diff)<0.01:
			if counter == 0:
				dist = computed[key][0][i]
				tmp = dist
				counter+=1

			elif counter==1:
				#compare up to three elments
				dist = computed[key][0][i]
				if (dist-tmp)>3.0 :
			            #this is an error and we found the correct distance
			            match[key] = dist
			            counter = 0 
				    break
			        else:
			    	    tmp = dist
			    	    counter+=1
			elif counter==2:
				#compare up to three elments
				dist = computed[key][0][i]
				if (dist-tmp)>3.0 :
			            #this is an error and we found the correct distance
			            match[key] = dist
			            counter = 0 
				    break
			        else:
			    	    tmp = dist
			    	    counter+=1

			elif counter==3:
				#now last opportunity to find a wrong element
				dist = computed[key][0][i]
				if (dist-tmp)>3.0 :
				    match[key]=dist
				    counter=0
				    break
				else:
				    #save the previous one
				    match[key]=tmp
				    counter = 0 
				    break

                        else:
                            continue
#sort the keys of match

keys = []
for key in match.keys():
	keys.append(key)

sort_keys  = natural_sort(keys)

#save everything onto one file

ofile = open("theta_dist.csv","a")
print("Appending final results to theta_dist.csv file...")
for keys in sort_keys:
	ofile.write("%.4f,%.4f,%.4f\n" % (float(bondlength),float(keys),float(match[keys])))
print("All right. Goodbye :) ")
