#AUGUST 2017 Stefano Bosisio
#Scrip to create a plot for each bondlength dipole to see the goodness of
#converged distance

import matplotlib.pyplot as plt
import matplotlib as matplotlib
import math,re,os,sys
import numpy as np
import seaborn as sbn 
import glob
sbn.set_style("whitegrid")

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
		deg = ifolder.split("/")[0].split("_")[0]

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



#MAIN
computed = extract_potential(glob.glob("*_deg/computed.dat"))
theory = extract_potential(glob.glob("*_deg/theory.dat"))

#create the plot for 0-180 every 30
angles = ["0","30","60","90","120","150","180"]


color = sbn.color_palette()

for key in angles:
	fig,ax =  plt.subplots(figsize=(10,5))
	ax.plot(computed[key][0],computed[key][1],color=color[0],marker="o",label="Computed %s" % key)
	ax.plot(theory[key][0],theory[key][1],color=color[1],marker="x",label="Theory %s" % key)

	ax.set_xlabel(r"r / $\AA$",fontsize=20)
	ax.set_ylabel(r"V(r) / kcal mol$^{-1}$",fontsize=20)
	ax.tick_params(labelsize=18)
	ax.legend(loc="best",fontsize=20)
	figname="Pot-%s.png" % key
	plt.savefig(figname,dpi=300)
	plt.clf()
