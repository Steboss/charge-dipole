#Aug2017 Stefano Bosisio
#Script to parse all results from the output_X_X directories

import sys,os
import math
import numpy as np
import time
import mdtraj as md
from parmed.amber import *
import multiprocessing
from multiprocessing import Pool, Process
from functools import partial
import glob

######
#MAIN#
######

inptf = glob.glob("output_*/*.csv")

BWRF= {}
BWRF_charge_charge = {}
RF_charge_dipole = {}
AVG_RF_dipole = {}
AVG_RF_int_dipole = {}
shells = [] # create a list to store the shells values, 
#we will use it to go through later for printing ordered shells
print("Reading all the output files...")
for files in inptf:
    print("Processing %s" %files)
    #first find the kind  of energy:
    ifile = open(files,"r").readlines()
    #recognize the kind of Coulombic potential and store it in the
    #appropriate dictionary
    energy_type = files.split("/")[1].split(".csv")[0]  
    for lines in ifile:
        shell = float(lines.split(",")[0])
        nwaters = float(lines.split(",")[1])
        nwatersstd = float(lines.split(",")[2])
        nrg = float(lines.split(",")[3])
        nrgstd = float(lines.split(",")[4])
        
        if shell in shells:
            pass
        else:
            shells.append(shell)
        if "BWRF_charge_charge" in energy_type:

            if shell in BWRF_charge_charge:
                BWRF_charge_charge[shell].append(nrg)
            else:
                BWRF_charge_charge[shell]=[]
                BWRF_charge_charge[shell].append(nrg)

        elif "AVG_RF_int_dipole" in energy_type:
            if shell in AVG_RF_int_dipole:
                AVG_RF_int_dipole[shell].append(nrg)
            else:
                AVG_RF_int_dipole[shell] = []
                AVG_RF_int_dipole[shell].append(nrg)

        elif "AVG_RF_dipole" in energy_type:
            if shell in AVG_RF_dipole:
                AVG_RF_dipole[shell].append(nrg)
            else:
                AVG_RF_dipole[shell] = []
                AVG_RF_dipole[shell].append(nrg)

        elif "RF_charge_dipole" in energy_type:
            if shell in RF_charge_dipole:
                RF_charge_dipole[shell].append(nrg)
            else:
                RF_charge_dipole[shell]=[]
                RF_charge_dipole[shell].append(nrg)
        else:
            if shell in BWRF:
                BWRF[shell].append(nrg)
            else:
                BWRF[shell]=[]
                BWRF[shell].append(nrg)
print("Created dictionary")
shells.sort()
#since the shells are common to all te dictionaries, we can cycle through 
#BWRF.keys() and compute the average and std err
outputfolder = "output_average"
if  not os.path.exists(outputfolder):
    os.makedirs(outputfolder)
print("Writing files...")
BWRF_file = open("%s/BWRF.csv"%outputfolder,"w")#BWRF results
BWRF_charge_charge_file = open("%s/BWRF_charge_charge.csv"%outputfolder,"w")#BWRF_charge_charge
RF_charge_dipole_file = open("%s/RF_charge_dipole.csv"%outputfolder,"w")#RF_charge_dipole
AVG_RF_dipole_file = open("%s/AVG_RF_dipole.csv"%outputfolder,"w") # AVG_RF dipole
AVG_RF_int_dipole_file = open("%s/AVG_RF_int_dipole.csv" % outputfolder,"w")#AVG_RF_dipole integrated

for key in shells:
    #output: shell, nrg, stddev
    BWRF_file.write("%.4f,%.8f,%.8f\n" % (key,np.mean(BWRF[key]), np.std(BWRF[key])))
    BWRF_charge_charge_file.write("%.4f,%.8f,%.8f\n" %(key,np.mean(BWRF_charge_charge[key]),np.std(BWRF_charge_charge[key])))
    RF_charge_dipole_file.write("%.4f,%.8f,%.8f\n" %(key,np.mean(RF_charge_dipole[key]),np.std(RF_charge_dipole[key])))
    AVG_RF_dipole_file.write("%.4f,%.8f,%.8f\n"%(key,np.mean(AVG_RF_dipole[key]),np.std(AVG_RF_dipole[key])))
    AVG_RF_int_dipole_file.write("%.4f,%.8f,%.8f\n" %(key,np.mean(AVG_RF_int_dipole[key]),np.std(AVG_RF_int_dipole[key])))

BWRF_file.close()
BWRF_charge_charge_file.close()
RF_charge_dipole_file.close()
AVG_RF_dipole_file.close()
AVG_RF_int_dipole_file.close()
print("Everything's done :)")