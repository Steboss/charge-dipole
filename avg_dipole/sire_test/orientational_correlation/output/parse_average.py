#Aug2017 Stefano Bosisio
#Script to parse all results from the output_X_X directories
#Usage:
#python parse_average.py   for NO reaction field calculations
#python parse_average.py RF for RF calculations

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

#for each input output_ folder we need to split the shell width

rf = sys.argv[1]  # options RF and NO
inptf = glob.glob("output_*/*.csv")

#create sorted list for the width  and store the output foldre names
widths = []
outputs = []
ifiles = []
for inpt in inptf:
    width = float(inpt.split("_")[3].split("/")[0])
    #if width ==1.0 :
    #    continue
    if width in widths:
        pass
    else:
        widths.append(width)
    output = inpt.split("_")
    merger = "_".join([output[0],output[1],output[2]])
    if merger in outputs:
        pass
    else:
        outputs.append(merger)
    ifile = inpt.split("/")[1]
    if ifile in ifiles:
        pass
    else:
        ifiles.append(ifile)

widths.sort()

#now initialise all the dictionaries- remember to free them after each width has been done
#angles.csv  AVG_RF_dipole.csv  AVG_RF_int_dipole.csv  BWRF.csv  n_waters.csv  README  RF_charge_dipole.csv

BWRF = {} #-->BWRF.csv
BWRF_std = {}
RF_charge_dipole = {} #-->RF_charge_dipole.csv
RF_charge_dipole_std = {}
AVG_RF_dipole = {} #-->AVG_RF_dipole.csv
AVG_RF_dipole_std = {}
AVG_RF_int_dipole = {} #--> AVG_RF_int_dipole.csv
AVG_RF_int_dipole_std = {}
n_waters = {}#--> n_waters.csv
n_waters_std = {}

shell_list = []
#we will use it to go through later for printing ordered shells


#ccle throught each output for each width 
print("Reading all the output files...")

if  rf =="RF":
    print("reaction field calculations...")
    bwrf_keyword = "BWRF"
    rf_charge_dipole_keyword = "RF_charge_dipole"
    avg_rf_dipole_keyword = "AVG_RF_dipole"
    avg_rf_int_dipole_keyword = "INT_dipole"
elif rf=="NO":
    print("No reaction field calculations...")
    bwrf_keyword = "/Coulomb"
    rf_charge_dipole_keyword = "/charge_dipole.csv"
    avg_rf_dipole_keyword = "/avg_rot_charge_dipole.csv"
    avg_rf_int_dipole_keyword = "/avg_rot_charge_dipole_int.csv"
else:
    print("Wrong RF choice. It's possible RF - for reaction field calculations- or NO - for non reaction field calcualtions...")
    sys.exit(-1)


for width in widths:
    print("Processing width %.2f" % width)
    for output in outputs:
        print("Folder %s" % output)
        for ifile in ifiles:
            
            foldername = "_".join([output,"%.2f"%width])
            filename = "/".join([foldername,ifile])
            readfile = open(filename,"r").readlines()
            energy_type = filename.split("/")[1].split(".csv")[0]
            #recognize what file is
            if bwrf_keyword in filename:
                #BWRF
                #for the moment  angle.csv has no shell, do it 
                for lines in readfile:
                    shell = float(lines.split(",")[0])
                    
                    nrg = float(lines.split(",")[1])
                    stdnrg = float(lines.split(",")[2])
                    if shell in BWRF : 
                        BWRF[shell].append(nrg)
                        BWRF_std[shell]  = math.sqrt(BWRF_std[shell]**2 + stdnrg**2)
                    else:
                        BWRF[shell]=[]
                        BWRF[shell].append(nrg)
                        BWRF_std[shell] = stdnrg
                        shell_list.append(shell)
            

            elif rf_charge_dipole_keyword in filename:
                #RF_charge_dipole
                for lines in readfile:
                    shell = float(lines.split(",")[0])
                    nrg = float(lines.split(",")[1])
                    stdnrg = float(lines.split(",")[2])
                    if shell in RF_charge_dipole:
                        RF_charge_dipole[shell].append(nrg)
                        RF_charge_dipole_std[shell] = math.sqrt(RF_charge_dipole_std[shell]**2 + stdnrg**2)
                    else:
                        RF_charge_dipole[shell]=[]
                        RF_charge_dipole[shell].append(nrg)
                        RF_charge_dipole_std[shell] = stdnrg

            elif avg_rf_dipole_keyword in filename:
                for lines in readfile:
                    shell = float(lines.split(",")[0])
                    nrg = float(lines.split(",")[1])
                    stdnrg = float(lines.split(",")[2])
                    if shell in AVG_RF_dipole:
                        AVG_RF_dipole[shell].append(nrg)
                        AVG_RF_dipole_std[shell] = math.sqrt(AVG_RF_dipole_std[shell]**2 + stdnrg**2)
                    else:
                        AVG_RF_dipole[shell] = []
                        AVG_RF_dipole[shell].append(nrg)
                        AVG_RF_dipole_std[shell] = stdnrg

            elif avg_rf_int_dipole_keyword in filename:
                for lines in readfile:
                    shell =  float(lines.split(",")[0])
                    nrg = float(lines.split(",")[1])
                    stdnrg = float(lines.split(",")[2])
                    if shell in AVG_RF_int_dipole:
                        AVG_RF_int_dipole[shell].append(nrg)
                        AVG_RF_int_dipole_std[shell] = math.sqrt(AVG_RF_int_dipole_std[shell]**2 + stdnrg**2)
                    else:
                        AVG_RF_int_dipole[shell]=[]
                        AVG_RF_int_dipole[shell].append(nrg)
                        AVG_RF_int_dipole_std[shell] = stdnrg

            elif "n_waters" in filename:
                for lines in readfile:
                    shell = float(lines.split(",")[0])
                    n_wat = float(lines.split(",")[1])
                    stdwat = float(lines.split(",")[2])
                    if shell in n_waters:
                        n_waters[shell].append(n_wat)
                        n_waters_std[shell] =math.sqrt(n_waters_std[shell]**2 + stdwat)
                    else:
                        n_waters[shell]=[]
                        n_waters[shell].append(n_wat)
                        n_waters_std[shell]=stdwat
                    

            else:
                continue



    #once all the outputs are finished save the file
    #compute the average
    print("Saving file")

    #sort the shells
    shell_list.sort()
    
    #since the shells are common to all te dictionaries, we can cycle through 
    #BWRF.keys() and compute the average and std err
    outputfolder = "output_average_%.2f" % width
    if  not os.path.exists(outputfolder):
        os.makedirs(outputfolder)
    #print("Writing files...")
    BWRF_file = open("%s/BWRF.csv"%outputfolder,"w")#BWRF results
    RF_charge_dipole_file = open("%s/RF_charge_dipole.csv"%outputfolder,"w")#RF_charge_dipole
    AVG_RF_dipole_file = open("%s/AVG_RF_dipole.csv"%outputfolder,"w") # AVG_RF dipole
    AVG_RF_int_dipole_file = open("%s/AVG_RF_int_dipole.csv" % outputfolder,"w")#AVG_RF_dipole integrated
    waters_file = open("%s/waters.csv"%outputfolder,"w")

    for key  in shell_list:
        #output: shell, nrg, stddev
        BWRF_file.write("%.4f,%.8f,%.8f\n" % (key,np.mean(BWRF[key]), BWRF_std[key]))
        RF_charge_dipole_file.write("%.4f,%.8f,%.8f\n" %(key,np.mean(RF_charge_dipole[key]),RF_charge_dipole_std[key]))
        AVG_RF_dipole_file.write("%.4f,%.8f,%.8f\n"%(key,np.mean(AVG_RF_dipole[key]),AVG_RF_dipole_std[key]))
        AVG_RF_int_dipole_file.write("%.4f,%.8f,%.8f\n" %(key,np.mean(AVG_RF_int_dipole[key]),AVG_RF_int_dipole_std[key]))
        waters_file.write("%.4f,%.4f,%.4f\n"%(key,np.mean(n_waters[key]),n_waters_std[key]))

    BWRF_file.close()
    RF_charge_dipole_file.close()
    AVG_RF_dipole_file.close()
    AVG_RF_int_dipole_file.close()
    waters_file.close()
    #print("Width %.2f done" % width)

    #re-initialize the dictionaries
    BWRF = {} #-->BWRF.csv
    BWRF_std = {}
    RF_charge_dipole = {} #-->RF_charge_dipole.csv
    RF_charge_dipole_std = {}
    AVG_RF_dipole = {} #-->AVG_RF_dipole.csv
    AVG_RF_dipole_std = {}
    AVG_RF_int_dipole = {} #--> AVG_RF_int_dipole.csv
    AVG_RF_int_dipole_std = {}
    n_waters = {}#--> n_waters.csv
    n_waters_std = {}
    shell_list = []

print("Everything 's done")
