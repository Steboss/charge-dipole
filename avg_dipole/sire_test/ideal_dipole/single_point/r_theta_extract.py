#AUG 2017 Stefano Bosisio
#Script to extract theta, r and energy for each dipole length
#This script must be placed in all the l_* folders 
#At the end an output file will be created as:
#r,theta, energy_theory, energy_computed
#Usage: cd l_*; | python r_theta_extract.py
import sys
import os
import glob


icomp = glob.glob("*_deg/computed.dat")
ipred = glob.glob("*_deg/theory.dat")

ofile = open("r_theta.csv","w")
computed = {}
theory = {}
r_list = []
deg_list = []

for comp in icomp:
    #for each file read the degree
    deg = float(comp.split("/")[0].split("_")[0])
    if deg in deg_list:
        pass
    else:
        deg_list.append(deg)
    #open the file
    ifile = open(comp,"r").readlines()
    #now extract r and energy for each file
    for lines in ifile:
        r = float(lines.split(",")[0])
        nrg = float(lines.split(",")[1])
        if r in r_list:
            pass
        else:
            r_list.append(r)
        #now create the dictionary
        key = "%.4f-%.4f" % (r,deg)
        if key in computed:
            computed[key].append(nrg)
        else:
            computed[key]=[]
            computed[key].append(nrg)
    


#do the same for theory, the keys are the same
for pred in ipred:
    deg = float(pred.split("/")[0].split("_")[0])
    ifile = open(pred,"r").readlines()
    for lines in ifile:
        r = float(lines.split(",")[0])
        nrg = float(lines.split(",")[1])
        key = "%.4f-%.4f" % (r,deg)
        if key in theory:
            theory[key].append(nrg)
        else:
            theory[key]=[]
            theory[key].append(nrg)
    



#sort the lists
r_list.sort()
deg_list.sort()

for r in r_list:
    for deg in deg_list:
        key = "%.4f-%.4f" % (r,deg)
        #print("distance %.4f angle %.4f val %.4f" % (r,deg,theory[key][0],computed[key][0]))
        ofile.write("%.4f,%.4f,%.4f,%.4f\n" % (r,deg,theory[key][0],computed[key][0]))



