#AUG 2017 Stefano Bosisio
#script to collect the ion-dipole distance and the Ntimes (numer of particles N)
#the difference between theoretical energy adn computed energy. 
#This script collects initially all the *_deg/computed.dat/theory.dat files
#For each value of r it stores the energy for each lambda
#producing a dictionary: r-deg = nrg_theory, nrg_computed
#then for each r-deg pairs it computed the difference between theory and prediction
#taking the absolute value and averaging this difference over all the degrees
#for each average difference we can compute how many particle we have on the surface,
#given by the radius r and the water radius, 2.0 A
#multiply the numer of water molecules on the surface, N, and the average difference
#
#At the end, it write a file, r_theta.csv, with: r, N*<DU> 
#in this way we can know which value of r has the lowest, spherically weighted,
#difference w.r.t the theory, and find the correct dipole bond length (L_* folder)
#and the correct distance of agreement between theory and predictions
#Usage:
#python r_theta_extract.py
#output file:
#L, r, N<DU>, std.dev(N<DU>)
import sys
import os
import glob
import math

ifolders = glob.glob("../l_*")
ofile = open("r_theta.csv","w") #this will be the output file
starting_folder = os.getcwd()
bondlength = []

for folder in ifolders:
    bond = float(folder.split("/")[1].split("_")[1])#this is the dipole bond length
    bondlength.append(bond)
    
bondlength.sort()

for bond in bondlength:
    
    print("Processing %.3f" % bond)
    #change directory
    os.chdir("../l_%.3f" % bond)
    #do the routine

    icomp = glob.glob("*_deg/computed.dat") #take the computed energy 
    itheory = glob.glob("*_deg/theory.dat") #take the theory energy files

    computed = {} # dictionary for computed 
    theory = {} #dictionary for theory
    r_list = [] #a list for the distances  of the ion w.r.t the dipole, important sort
    deg_list = [] # degree list

    print("Processing computed energies...")
    for comp in icomp: #process the computed files
        #for each file read the degree
        deg = float(comp.split("/")[0].split("_")[0])
        if deg in deg_list:
            pass
        else:
            #print(deg)
            deg_list.append(deg) #append into deg list if not present
        #open the file
        ifile = open(comp,"r").readlines() #open  each computed.dat file
        #now extract r and energy for each file
        for lines in ifile:
            r = float(lines.split(",")[0]) #take the distance ion-dipole 
            nrg = float(lines.split(",")[1]) #take the energy
            if r in r_list:
                pass
            else:
                r_list.append(r) #append distance if not present
            #now create the dictionary
            key = "%.4f-%.4f" % (r,deg)  # create key like distance-degree ,e.g. 3.0000 - 0.0000
            if key in computed:
                computed[key].append(nrg) 
            else:
                computed[key]=[]
                computed[key].append(nrg)
        

    #do the same for theory, the keys are the same
    print("Processing the theoretical predictions...")
    for pred in itheory:
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
        

    #sort the lists so we can write ordered outputfile
    r_list.sort()
    deg_list.sort()

    avg = 0.0 #intiialise an average
    counter = 0 #initialise a counter
    print("Creating output file r_theta.csv...")
    for r in r_list:
        N_particles = math.ceil((r**2)/(4.0)) #particle on the sruface
        #this is the same of 4pi*r**2/4pi*2**2, where 2 is the radius of a water mol
        for deg in deg_list:
            key = "%.4f-%.4f" % (r,deg) #create the key and cycle through all the degrees for a specific r

            #now process each dictionary and give as output a file with 
            theory_nrg = theory[key][0] #theoretical energy
            comp_nrg = computed[key][0] #computed energy
            #compute the difference in energy, as an absolute value
            diff = (abs(theory_nrg - comp_nrg))
            #compute a running average and variances
            if counter==0:
                avg = diff
                var = 0 
                counter+=1
            else:
                tmp = avg
                avg = avg + (diff - (avg))/counter
                var = var + (diff - (tmp))*(diff-avg)
                counter+=1
        #now convert the variance into a standard deviation
        stddev = math.sqrt(var/(counter-1))
        #now tohave the average energy and std.dev on the sphere multiply the value by N
        en_sphere = N_particles*avg
        en_sphere_dev = N_particles*stddev
        #now write in the file
        #r,nrg
        ofile.write("%.4f,%.4f,%.8f,%.8f\n" % (bond,r,en_sphere,en_sphere_dev))
        #reset 
        avg = 0.0
        var = 0.0
        stddev = 0.0
        counter = 0
    os.chdir(starting_folder)

    print("Done!:)")




