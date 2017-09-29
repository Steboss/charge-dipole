#AUG 2017 Stefano Bosisio
#Test on the reweigthing process for charging correction
#in this test I want to compute the amount of free energy corrections  by considering a regione within
#the simulated cubic box
#TODO: ADD TO SIRE SCRIPT, so infromation can be retrieved automatically
#Usage: ~/sire_pbc/sire.app/bin/python reweigthing_test.py
import os,sys, random
import math
import numpy as np
import mdtraj as md
from Sire.Tools.OpenMMMD import *
from Sire.Tools import Parameter, resolveParameters

dielectric = 82.0
one_over_four_pi_eps0=332.0637090025476

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


def BarkerWatts(distance,water_charge,rcutoff):

    nacharge = 1.0  #give the ion charge in e

    krf = (1/rcutoff**3)*((dielectric-1)/(2*dielectric+1)) #krf term
    pot = one_over_four_pi_eps0*(nacharge*water_charge)*((1/distance) + krf*distance**2)


    return pot

def BarkerWattsShifted(distance,water_charge,rcutoff):

    nacharge = 1.0  #give the ion charge in e

    krf = (1/rcutoff**3)*((dielectric-1)/(2*dielectric+1)) #krf term
    crf = (1/rcutoff)*((3*dielectric)/(2*dielectric + 1))
    pot = one_over_four_pi_eps0*(nacharge*water_charge)*((1/distance) + krf*distance**2 - crf)


    return pot


def computePotentials(frame,framenumb,top,rcutoff,box_edges):
    #Compute the BWRF for short cutoff, the one used in the simualtion
    #for big cutoff, at each iteraction detect the smallest box_edge and
    #compute the minimum images  convection
    nrgsBW = [] #list to store all the energy valus, sum up at the end BWRF
    nrgsBWshift = []
    #min_edge = min(box_edges[0])
    #print("Min edge length %.4f" % min_edge)


    coordinates = frame[0]
    #number of atoms
    natoms = len(coordinates)
    #sodium coodinates, should always be the first ones
    na_coords = (coordinates[0][0],coordinates[0][1],coordinates[0][2])

    at_count = 0 
    #now compute the ion-atoms
    for res in top.residues:
        if "HOH" in res.name:

            for at in res.atoms:
                #compute the com
                atname = at.name
                atidx = at.index
                #for each atom compute the distance
                dist = compute_distance(na_coords,coordinates[atidx])
                #now it distance is grater than > min_edge/2.0 rescale it
                if dist > (rcutoff/2.0):
                    dist = dist - rcutoff/2.0
                    if dist < rcutoff:
                        at_count+=1
                        
                        #compute the interaction
                        if "H" in atname:
                            atcharge = 0.417
                        elif "O" in atname:
                            atcharge = -0.834
                        else:
                            pass


                        #compute the BWRF energy here with the current cutoff#
                        potBW = BarkerWatts(dist,atcharge,rcutoff)
                        potBWshift = BarkerWattsShifted(dist,atcharge,rcutoff)
                        
                        nrgsBW.append(potBW)
                        nrgsBWshift.append(potBWshift)
                    else:
                        continue


        else:
            continue
    #print out
    #print("Cutoff %.4f  number of atoms %.2f" % (rcutoff,at_count))
	#sum up the contributions
    if at_count==0:
        nrgBW=sum(nrgsBW)
        nrgBWshift=sum(nrgsBWshift)
    else:
        nrgBW = (sum(nrgsBW))/at_count
        nrgBWshift = (sum(nrgsBWshift))/at_count
    return nrgBW,nrgBWshift


def getFreeEnergy(delta_nrgs):
    #print("temperature")
    ##print(temperature.val)
    free_nrg = FreeEnergyAverage(temperature.val)
    for nrg in delta_nrgs:
        free_nrg.accumulate(nrg)
    deltaG = free_nrg.average() * kcal_per_mol
    return deltaG

def resample(values):
    nvals = len(values)
    new_values = []
    for x in range(0,nvals):
        i = random.randint(0,nvals-1)
        new_values.append(values[i])
    return new_values


##############
#MAIN SCRIPT##
##############

topfile = sys.argv[1]  #take topology 
trajfile = sys.argv[2] #traj
mdtrajtop = md.load_prmtop(topfile) #this function allows the usage of SYSTEM.top
mdtrajdcd = md.open(trajfile,"r") #this function does not require huge memory allocation
nframes = len(mdtrajdcd) -1# this is necessary to compute the radius of cutoff

#import pdb; pdb.set_trace()
#now compute the  minimum box length

#compute the average radius of cutoff
minimum = 100000 #set a minumim (higher than the expected box length)
maximum = 0 # compute th emaximum radius
print("Computing the maximum box length")
#compute the minimum length to have a sphere inscribed a cube
#and the maximum which will be used to draw the sphere circumscribed the cube
for framenumber in range(0,nframes):#
    current, cell_lengths, angles = mdtrajdcd.read(n_frames=1) #for each frame extract the cell_lengths
    for length in cell_lengths[0]: #comparison
        if length< minimum:
            minimum = length
        if length> maximum:
            maximum = length
#use the maximum in order to encompass all the interactions
mdtrajdcd.seek(0) #reset to zero the frames

maximum = maximum*math.sqrt(3) #this is the max cutoff
#which can circumscribe da sphere around the "cubic" box
print("Max Radius of cutoff %.4f" % maximum)

#create a list of cutoff we want analyse 
#and add the max and min as well
sim_cutoffs= []
sim_cutoffs = np.linspace(2,maximum,20)
#sim_cutoffs.append(minimum)
print("All the cutoff to examine")
print(sim_cutoffs)

current_frame = 0

resultsBW = {}
resultsBWshift = {}
#create some dicts for delta nrgs
deltaBW = {}
deltaBWshift = {}

#now follow this way to compute 
#1) potentials avg
#2) avg DG
#block average of 5 snapshots each
#dictionaries:
dG = {}
dGshift = {}
dU = {}
dUshift = {}


while (current_frame <= nframes):
    #print ("Processing frame %s " % current_frame)
    #print ("CURRENT POSITION %s " % mdtrajdcd.tell() )
    frames_xyz, cell_lengths, cell_angles = mdtrajdcd.read(n_frames=1)
    #this can be computed once
    #nrgs_longBW,nrgs_longBWshift  = computePotentials(frames_xyz,current_frame,mdtrajtop,maximum,cell_lengths)

    for simc in sim_cutoffs:
        #print("Cutoff %.4f " % simc)
        #compute the BWRF potentials
        nrgs_shortBW,nrgs_shortBWshift = computePotentials(frames_xyz,current_frame,mdtrajtop,simc,cell_lengths)
        #print("Total potential %.4f"  % nrgs_shortBW)
        #print (system_longc.energy())
        #delta_nrgBW = (nrgs_longBW - nrgs_shortBW)
        #delta_nrgBWshift = (nrgs_longBWshift - nrgs_shortBWshift)

        if simc in resultsBW:
            #deltaBW[simc].append(delta_nrgBW)
            #deltaBWshift[simc].append(delta_nrgBWshift)
            resultsBW[simc].append(nrgs_shortBW)
            resultsBWshift[simc].append(nrgs_shortBWshift)
        else:
            #deltaBW[simc] = [delta_nrgBW]
            #deltaBWshift[simc] = [delta_nrgBWshift]
            resultsBW[simc] = [nrgs_shortBW]
            resultsBWshift[simc]=[nrgs_shortBWshift]

    if current_frame%10==0.0:
        #compute free energies
        for simc in sim_cutoffs:
            deltaG  = getFreeEnergy(resultsBW[simc])
            deltaGshift = getFreeEnergy(resultsBWshift[simc])
            meanU = np.mean(resultsBW[simc])
            meanUshift = np.mean(resultsBWshift[simc])
            #now store in the dicts:
            if simc in dG: 
                dG[simc].append(deltaG.value())
                dGshift[simc].append(deltaGshift.value())
                dU[simc].append(meanU)
                dUshift[simc].append(meanUshift)
            else:
                dG[simc]=[deltaG.value()]
                dGshift[simc] = [deltaGshift.value()]
                dU[simc] = [meanU]
                dUshift[simc] = [meanUshift]
        #now reset the dictionary
        resultsBW= {}
        resultsBWshift = {}
        current_frame+=1
    else:
        current_frame += 1

print("#######################################################################")
#print(dU)
#create an output file
ofile = open("electostratics.dat","w")

for simc in sim_cutoffs:
    #########DG############
    meandG = 0.0
    vardG = 0.0
    count_dG = 0
    for dGs in dG[simc]:
        if abs(dGs)!=float('Inf'):
            if count_dG ==0:
                vardG = 0
                meandG += dGs
                count_dG +=1
            else:
                tmp = meandG
                meandG = meandG + (dGs - meandG)/count_dG
                vardG = vardG + (dGs - tmp)*(dGs - meandG)
                count_dG+=1
    #compute the std dev
    if count_dG ==0:
        meandG = 999.0
        stddev = 999.0
    else:
        stddev = math.sqrt(vardG/(count_dG+1))

    dG[simc] = meandG,stddev
    #########DGshift############
    meandGshift = 0.0
    vardGshift = 0.0
    count_dGshift = 0
    for dGsh in dGshift[simc]:
        if abs(dGsh) !=float('Inf') :
            if count_dGshift ==0:
                vardGshift = 0
                meandGshift += dGsh
                count_dGshift +=1
            else:
                tmp = meandGshift
                meandGshift = meandGshift + (dGsh - meandGshift)/count_dGshift
                vardGshift = vardGshift + (dGsh - tmp)*(dGsh - meandGshift)
                count_dGshift+=1
    #compute the std dev
    if count_dGshift == 0:
        meandGshift = 999.0
        stddevdGshift = 999.0
    else:
        stddevdGshift = math.sqrt(vardGshift/(count_dGshift+1))

    dGshift[simc] = meandGshift,stddevdGshift
    #########Upot############
    meanUpot = 0.0
    varUpot = 0.0
    count_Upot = 0
    #print(np.mean(dU[simc]))
    for dUs in dU[simc]:
        if abs(dUs)!=float('Inf'):
            if count_Upot ==0:
                varUpot = 0
                meanUpot += dUs
                count_Upot +=1
            else:
                tmp = meanUpot
                meanUpot = meanUpot + (dUs - meanUpot)/count_Upot
                varUpot = varUpot + (dUs - tmp)*(dUs - meanUpot)
                count_Upot+=1
    #compute the std dev
    if count_Upot == 0:
        meanUpot = 999.0
        stddevUpot = 999.0
    else:
        stddevUpot = math.sqrt(varUpot/(count_Upot+1))

    dU[simc] = meanUpot,stddevUpot

    #########DGshift############
    meandUshift = 0.0
    varUshift = 0.0
    count_dUshift = 0
    for dUsh in dUshift[simc]:
        if abs(dUsh) !=float('Inf') :
            if count_dUshift ==0:
                varUshift = 0
                meandUshift += dUsh
                count_dUshift +=1
            else:
                tmp = meandUshift
                meandUshift = meandUshift + (dUsh - meandUshift)/count_dUshift
                varUshift = varUshift + (dUsh - tmp)*(dUsh - meandUshift)
                count_dUshift+=1
    #compute the std dev
    if count_dUshift ==0:
        meandUshift =999.0
        stddevUshift = 999.0
    else:
        stddevUshift = math.sqrt(varUshift/(count_dUshift+1))

    dU[simc] = meandUshift,stddevUshift

    print("Rc %.4f, dU %.4f+/-%.4f, dG %.4f+/-%.4f dUshift %.4f+/-%.4f, dGshift %.4f+/-%.4f" % (simc,dU[simc][0],dU[simc][1],dG[simc][0],dG[simc][1],dUshift[simc][0],dUshift[simc][1],dGshift[simc][0],dGshift[simc][1]))
    ofile.write("%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n" %(simc,dU[simc][0],dU[simc][1],dG[simc][0],dG[simc][1],dUshift[simc][0],dUshift[simc][1],dGshift[simc][0],dGshift[simc][1]))
        


