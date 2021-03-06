#Aug2017 Stefano Bosisio
#Master script to run multiple times orientation-distribution_function.py
#Usage:
#python master.py traj*  TOP  CRD 

import sys,os
import math
import numpy as np
import time
import mdtraj as md
from parmed.amber import *
import multiprocessing
from multiprocessing import Pool, Process
from functools import partial

######
#MAIN#
######

#we know that the total number of frames is 20'000
#so we start with 8 process
try:
    traj = sys.argv[1]
    top = sys.argv[2]
    width = sys.argv[3]
    reactionfield = sys.argv[4]  # options: RF or NO
except:
    traj ="traj000000001.dcd"
    top = "SYSTEM.top"
    width = 1.0
    reactionfield = "RF"

nproc = 8
nframes = 20000
print("Dividing %d frames into %d chunkcs" %(nframes,nproc))
chunks = np.linspace(0,nframes,nproc+1)

chunks_list = []

for i,val in enumerate(chunks,0):
    if i==0:
        start = int(val)
        end = int(chunks[i+1])
    elif i==len(chunks)-1:
        break
    else:
        start = int(val+1)
        end = int(chunks[i+1])
    chunks_list.append([start,end])


print("Sending data to orientation_distribution_function.py script...")
for i,val in enumerate(chunks_list,0):
    start = val[0]
    end   = val[1]
    ok_file = "foo_%d.out" % i
    err_file = "foo_%d.err" % i
    if reactionfield == "RF":
        cmd = "nohup python orientation_distribution_function_RF.py %s %s %s %d %d  > %s 2> %s </dev/null &" %(traj,top,width,start,end,ok_file,err_file)
        print(cmd)
        os.system(cmd)
    elif reactionfield=="NO":
        cmd = "nohup python orientation_distribution_function.py %s %s %s %d %d > %s 2> %s </dev/null &" % (traj,top,width,start,end,ok_file,err_file)
        print(cmd)
        os.system(cmd)
    else:
        print("ERROR! The reaction field methods is wrong!")
        print("Choices: RF for reaction field or NO for no reaction field")
        sys.exit(-1)



