import os,sys


nots = []
for root, dirs, files in os.walk(os.getcwd()):
    if "l_" in root:
        if len(root.split("/")) == 10 :
            if "inputfiles" in root or ".ipynb_checkpoints" in root: 

                continue
            
            else:
                nots.append(root.split("/")[-1])
#print(nots)
for val in nots:
    cmd = "nohup ipython notebook %s/comparison_rf_theory.ipynb >/dev/null 2>&1 &" % val
    print(cmd)
    os.system(cmd)
