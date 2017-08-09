#Which is the minimum distance for each bond length?
import numpy as np

ifile =np.loadtxt("theta_dist.csv",delimiter=",",usecols=[0,1,2])
print(len(ifile[:,0]))
#
counter = 0 
distances = []

for val in ifile:
	bondlength = val[0]
	distance = val[2]

	if counter==0:
		tmp = bondlength
		counter+=1
		
	if tmp!=bondlength:
		#compute the mean
		print("Minimum distance for %.4f is %.4f" % (tmp,min(distances)))
		distances = []
		tmp = bondlength
		distances.append(distance)

	distances.append(distance)
