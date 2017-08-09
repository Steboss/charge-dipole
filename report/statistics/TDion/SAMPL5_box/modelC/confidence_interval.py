#Script to compute the error bounds in R2 and MUE
import numpy as np
import scipy.stats
import glob
import sys


def compute_intervals(absolute, relative):
    #deal with the csv files
    abs_val, rel_val = read_subdata(absolute, relative)
    #now compute the errors
    compute_errors(abs_val, rel_val)

def read_subdata(absolute, relative):
    #open the absolute
    f = open(absolute, 'r')
    lines = f.readlines()
    abs_val = []
    for l in lines:
        val = float(l.split(",")[1])
        err = float(l.split(",")[2])
        abs_val.append(np.array([val,err]))
    f.close()

    abs_val = np.array(abs_val)
    #now do the same with the relative
    f = open(relative, 'r')
    lines = f.readlines()
    rel_val = []
    for l in lines:
        val = float(l.split(",")[1])
        err = float(l.split(",")[2])
        rel_val.append(np.array([val,err]))
    f.close()
    rel_val = np.array(rel_val)
    #print values
    print("Absolute free energies...")
    print(abs_val)
    print("Relative free energies...")
    print(rel_val)
    return (abs_val, rel_val)


def compute_errors(absolute, relative):
    r_dist = []
    mue_dist = []
    rsqrt_dist = []
    #now we are going to bootstrap each relative values, by using its error
    for i in range(1000):
        #here we take our relative values along with their erros and we
        #extract a new value by taking a random gaussian number between val-err
        relative_new = perturb(relative)
        #compute R2
        rsqrt,r = calculate_R2(relative_new[:,0],absolute[:,0])
        #compute MUE
        mue = calculate_mue(relative_new[:,0],absolute[:,0])
        mue_dist.append(mue)
        r_dist.append(r)
        rsqrt_dist.append(rsqrt)

    rsqrt_dist = np.array(rsqrt_dist)
    r_dist = np.array(r_dist)
    mue_dist = np.array(mue_dist)
    #now we need to know what are the max error bounds for r, r2 and mue
    #based on a 95% confidence interval
    (r_min, r_mean, r_max) = get_error_bounds(r_dist)
    (rsqrt_min, rsqrt_mean, rsqrt_max) = get_error_bounds(rsqrt_dist)
    (mue_min, mue_mean, mue_max) = get_error_bounds(mue_dist)
    print ("Observables are:")
    print ("R: %f < %f < %f " % (r_min, r_mean, r_max))
    print ("Rsqrt: %f < %f < %f " % (rsqrt_min, rsqrt_mean, rsqrt_max))
    print ("mue: %f < %f < %f " % (mue_min, mue_mean, mue_max))

    #write everything on a files
    outputfile = open("errorbounds.dat","w")
    outputfile.write("R: %f < %f < %f\n" % (r_min, r_mean, r_max))
    outputfile.write("Rsqrt: %f < %f < %f\n" % (rsqrt_min, rsqrt_mean, rsqrt_max))
    outputfile.write("mue: %f < %f < %f\n" % (mue_min, mue_mean, mue_max))
    outputfile.close()


def perturb(data):
    r""" generates new set of data based on gauss distribution
    Parameters
    ----------
    data : nd.array(shape(datapoints,2))
        first column holding actual data, second error on data

    """
    repeat = np.zeros(np.shape(data))

    count = 0
    for d in data:
        val = d[0]
        err = d[1]
        if err != 0.0:
            #print(val,err)
            val2 = np.random.normal(val, err)
        else:
            val2 = val
        repeat[count][0] = val2
        repeat[count][1] = err
        count = count + 1

    return repeat


def calculate_R2 ( series1, series2 ):
    r_value,p = scipy.stats.pearsonr(series1,series2)

    return r_value**2, r_value


def calculate_mue( series1, series2 ):

    sumdev = 0.0
    for x in range(0,len(series1)):
        sumdev += abs( series1[x] - series2[x] )
    sumdev /= len(series1)

    #print sumdev
    return sumdev


def get_error_bounds(data):
    #we will have a distribution of values for r,r2 and mue, thus from here we
    #can assess the max and min error based on 95% CI
    n = data.shape[0]
    #take the mean
    mean = np.mean(data)
    data_sorted = np.sort(data)
    #now from the distribution, sorted, compute the 5% and 95% values - basically
    #we are splitting the distribution at 5% and 95%
    x1= int(np.floor(0.05*n))
    x2= int(np.ceil(0.95*n))
    #so from here we will have the lower and upper bound
    lower_bound = data_sorted[x1]
    upper_bound = data_sorted[x2]
    return (lower_bound, mean, upper_bound)



if __name__ == '__main__':
    print ("========Analysis==========")
    print ("++++++++++++++Protocol+++++++++++++++++")
    #take abs and relative csv inputfiles
    absolute = sys.argv[1]
    relative = sys.argv[2]
    #go!
    compute_intervals(absolute, relative)
