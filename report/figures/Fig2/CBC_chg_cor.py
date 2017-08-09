#Plot of CBC molecules with charging corrections
import matplotlib.pyplot as plt
from numpy import *
import math, os, sys, re
import seaborn as sbn
sbn.set_style("whitegrid")

#select seaborn color palette to color all the bars
colors =sbn.color_palette()#("hls",16)
#now create a list of number to locate all the free energies
#we will have on the x-axis the molecule - to be modified in powerpoint -
ind = arange(10) -0.75 #add a -0.25 so we center the plot
#width of each bar
width = 0.25 #which is the same as the spacing
fig = plt.figure(figsize=[10,7])
ax = fig.add_subplot(111)

#BARS:

mod_c = [-6.57,1.24,1.46,-10.31,-14.43,-16.93,-13.55,-1.97,-14.61,-17.07]
err_c = [0.05, 0.05, 0.09, 0.06, 0.10, 0.13, 0.12, 0.35, 0.08, 0.10]
mod_d = [-0.71,6.79,19.53,-10.21,-12.51,-10.69,-23.54,-3.15,-14.00,-16.86]
err_d = [0.04,0.05,0.10,0.51,0.33,0.02,0.95,0.93,0.38,0.73]
exp = [-5.80,-2.50,-4.2,-7.54,-8.53,-8.64,-5.17,-6.17,-7.40,-10.35]
err_exp = [0.03,0.07,0.03,0.03,0.05,0.05,0.02,0.04,0.02,0.02]

for i,val in enumerate(mod_c,0):
    print("Computing bar %d" % i)
    bar_c = ax.bar(ind[i],val,width,color=colors[0],yerr=err_c[i],\
                    error_kw = dict(elinewidth=2,ecolor="black"))
    bar_d = ax.bar(ind[i]+width,mod_d[i],width,color=colors[1],yerr=err_d[i],\
                    error_kw = dict(elinewidth=2,ecolor="black"))
    bar_exp = ax.bar(ind[i]+2*width,exp[i],width,color=colors[2],yerr=err_exp[i],\
                    error_kw = dict(elinewidth=2,ecolor="black"))


#ax.set_ylim()
#plt.yticks(arange(start,end,pace))
#set the fontsize of the y-axis ticks
plt.yticks(fontsize=15)
ax.set_ylabel("Computed $\Delta G^\circ_\mathrm{bind}$ / kcal mol$^{-1}$", fontsize=20)

plt.xticks(ind+0.37)
tick_name = ["G1","G2","G3","G4","G5","G6","G7","G8","G9","G10"]
#assign the ticks name
xTickMarks=[str(x) for x in tick_name]
#set the size
xtickNames = ax.set_xticklabels(xTickMarks)
#ax.xaxis.tick_top()
#remember to first put the labels on top of the plot and then fontsize
#otherwise it won't work- it may be a bug in matplotlib
plt.xticks(fontsize=20)
#now set the ab,c,d,e,f,g,h on top of the plot
#ax.xaxis.set_label_position("top")
#grid on plot
plt.grid(True)
#plt.tight_layout()
#legend here we need to fix the labels
ax.legend((bar_c,bar_d,bar_exp),("$\Delta G^\circ_\mathrm{bind}$","$\Delta G^\circ_\mathrm{bind-cor}$","$\Delta G^\circ_\mathrm{EXP}$"),fontsize=20,bbox_to_anchor=(0.95, 1.10),\
ncol=3, fancybox=True,)

plt.savefig("sire_histogram.pdf")#,dpi=300,transparent=True)
#plt.savefig("sire_histogram.pdf",dpi=300)
