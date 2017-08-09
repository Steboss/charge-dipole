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

mod_c = [-8.52,-1.93,-6.00,-11.28,-13.77,-17.37,-15.45,-4.31,-15.45,-17.93]
err_c = [0.02,0.12,0.35,0.59,0.13,0.07,0.92,0.81,0.52,0.04]
mod_d = [-7.98,-1.83,-4.52,-11.24,-13.57,-16.67,-14.79,-4.09,-15.16,-17.57]
err_d = [0.03,0.01,0.24,0.61,0.14,0.40,1.09,0.73,0.46,0.04]
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
