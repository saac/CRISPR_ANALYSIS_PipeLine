import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.image as mpimg
#import matplotlib.animation as animation
#import time
import pandas as pd

import os
import sys

archivo = str(sys.argv[1])
nombrearchivo = str(sys.argv[3])

organism_abundance = np.genfromtxt(archivo, skip_header=1)

nombre = nombrearchivo.split("-abundance.txt")

figname = nombre[0]
print figname

label = np.unique(organism_abundance[:,2]).astype(int);
stacked_plot = np.zeros([int(sys.argv[2]), label.shape[0]])
new_timeRecord = 1;
timesOfRecord = np.array([])
time = np.array([])
for i in range(organism_abundance.shape[0]):
    if organism_abundance[i,0] == new_timeRecord:
        timesOfRecord = np.concatenate((timesOfRecord, [new_timeRecord]))
        time = np.concatenate((time, [organism_abundance[i,1]]))
        new_timeRecord = new_timeRecord + 1
    xlabel = int(organism_abundance[i,0] - 1)
    ylabel = np.argmax(organism_abundance[i,2] == label)
    stacked_plot[xlabel, ylabel] = organism_abundance[i,3]

    
t1 = 0
t2 = int(sys.argv[2])

beg = int(sys.argv[4])
end = int(sys.argv[5])

df2 = pd.DataFrame(stacked_plot[t1:t2])
df2.plot.area(stacked=True, legend=False, linewidth=0);
plt.title(figname+" Abundance")
plt.xlim(beg, end)
plt.savefig(str(beg)+"_to_"+str(end)+"_hrs_"+figname+"-abundance.png",dpi=300)
#plt.show()
plt.close()


