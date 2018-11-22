
import sys
import math as mth
import numpy as np
#import matplotlib.pyplot as plt

from pylab import *

#python TimeSeries.py mu5e-7_initialDiffDp2_S8P10_time-series-data.txt


   
    
def TimeSeries(datafile):
    
    data = open(datafile).readlines()[1:]
    
    T = []
    Vd = []
    Bd = []
    PDI = []
    
    for row in data:
        row=row.split(" ")
        #print row[0], row[4], row[5], row[7]
        T.append(float(row[0]))
        Bd.append(float(row[4]))
        Vd.append(float(row[5]))
        PDI.append(float(row[8]))
        #print i
    

    return T, Vd, Bd, PDI
    #return T, Vd, Bd

    


####################################################################################

#python TimeSeries.py mu5e-7_initialDiffDp2_S8P10_time-series-data.txt

data_file = str(sys.argv[1])
#T, Vd, Bd = TimeSeries(data_file)
T, Vd, Bd, PDI = TimeSeries(data_file)

#print len(T)
#print len(Vd)
#print len(Bd)
#print len(PDI)

archivo = data_file.split(".txt")
outputname = archivo[0]+".png"
#print outputname

subplot(311) # Grafica uno
plot(T, Vd, 'r', label='Virus')
grid(True)
legend(loc='upper right')
title('Density and PDI')
ylabel('Density')

subplot(312) # Grafica dos
plot(T, Bd, 'b', label='Bacteria')
grid(True)
legend(loc='upper right')
ylabel('Density')

subplot(313) # Grafica dos
plot(T, PDI, 'g')
grid(True)
#legend(loc='upper left')
xlabel('time $t$ (hr)')
ylabel('PDI')



###################################################################################################################

#subplot(411) # Grafica uno
##plot(T, Vd, 'r', label='Virus')
#plot(T, Vd, 'r')
##plot(T, Bd, 'b', label='Bacteria')
#grid(True)
#legend(loc='upper right')
##title('Density ')
#ylabel('Virus Density')

#subplot(412) # Grafica dos
#plot(T, Bd, 'b')
##plot(T, Hb, 'r', label='Bacteria')
#grid(True)
#legend(loc='upper right')
#ylabel('Bacteria Density')

#subplot(413) # Grafica dos
#plot(T, Bd, label='Bacteria')
#plot(T, Vd, label='Virus')
#grid(True)
#legend(loc='upper left')
#ylabel('Virus vs Bacteria Density')

#subplot(414) # Grafica dos
#plot(T, PDI, 'g')
#grid(True)
#legend(loc='upper left')
#xlabel('time $t$ (hr)')
#ylabel('PDI')

####################################################

#subplot(414) # Grafica dos
#plot(T, Vd, 'g')
#grid(True)
#legend(loc='upper left')
#xlabel('time $t$ (hr)')
#ylabel('PDI')


savefig(outputname,dpi=300)
#show()
close()

