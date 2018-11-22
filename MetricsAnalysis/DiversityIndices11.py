
import sys
import math as mth
import numpy as np
#import matplotlib.pyplot as plt

from pylab import *

#python DiversityIndices9.py mu3e-7_initialDiffDp20_S8P10_data-phage.txt mu3e-7_initialDiffDp20_S8P10_data-bact.txt 


def DiversityIndices0(datafile):
    
    nombre = datafile.split(".txt")
    
    output = open(str(nombre[0])+'_DiversityIndices.txt',"w")
    output.write("t \t Richness \t Shannon Index \n")
    
    dataset = np.genfromtxt(datafile, skip_header=1)
    dataslide = []
    
    t = 1
    
    for i in range(len(dataset)):
        if dataset[i][0] != t:
            R = len(dataslide)
            H = Shannon(dataslide)
            output.write("%s \t %d \t %f \n" % (str(int(t)), R, H))
            dataslide = []
            t=t+1
            
        dataslide.append(dataset[i])
        
    if dataset[i][0] == dataset[len(dataset)-1][0]:
        R = len(dataslide)
        H = Shannon(dataslide)
        output.write("%s \t %d \t %f \n" % (str(int(t)), R, H))
        dataslide = []

    output.close()
    
    
def DiversityIndices(datafile):
    
    t_list = []
    R_list = []
    H_list = []
    E_list = []
    
    nombre = datafile.split(".txt")
    output = open(str(nombre[0])+'_DiversityIndices.txt',"w")
    output.write("t \t Richness \t Shannon Index \t Evenness \n")
    
    dataset = np.genfromtxt(datafile, skip_header=1)
    dataslide = []
    
    t = 1
    
    for i in range(len(dataset)):
        if dataset[i][0] != t:
            
            R = len(dataslide)
            H = Shannon(dataslide)
            
            E = (H/(mth.log(R)))
            
            t_list.append(t)
            R_list.append(R)
            H_list.append(H)
            E_list.append(E)
            
            output.write("%s \t %d \t %f \t %f \n" % (str(int(t)), R, H, E))
            
            dataslide = []
            t=t+1
            
        dataslide.append(dataset[i])
        
    if dataset[i][0] == dataset[len(dataset)-1][0]:
        
        R = len(dataslide)
        H = Shannon(dataslide)
        E = (H/(mth.log(R)))
        
        t_list.append(t)
        R_list.append(R)
        H_list.append(H)
        E_list.append(E)
        
        output.write("%s \t %d \t %f \t %f \n" % (str(int(t)), R, H, E))
        
        dataslide = []

    output.close()
    
    return t_list, R_list, H_list, E_list


def Shannon(dataslide):
    P=0
    H=0
    for j in dataslide:
        P=P+j[3]
    for j in dataslide:
        p = j[3]/P
        H=H+(p*(mth.log(p)))
    if (H!=0):H=(-1)*H
    
    return H



####################################################################################

#python DiversityIndices9.py mu5e-7_initialDiffDp2_S8P10_data-phage.txt mu5e-7_initialDiffDp2_S8P10_data-bact.txt
#python DiversityIndices9.py mu5e-7_initialDiffDp2_S8P10_data-phage.txt mu5e-7_initialDiffDp2_S8P10_data-bact.txt

#data_file = str(sys.argv[1])
#DiversityIndices0(data_file)

virus_file = str(sys.argv[1])
bacteria_file = str(sys.argv[2])

archivo = bacteria_file.split(".txt")
nombre = archivo[0].split("_")
output = nombre[0]+nombre[1]+nombre[2]+nombre[3]

outputname = output+"_DiversityIndices.png"
#print outputname

Tv,Rv,Hv, Ev = DiversityIndices(virus_file)
Tb,Rb,Hb, Eb = DiversityIndices(bacteria_file)


subplot(311) # Grafica uno
plot(Tv, Rv, '.', label='Virus')
plot(Tb, Rb, '.', label='Bacteria')
grid(True)
legend(loc='upper left')
title('Diversity Indices '+output)
ylabel('Richness  $S$')

subplot(312) # Grafica dos
plot(Tv, Hv, label='Virus')
plot(Tb, Hb, 'r', label='Bacteria')
grid(True)
legend(loc='upper left')
#xlabel('time $t$ (hr)')
ylabel('ShannonIndex  $H$')

subplot(313) # Grafica dos
plot(Tv, Ev, 'b', label='Virus')
plot(Tb, Eb, 'r', label='Bacteria')
grid(True)
legend(loc='upper left')
xlabel('time $t$ (hr)')
ylabel('Evenness  $E_H$')


savefig(outputname,dpi=300)
#show()
close()

