
import sys
import math as mth
import numpy as np

import matplotlib.pyplot as plt

from pylab import *


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


def DiversityIndices2(datafile,t1,t2):
    
    t_list = []
    R_list = []
    H_list = []
    E_list = []

    t0 = t1
    
    dataset = np.genfromtxt(datafile, skip_header=1)
    dataslide = []
    
    tmp1=0
    tmp2=0
    
    for i in range(1,len(dataset)):
        if dataset[i][0] == t1:
            tmp1 = i+2
            break;
    for i in range(tmp1,len(dataset)):            
        if dataset[i][0] == t2:
            tmp2 = i+2
            break;
        
    t1 = tmp1
    t2 = tmp2
    
    for i in range(t1,t2):        
        if float(dataset[i][0]) != float(t0):
            
            R = len(dataslide)
            H = Shannon(dataslide)
            
            E = (H/(mth.log(R)))
            
            t_list.append(t0)
            R_list.append(R)
            H_list.append(H)
            E_list.append(E)
            
           
            dataslide = []
            t0=t0+1
            
        dataslide.append(dataset[i])
        
    if dataset[i][0] == dataset[len(dataset)-1][0]:
        
        R = len(dataslide)
        H = Shannon(dataslide)
        E = (H/(mth.log(R)))
        
        t_list.append(t0)
        R_list.append(R)
        H_list.append(H)
        E_list.append(E)
        
       
        dataslide = []

    
    return t_list, R_list, H_list, E_list



def DiversityIndicesPlot(datafile,t1,t2):
    
    data = open(datafile).readlines()[1:]
    
    T = []
    R = []
    H = []
    E = []
    
    #for row in data:
        #row=row.split("\t")
        ##print row[0], row[1], row[2], row[3]
        #T.append(float(row[0]))
        #R.append(float(row[1]))
        #H.append(float(row[2]))
        #E.append(float(row[3]))
        ##print i
        
    for i in range(t1-1,t2):
        #print data[i]
        row=data[i].split("\t")
        #print row[0], row[1], row[2], row[3]
        T.append(float(row[0]))
        R.append(float(row[1]))
        H.append(float(row[2]))
        E.append(float(row[3]))
        ##print i        
        
        
    #for i in range(0,150):
        #print T[i]
        
    return T, R, H, E



#def PlotData(Tempo,Values,t1,t2):
    
    #V = []
    #T = []
    
    #for t in range(t1-1,t2):
        #T.append(float(Tempo[t]))
        #V.append(float(Values[t]))
    

    #plot(T, V, '.', label='Virus')
    ##plot(Tb, Rb, '-', label='Bacteria')
    #grid(True)
    #legend(loc='upper left')
    #title('Diversity Indices ')
    #ylabel('Richness  $S$')
    
    #show()
    #close()

####################################################################################

##python DiversityIndices9.py mu5e-7_initialDiffDp2_S8P10_data-phage.txt mu5e-7_initialDiffDp2_S8P10_data-bact.txt
##python DiversityIndices9.py mu5e-7_initialDiffDp2_S8P10_data-phage.txt mu5e-7_initialDiffDp2_S8P10_data-bact.txt

##data_file = str(sys.argv[1])
##DiversityIndices0(data_file)


##################################

#data_file = str(sys.argv[1])
#T, R, H, E = DiversityIndicesPlot(data_file,100,200)
#print T


##virus_file = str(sys.argv[1])
##Tv,Rv,Hv, Ev = DiversityIndices(virus_file)

##PlotData(Tv,Rv,1,2500)

##############################################


#virus_file = str(sys.argv[1])
#bacteria_file = str(sys.argv[2])

#archivo = bacteria_file.split(".txt")
#nombre = archivo[0].split("_")
#output = nombre[0]+nombre[1]+nombre[2]+nombre[3]

#outputname = output+"_DiversityIndices.png"
##print outputname

#Tv,Rv,Hv, Ev = DiversityIndices(virus_file)
#Tb,Rb,Hb, Eb = DiversityIndices(bacteria_file)


#subplot(311) # Grafica uno
#plot(Tv, Rv, '.', label='Virus')
#plot(Tb, Rb, '.', label='Bacteria')
#grid(True)
#legend(loc='upper left')
#title('Diversity Indices '+output)
#ylabel('Richness  $S$')

#subplot(312) # Grafica dos
#plot(Tv, Hv, label='Virus')
#plot(Tb, Hb, 'r', label='Bacteria')
#grid(True)
#legend(loc='upper left')
##xlabel('time $t$ (hr)')
#ylabel('ShannonIndex  $H$')

#subplot(313) # Grafica dos
#plot(Tv, Ev, 'b', label='Virus')
#plot(Tb, Eb, 'r', label='Bacteria')
#grid(True)
#legend(loc='upper left')
#xlabel('time $t$ (hr)')
#ylabel('Evenness  $E_H$')


#savefig(outputname,dpi=300)
##show()
#close()

