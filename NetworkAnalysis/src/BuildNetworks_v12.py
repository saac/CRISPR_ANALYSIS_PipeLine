import numpy as np
import sys

def Snapshoot(datafile,t):
    dataset = np.genfromtxt(datafile, skip_header=1)
    dataslide = []
    for i in range(len(dataset)):
        if dataset[i][0] == t:
            dataslide.append(dataset[i])
    
    return dataslide


######################################################  



def BipartieNetwork(virus_file,bacteria_slide,virus_slide,slide):

    nombre = virus_file.split("data-phage.txt")
    output = open('Bipartite_Network_'+str(nombre[0])+'Time_'+str(slide)+'.txt',"w")

    for i in bacteria_slide:
        B,bstr = Spacers(i)
        #print "B = "+str(bstr)+"  V = "+str(vstr)+"\t"+"mathch = ",len(tmp),"  ",tmp
        for j in virus_slide:
            V,vstr = Protospacers(j)
            tmp = np.intersect1d(B, V)
            if len(tmp) > 0:
                #print "B = "+str(bstr)+"  V = "+str(vstr)+"\t"+"mathch = ",len(tmp),"  ",tmp
                output.write("B_%s\tV_%s\t%s\n" % (str(int(bstr)), str(int(vstr)),str(len(tmp))))
            
    output.close()
    
    
    
def BipartieMatrixT(virus_file,bacteria_slide,virus_slide,slide):

    nombre = virus_file.split("data-phage.txt")
    output = open('Bipartite_MATRIX_'+str(nombre[0])+'Time_'+str(slide)+'.txt',"w")
    
    output.write("\t")    
    for j in virus_slide:
        V,vstr = Protospacers(j)
        output.write("V_%s\t" % (str(int(vstr))))  
    output.write("\n")

    for i in bacteria_slide:
        B,bstr = Spacers(i)
        output.write("B_%s\t" % (str(int(bstr))))
        for j in virus_slide:
            V,vstr = Protospacers(j)            
            tmp = np.intersect1d(B, V)
            
            if len(tmp) > 0:
                output.write("%s\t" % (str(len(tmp))))
            else:
                output.write("0\t")
        output.write("\n")
            
    output.close()
    
    
def BipartieMatrix(virus_file,bacteria_slide,virus_slide,slide):

    nombre = virus_file.split("data-phage.txt")
    output = open('Bipartite_MATRIX_'+str(nombre[0])+'Time_'+str(slide)+'.txt',"w")
    
    output.write("\t")    
    for i in bacteria_slide:
        B,bstr = Spacers(i)
        output.write("B_%s\t" % (str(int(bstr))))
    output.write("\n")

    for i in virus_slide:
        V,vstr = Protospacers(i)
        output.write("V_%s\t" % (str(int(vstr))))
        for j in bacteria_slide:
            B,bstr = Spacers(j)           
            tmp = np.intersect1d(B, V)
            
            if len(tmp) > 0:
                output.write("%s\t" % (str(len(tmp))))
            else:
                output.write("0\t")
        output.write("\n")
            
    output.close() 
    
    
def BipartieInfectionMatrix(virus_file,bacteria_slide,virus_slide,slide):

    nombre = virus_file.split("data-phage.txt")
    output = open('BipartieInfection_MATRIX_'+str(nombre[0])+'Time_'+str(slide)+'.txt',"w")
    
    output.write("\t")
    TotBabun = 0.0
    #TotVabun = 0.0
    for i in bacteria_slide:
        B,bstr = Spacers(i)
        Babun = float(i[3])
        TotBabun = TotBabun + Babun
        output.write("B_%s\t" % (str(int(bstr))))
    output.write("\n")
    
    #for i in virus_slide:
        #Vabun = float(i[3])
        #TotVabun = TotVabun + Vabun
    
    #Vmax = MaxAbundance(virus_slide)
    #print Vmax    

    for i in virus_slide:
        V,vstr = Protospacers(i)
        Vabun = float(i[3])
        output.write("V_%s\t" % (str(int(vstr))))
        
        for j in bacteria_slide:
            B,bstr = Spacers(j)
            Babun = float(j[3])
            tmp = np.intersect1d(B, V)           
            
            if len(tmp) == 0:
                w = (Babun*Vabun)/TotBabun      
                #w = (Babun*Vabun)/(TotBabun*Vmax)
                #w = (Babun*Vabun)/(TotBabun*TotVabun)
                output.write("%s\t" % (str(w)))
                #output.write("B_%s\tV_%s\t%s\n" % (str(int(bstr)), str(int(vstr)), str(w)))                
            else:
                output.write("0\t")
        output.write("\n")
            
    output.close()     
    
    
    
def ProtospacersByVirus(virus_file,virus_slide,slide):

    nombre = virus_file.split("data-phage.txt")
    #output = open('Protospacers-by-virus_'+str(nombre[0])+'Time_'+str(slide)+'.csv',"w")
    output = open('Protospacers-by-virus_'+str(nombre[0])+'Time_'+str(slide)+'.txt',"w")
    
    Pts = []
    
    output.write("\t")    
    for i in virus_slide:
        V,vstr = Protospacers(i)
        for j in V:
            if (str(int(j)) not in Pts):
                Pts.append(str(int(j)))
                output.write("Ps_%s\t" % (str(int(j))))
                #output.write("%s \t" % (str(int(j))))
    output.write("\n")

    for i in virus_slide:
        V,vstr = Protospacers(i)
        output.write("V_%s\t" % (str(int(vstr))))
        for j in Pts:
            for k in V:
                if (str(int(k)) == j):
                    output.write("1\t")
                    break;
                else:
                    if (k == V[-1]):
                        output.write("0\t")
        output.write("\n")
            
    output.close()          
    
    
  
def SpacersByBacteria(bacteria_file,bacteria_slide,slide):

    nombre = bacteria_file.split("data-bact.txt")
    #output = open('Spacers-by-bacteria_'+str(nombre[0])+'Time_'+str(slide)+'.csv',"w")
    output = open('Spacers-by-bacteria_'+str(nombre[0])+'Time_'+str(slide)+'.txt',"w")
    
    Sps = []
    
    output.write("\t")    
    for i in bacteria_slide:
        B,bstr = Spacers(i)
        for j in B:
            if ((str(int(j)) not in Sps) and ((str(int(j)) != "-1"))):
                Sps.append(str(int(j)))
                output.write("Sp_%s\t" % (str(int(j))))
                #output.write("%s \t" % (str(int(j))))
    output.write("\n")        

    #print Sps

    for i in bacteria_slide:
        B,bstr = Spacers(i)
        output.write("B_%s\t" % (str(int(bstr))))
        #for j in range(len(Sps)):
        for j in Sps:
            for k in B:
                #bol = False
                #if (str(int(k)) == Sps[j]):
                #print "ProrSpacer = ",k,"\t","Spacer = ",j
                #print "last = ",B[-1]
                if (str(int(k)) == j):
                    output.write("1\t")
                    break;
                    #bol = True
                else:
                    #print k,"  ",B[-1]
                    #if (k == len(B)+1):
                    if (k == B[-1]):
                        output.write("0\t")
                        break;
        output.write("\n") 
    output.close()     
    
    
  
    
def TripartieNetwork(virus_file,bacteria_slide,virus_slide,slide):

    nombre = virus_file.split("data-phage.txt")
    output = open('Tripartite_Network_'+str(nombre[0])+'Time_'+str(slide)+'.txt',"w")
    
    Pts = []
    Sps = []

    for i in bacteria_slide:
        B,bstr = Spacers(i)
        for j in B:
            if ((str(int(j)) not in Sps) and ((str(int(j)) != "-1"))):
                Sps.append(str(int(j)))
        
    for i in virus_slide:
        V,vstr = Protospacers(i)
        for j in V:
            if (str(int(j)) not in Pts):
                Pts.append(str(int(j)))
        
    SpacersSet = np.union1d(Sps, Pts)
        
   
    for i in bacteria_slide:
        B,bstr = Spacers(i)
        I = np.intersect1d(B, SpacersSet)
        #print "B_"+str(bstr), I, len(I)
        for j in I:
            #print "B_"+str(bstr), j
            output.write("B_%s\tSP_%s\n" % (str(bstr), str(j)))
        #print "----------"
    
    for i in virus_slide:
        V,vstr = Protospacers(i)
        I = np.intersect1d(V, SpacersSet)
        #print "V_"+str(bstr), I, len(I)
        for j in I:
            #print "V_"+str(vstr), j
            output.write("V_%s\tSP_%s\n" % (str(vstr), str(j)))
        #print "----------"    
    
    output.close()
        
    
  
 ######################################################   

def SimilarityBacteriaNetwork(bacteria_file,bacteria_slide,slide):

    nombre = bacteria_file.split("data-bact.txt")
    output = open('Similarity-Bacteria_Network_'+str(nombre[0])+'Time_'+str(slide)+'.txt',"w")

    for i in bacteria_slide:
        B1,source_str = Spacers(i)
        for j in bacteria_slide:
            B2,target_str = Spacers(j)
            #print source_str, target_str
            if (source_str != target_str):
                tmp1 = np.intersect1d(B1, B2)
                index = np.argwhere(tmp1==-1)
                tmp = np.delete(tmp1, index)
                if len(tmp) > 0:
                    #print str(source_str)+"\t"+str(target_str)+"\t"+"mathch = ",len(tmp),"  ",tmp
                    output.write("%s\t%s\t%s\n" % ("B_"+str(int(source_str)), "B_"+str(int(target_str)),str(len(tmp))))
            
    output.close()
    


def SimilarityVirusNetwork(virus_file,virus_slide,slide):

    nombre = virus_file.split("data-phage.txt")
    output = open('Similarity-Phage_Network_'+str(nombre[0])+'Time_'+str(slide)+'.txt',"w")

    for i in virus_slide:
        V1,source_str = Protospacers(i)
        for j in virus_slide:
            V2,target_str = Protospacers(j)
            #print source_str, target_str
            if (source_str != target_str):
                tmp1 = np.intersect1d(V1, V2)
                index = np.argwhere(tmp1==-1)
                tmp = np.delete(tmp1, index)
                if len(tmp) > 0:
                    #print str(source_str)+"\t"+str(target_str)+"\t"+"mathch = ",len(tmp),"  ",tmp
                    output.write("%s\t%s\t%s\n" % ("V_"+str(int(source_str)), "V_"+str(int(target_str)),str(len(tmp))))
            
    output.close()


######################################################


def UNSimilarityVirusNetwork(virus_file,virus_slide,slide):

    nombre = virus_file.split("data-phage.txt")
    output = open('Un-Similarity-Phage_Network_'+str(nombre[0])+'Time_'+str(slide)+'.txt',"w")

    for i in virus_slide:
        V1,source_str = Protospacers(i)
        for j in virus_slide:
            V2,target_str = Protospacers(j)
            #print source_str, target_str
            if (source_str != target_str):
                A=set(V1)
                B=set(V2)
                #print A
                #print B
                if (-1 in A): A.remove(-1)
                if (-1 in B): B.remove(-1)
                tmp = (A.difference(B))
                if len(tmp) > 0:
                #if len(tmp) == 0:
                    #print str(source_str)+"\t"+str(target_str)+"\t"+"mathch = ",len(tmp),"  ",tmp
                    output.write("%s\t%s\t%s\n" % ("V_"+str(int(source_str)), "V_"+str(int(target_str)),str(len(tmp))))
            
    output.close()



def InfectionNetwork(virus_file,bacteria_slide,virus_slide,slide):

    nombre = virus_file.split("data-phage.txt")
    output = open('Infection_Network_'+str(nombre[0])+'Time_'+str(slide)+'.txt',"w")
    
    TotBabun = 0.0
    for i in bacteria_slide:
        Babun = float(i[3])
        TotBabun = TotBabun + Babun
    
    #print TotBabun
    
    for i in bacteria_slide:
        B,bstr = Spacers(i)
        Babun = float(i[3])
        for j in virus_slide:
            V,vstr = Protospacers(j)
            tmp = np.intersect1d(B, V)
            Vabun = float(j[3])
            if len(tmp) == 0:
                w = (Babun*Vabun)/TotBabun
                #output.write("B_%s\tV_%s\t%s,%s\n" % (str(int(bstr)), str(int(vstr)), str(Babun), str(Vabun)))
                output.write("B_%s\tV_%s\t%s\n" % (str(int(bstr)), str(int(vstr)), str(w)))                
            
    output.close()

######################################################


#def Protospacers(Vslide):
    #ptspacers = []
    #for i in range(len(Vslide)):
        #l = len(Vslide[i])
        #for j in range(4,l):
            #ptspacers.append(Vslide[0][j])
   
    #return ptspacers
    
def Protospacers(virion):
    ptspacers = []
    vir = int(virion[2])
    l = len(virion)
    for j in range(4,l):
        ptspacers.append(int(virion[j]))
   
    return ptspacers, vir 
                
           
def Protospacers2(virusarray,strain):
    ptspacers = []
    #for i in range(len(virusarray)):
    l = len(virusarray[strain])
    for j in range(4,l):
        ptspacers.append(int(Vslide[strain][j]))
   
    return ptspacers           
            
            

def Spacers(strain):
    spcrs = []
    l = len(strain)
    strn = int(strain[2])
    for j in range(4,l):
        spcrs.append(int(strain[j]))
   
    return spcrs, strn
            
            
def MaxAbundance(slide):
    abuns = []
    for i in slide:
        abuns.append(float(i[3]))
        
    maxA = np.nanmax(abuns)
    
    return maxA
                              
            
##################################################

#print "Hola"

virus_file = str(sys.argv[1])
bacteria_file = str(sys.argv[2])
slide = int(sys.argv[3])

virus_slide = Snapshoot(virus_file,slide)
bacteria_slide = Snapshoot(bacteria_file,slide)

print "----------- Bipartite: -----------------"
BipartieNetwork(virus_file,bacteria_slide,virus_slide,slide)
print "----------- Similarity Bacteria: -----------------"
SimilarityBacteriaNetwork(bacteria_file,bacteria_slide,slide)
print "----------- Similarity Virus: -----------------"
SimilarityVirusNetwork(virus_file,virus_slide,slide)
print "----------- Tripartite: -----------------"
TripartieNetwork(virus_file,bacteria_slide,virus_slide,slide)

print "#########################################################"

print "----------- Bipartie Matrix: -----------------"
BipartieMatrix(virus_file,bacteria_slide,virus_slide,slide)
print "----------- Protospacers by Virus Matrix: -----------------"
ProtospacersByVirus(virus_file,virus_slide,slide)
print "----------- Spacers by Bacteria Matrix: -----------------"
SpacersByBacteria(bacteria_file,bacteria_slide,slide)

print "#########################################################"

print "----------- DisSimilarity Virus: -----------------"
UNSimilarityVirusNetwork(virus_file,virus_slide,slide)


print "----------- Infection Network: -----------------"
InfectionNetwork(virus_file,bacteria_slide,virus_slide,slide)
print "----------- Bipartie Infection Matrix: -----------------"
BipartieInfectionMatrix(virus_file,bacteria_slide,virus_slide,slide)