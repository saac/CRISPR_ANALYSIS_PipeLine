# -*- coding: utf-8 -*-

import sys
import networkx as nx
import pickle

# RUN # python Net2Infomap8.py RedPrueba.sif tab nIn


path = sys.argv[1]
sep = sys.argv[2]
form = sys.argv[3]

D = False
#direction = sys.argv[5]
direction = sys.argv[4]  ## o corte para subcoms!!
if (direction == 'D'): D = True

diccionario = {}
diccionarioB = {}
genesDIC = []


#################################### GENERA Red ##################################################

if (D):
  G = nx.DiGraph()
  print "Componentes 'solo para redes  No-Dirigidas!'"
else:
  G = nx.Graph()

name = path.split('.')
output = name[0]+'.net'
SALIDA = open(output,"w")

if (sep == "cma"):
  sep = ","
elif (sep == "tab"):
  sep = "\t"
elif (sep == "spa"):
  sep = " "  

interacciones = open(path).readlines()
#print interacciones[0]

k = 1
tmp="---"
for i in interacciones:
    genes  = i.split(sep)

    a=genes[0].strip()
    b=genes[1].strip()
    c=genes[2].strip()
    
    #print a,b,c
  
    if (form == "nIn"):
      if (genes[0].strip() not in genesDIC): genesDIC.append(genes[0].strip())
      if (genes[2].strip() not in genesDIC): genesDIC.append(genes[2].strip())
      
    elif (form == "nnI"):      
      if (genes[0].strip() not in genesDIC): genesDIC.append(genes[0].strip())
      if (genes[1].strip() not in genesDIC): genesDIC.append(genes[1].strip())  
  
#################################### puntoNET ################# 
SALIDA.write("%s %s \n" % ("*Vertices", len(genesDIC)))


i = 1
for n in genesDIC:
    diccionario[n] = i
    diccionarioB[i] = n
    SALIDA.write('%s "%s" \n' % (i,n))
    i+=1

#print genesDIC
#print G.nodes()

#for i in diccionario:
#print diccionario
    
for i in interacciones:
    campos  = i.split(sep)
    
    #print campos

    a=campos[0].strip()
    b=campos[1].strip()
    c=campos[2].strip()
    
       
    if (form == "nIn"):
      if (b == "pp"):
	b = 1
      else:
	b=float(b)
	
      if (diccionario[a] == diccionario[c]): print a,c	
      if ((a in diccionario and c in diccionario) and (diccionario[a] != diccionario[c]) and 
	  ((G.has_edge(diccionario[a],diccionario[c]) == False) or (G.has_edge(diccionario[c],diccionario[a]) == False))):
      #print G.has_edge(diccionario[a],diccionario[c])
      #print G.has_edge(diccionario[c],diccionario[a])
      #print diccionario[a],diccionario[c]
	G.add_edge(diccionario[a],diccionario[c],w=b)
    
    elif (form == "nnI"):
      if (c == "pp"):
	c = 1
      else:
	c=float(c)
	
	#print c
	
	
	
      if (diccionario[a] == diccionario[b]): print a,b	
      if ((a in diccionario and b in diccionario) and (diccionario[a] != diccionario[b]) and 
	  ((G.has_edge(diccionario[a],diccionario[b]) == False) or (G.has_edge(diccionario[b],diccionario[a]) == False))):
      #print G.has_edge(diccionario[a],diccionario[c])
      #print G.has_edge(diccionario[c],diccionario[a])
      #print diccionario[a],diccionario[c]
	G.add_edge(diccionario[a],diccionario[b],w=c)	
	
#print G.nodes()	
	  
if (not D):
  #G = G.to_undirected()
#print G.nodes()

#################################### Componentes ##################################################

  if not nx.is_connected(G):
    L = sorted(nx.connected_components(G))
    a = len(L)
    b = len(L[0])
    c = len(L[a-1])
    print "Se tienen",a,"componentes!"
    print "El mayor es de",b,"nodos."
    print "El menor es de",c,"nodos."
    tamano = int(raw_input("Ingresa el tamaño minimo de los componentes a considerar  "))
    comp_num = 1
    for g in nx.connected_components(G):
      #if len(g) >= int(sys.argv[4]):
      #if len(g) >= int(sys.argv[5]):
      if len(g) >= int(tamano):      
	H = G.subgraph(g)
	dictEdges = {}
	comp_name = "Componente"+str(comp_num)
	OUT = open(comp_name+'.net',"w")
	OUT.write("%s %s \n" % ("*Vertices", H.number_of_nodes()))
	k = 1

	for i in H.nodes():
	  OUT.write('%s "%s" \n' % (k,diccionarioB[i]))
	  dictEdges[i] = k
	  k = k+1
	  
	OUT.write("%s %s \n" % ("*Edges", H.number_of_edges()))
	for e in H.edges():
	    ed = H.get_edge_data(*e)
	    OUT.write("%s %s %s \n" % (dictEdges[e[0]],dictEdges[e[1]],ed['w']))
	    
	pick = comp_name+'.pickle'
	pickle.dump(H, open(pick, 'w'))

	comp_num = comp_num + 1
  
################################################################################################################ 
  
#################################### puntoNET ##################################################  
SALIDA.write("%s %s \n" % ("*Edges", G.number_of_edges()))
for e in G.edges():
    ed = G.get_edge_data(*e)
    SALIDA.write("%s %s %s \n" % (e[0],e[1],ed['w']))
    
picklename = name[0]+'.pickle'
pickle.dump(G, open(picklename, 'w'))