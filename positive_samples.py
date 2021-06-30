from openbabel import *
import pybel
import pandas as pd
import csv
from rdkit import Chem
import numpy as np
import os,glob
import math
import random

dt=pd.read_csv('Nuclear_Receptor_BinaryRelation.csv');
#sol=pd.read_csv('Nuclear_Receptor.csv');
data=[]
smiles=[]
filename = 'Nuclear_Receptor.csv'
with open(filename, 'r') as csvfile:
    sol= csv.reader(csvfile)
    for row in sol:
            #print(row[0],row[1],row[2])
        data.append(row)
print(len(data))
mat3=np.array(data)
print(len(mat3))
smiles=[]
count=0
for i in dt.Drug:
#    count=1
#    print(" ------------ ")
    
    for j in mat3:
#        print(count," --> ",i," == ",j[0])
#        count=count+1
        if i==j[0]:
            #print(i," == ",j[0])
            
            #print(count,j[0])
            smiles.append(j[2])
            break
#print("sms",len(smiles))
mols=[pybel.readstring("smi",x)for x in smiles]
fps1=[x.calcfp() for x in mols]
#for i in fps1:
#    print(i.bits)
#print(len(fps1))


X1=[]
for i in dt.Target:
    #print(i)
    fname=os.getcwd()+"\\Nuclear receptor-PSSM\\"+i+".pssm"
    #print(fname)
    f=open(fname,"r")
    con=f.read()
    f.close()
    arr=[]
    data=[]
    content=con.split("\n");

    for row in content:
        row.strip()
        row=row.replace("-"," -")
        row=row.replace("    "," ")
        row=row.replace("   ", " ")
        row=row.replace("  ", " ")
        data.append([x for x in row.strip().split(" ")])
    mat1=[]
    for i in range(3,len(data)-8):
    #for i in range(3,22):
            mat1.append(data[i][2:22])

    mat=np.array(mat1);
    #print(mat)
    m=len(mat1)
    n=len(mat1[0])
    # print(m,n)
    #count=count+1
    # print(count)
    #if count==2:
    #   print(mat)
    for i in range(m):
        for j in range(n):
            # if(len(mat[i][j])>3):
            #     print(count, "---")
            #     print(len(mat[i][j]))
            mat1[i][j]=1/(1+math.exp(-1*float(mat[i][j])))

    #display values normalised pssm to range (0-1)
    #for i in range(m):
    #    for j in range(n):
    #        print(mat1[i][j],end=" ")
    #    print()

#------Mean of column values----------
    mean=[]
    for i in range(n):
        sum=0
        for j in range(m):
            sum=sum+mat1[j][i]
        mean.append(sum/m)
    # print("Mean",mean)

#------PsePSSM----------
    l=m
    #l=3
    lamda=1
    g=[]
    sum=0
    result=[]

    for lamda in range (1,11):
    # g.clear()
        for i in range(n):
            sum=0
            for j in range(l-lamda):
                #print(mat1[j][i],mat1[j+lamda][i])
                sum=sum+((mat1[j][i]-mat1[j+lamda][i])**2)
            sum=sum/(l-lamda)
            g.append(sum)

    # print(g)
    # print(len(g))
    result=mean+g
    # print("----------result-------")
    # print(result)
    # print(len(result))
    X1.append(result)
# print(len(result),len(X1[-1]))
# print(result)
# print(X1)
for i in range(5):
    print("------",i,"------")
    print(X1[i])
print(len(X1))

# X1 - Target Positive Samples
# fps1 - Drug in positive sample


#------------Drug - Target Interaction + ve Samples---------------

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from networkx.algorithms import bipartite
from operator import add


g1=nx.Graph()

g1.add_nodes_from(dt.Drug,bipartite='Target')
g1.add_nodes_from(dt.Target,bipartite='Drugs')

# print(customers[0:10])
# g1.add_nodes_from(customers[:10],bipartite='customers')
# g1.add_nodes_from(products[:5],bipartite='products')

positive_dataset=[]
# print(g1.nodes(data=True))
#g1.add_edges_from([(0, 'a')
# , (0, 'b'), (0, 'c'), (1, 'a'),(1, 'b'),(1,'c'),(2,'a')])
print(len(X1),len(fps1))
for i in range(len(dt)):
    g1.add_edge(dt.Drug[i],dt.Target[i])
print()
#print(dataset)
#g1.add_edges_from((row['customers'], row['products']) for idx, row in df.iterrows()])
#print(nx.bipartite.maximum_matching(g1))
#print(len(nx.bipartite.maximum_matching(g1)))
print(nx.is_connected(g1))
#A, B = bipartite.sets(g1)
#pos = dict()
#pos.update( (n, (1, i)) for i, n in enumerate(A) ) # put nodes from X at x=1
#pos.update( (n, (2, i)) for i, n in enumerate(B) ) # put nodes from Y at x=2
nx.draw(g1 ,with_labels=True)
plt.show()

#dataset=[[x,y.bits] for x in X for y in fps]
#print(dataset) 
x=[]
#print(len(dataset))
for i in range(len(X1)):
    #print(fps1[i].bits)
    x=[]
    x=[X1[i],fps1[i].bits]
    positive_dataset.append(x)
#print(positive_dataset[0:2])
#print(len(positive_dataset))

