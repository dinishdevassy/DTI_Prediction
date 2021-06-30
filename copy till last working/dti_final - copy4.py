# ---------------Target---------------

import os,glob
import numpy as np
import math
import random

#f=open("DrugBank_Approved protein names.txt","r")
#t=f.read()
#f.close()
#vals=t.split("\n")
#print(len(vals))

# print(os.getcwd())
X=[]
Y=[]
count=0
targetnames=[]
# fname=os.getcwd()+"\\DrugBank_Approved-PSSM\\spP45059.pssm"
#for fname in glob.glob(os.getcwd()+"\\DrugBank_Approved-PSSM\\*.pssm"):
for fname in glob.glob(os.getcwd()+"\\Nuclear receptor-PSSM\\*.pssm"):
    #print(fname)
    # print(fname[-13:-5])
    f1=fname[-12:-5]
    f1=f1.replace("\\","")
    targetnames.append(f1)
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

    #print(row)
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
    count=count+1
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
    # change line code below
    #X.append(str(result))
    X.append(result)
# print(len(X))
# print(len(result),len(X[-1]))
# print(result)
# print(X)
for i in range(5):
    print("------",i,"------")
    print(X[i])

#     print("-------------------- ",i," ----------------------")




#------------Drug---------------


filename = "drugsubset.fpt"

# Accumulate the union of all the fingerprint bits
union_bits = 0

# Extract only the hex value for each compound fingerprint
# This assumes the fingerprint file is small enough to fit into memory!
# (The parser was easier to write this way.)
count=1
fingerprint=[]
for block in open(filename).read().split("\n>"):
    union_bits = 0
    header, remainder = block.split("\n", 1)
    #print(header)
    if remainder.startswith("Possible superstructure"):
        _, remainder = remainder.split("\n", 1)

    # Remove spaces and newlines to get only a hex number
    fp_str = remainder.replace(" ", "").replace("\n", "")
    # Python handles 1024 bit long integers without a problem
    #print("-------------------")
    #print(count)
    #print(fp_str)
    fp = int(fp_str, 16)
    # Merge all the one bits together
    union_bits |= fp
    count+=1
    #print(union_bits)
# This should start "0x1" if the first 3 bits are 0
    print(hex(union_bits))
    print(len(hex(union_bits)[2:-1]), "bytes")
    fingerprint.append(hex(union_bits)[0:240])
    #print("-------------------")
print(fingerprint)
#print(len(fingerprint),len(fingerprint[0]))
#print( list(fingerprint[0]) )


#---------drug fps-----(for bipartite)
from openbabel import *
import pybel
import pandas as pd
from rdkit import Chem

d=[]
Y=[]

#smiles=['CCCC','CCCN','CC(Cl)(Cl)Cl']
sol=pd.read_csv('Nuclear_Receptor.csv');
smiles=[]
for i in sol.SMILES:
    smiles.append(i)
mols=[pybel.readstring("smi",x)for x in smiles]
fps=[x.calcfp() for x in mols]
print("-------fps[0]-------------")
print(fps[0])

for i in fps:
    print(i.bits)



#------------Drug - Target Interaction---------------

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from networkx.algorithms import bipartite


g=nx.Graph()

#print(len(sol.DrugId))
#print(sol.DrugId)
#print(len(targetnames))
#print(targetnames)

g.add_nodes_from(targetnames,bipartite='Target')
g.add_nodes_from(sol.DrugId,bipartite='Drugs')


dataset=[]
# print(g.nodes(data=True))
#g.add_edges_from([(0, 'a')
# , (0, 'b'), (0, 'c'), (1, 'a'),(1, 'b'),(1,'c'),(2,'a')])
for i in range(len(X)):
    #print (customers[i],"--",res1[i])
    for j in sol.DrugId:
        #print("asas",j.bits)
        g.add_edge(targetnames[i],j)
print()
#print(dataset)
#g.add_edges_from((row['customers'], row['products']) for idx, row in df.iterrows()])
print(nx.bipartite.maximum_matching(g))
print(len(nx.bipartite.maximum_matching(g)))
print(nx.is_connected(g))
A, B = bipartite.sets(g)
pos = dict()
pos.update( (n, (1, i)) for i, n in enumerate(A) ) # put nodes from X at x=1
pos.update( (n, (2, i)) for i, n in enumerate(B) ) # put nodes from Y at x=2
nx.draw(g,pos ,with_labels=True)
plt.show()

#dataset=[[x,y.bits] for x in X for y in fps]
#print(dataset) 
#print(len(dataset))
'''
##dset=[[x,y] for x in X for y in fingerprint]
#for i in dset:
#    l=list(i[-1])
#    i[0]=i[0]+l
#    
#for i in dset:
#    #print (i)
#   i.pop()
#dataset=[]
#for i in dset:
#    dataset.append(i[0])
'''
dataset =[]
for i in range(len(X)):
    for j in range(len(fingerprint)):
        k = X[i].copy()
        k.extend(list(fingerprint[j]))
        dataset.append(k)

#print(dataset)
#for i in dataset:
#    print(len(i))
#    print(i)
#print(len(dataset))
#print(len(X[0]),len(fingerprint[0]))
#A=np.array(dataset)
#print(A.shape)

#--------------   + ve Samples   --------------------

import numpy as np
filename = "positivedrugsubset.fpt"

# Accumulate the union of all the fingerprint bits
union_bits = 0

# Extract only the hex value for each compound fingerprint
# This assumes the fingerprint file is small enough to fit into memory!
# (The parser was easier to write this way.)
count=1
positivefingerprint=[]
for block in open(filename).read().split("\n>"):
    union_bits = 0
    header, remainder = block.split("\n", 1)
    #print(header)
    if remainder.startswith("Possible superstructure"):
        _, remainder = remainder.split("\n", 1)

    # Remove spaces and newlines to get only a hex number
    fp_str = remainder.replace(" ", "").replace("\n", "")
    # Python handles 1024 bit long integers without a problem
    #print("-------------------")
    #print(count)
    #print(fp_str)
    fp = int(fp_str, 16)
    # Merge all the one bits together
    union_bits |= fp
    count+=1
    #print(union_bits)
# This should start "0x1" if the first 3 bits are 0
    print(hex(union_bits))
    print(len(hex(union_bits)[2:-1]), "bytes")
    positivefingerprint.append(hex(union_bits)[0:240])







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
# fps1 - Drug in old positive sample
# positivefingerprint - Drug in positive samples

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

positivedataset =[]
for i in range(len(X1)):
    k = X1[i].copy()
    k.extend(list(positivefingerprint[i]))
    positivedataset.append(k)

#print(len(dataset))
#print(len(positivedataset))



#------- Negative Samples------------
#--convert list to set



print(len(dataset))
print("--------------------------")
print(len(positivedataset))
#print(dataset[0])
#print(positivedataset[0])
negativedataset=[]
c=0
for i in dataset:
    c=0
    c=positivedataset.count(i)
    #print(c)
    if(c==0):
        negativedataset.append(i)
print(len(negativedataset))


#---------------------LASSO-------------------------------
from sklearn.metrics import accuracy_score as vals
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn.linear_model import Lasso
from sklearn.linear_model import LinearRegression
import random
import pickle
import numpy as np
x=[]
y=[]
for i in range(len(dataset)):
    y.append(random.randrange(0, 2)) 
Y=np.array(y)
X=np.array(dataset)
print(X.shape)
print(Y.shape)

X_train,X_test,Y_train,Y_test=train_test_split(X,Y, test_size=0.2)
#print(len(X_train))
#print(len(X_test))
#print(len(Y_train))
#print(len(Y_test))

lasso=Lasso()
lasso.fit(X_train,Y_train)