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
# fname=os.getcwd()+"\\DrugBank_Approved-PSSM\\spP45059.pssm"
#for fname in glob.glob(os.getcwd()+"\\DrugBank_Approved-PSSM\\*.pssm"):
for fname in glob.glob(os.getcwd()+"\\Nuclear receptor-PSSM\\*.pssm"):
    #print(fname)
    # print(fname[-13:-5])
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
    X.append(str(result))
    #X.append(result)
# print(len(X))
# print(len(result),len(X[-1]))
# print(result)
# print(X)
for i in range(5):
    print("------",i,"------")
    print(X[i])

#     print("-------------------- ",i," ----------------------")




#------------Drug---------------
#pybel smiles to fp2
#https://open-babel.readthedocs.io/en/latest/UseTheLibrary/Python_Pybel.html

#drug dataset
#https://towardsdatascience.com/how-to-use-machine-learning-for-drug-discovery-1ccb5fdf81ad

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
#print(fps)



#from rdkit import DataStructs
#ms = [Chem.MolFromSmiles('[H][C@@]1(CC[C@@]2([H])[C@]3([H])CC=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)CCCC(C)C'), Chem.MolFromSmiles('CCO'),Chem.MolFromSmiles('COC')]
#fps = [Chem.RDKFingerprint(x) for x in ms]
#print(fps[0])
#for i in fps[0]:
#    print(i , end=" ")



#------------Drug - Target Interaction---------------

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from networkx.algorithms import bipartite


g=nx.Graph()

g.add_nodes_from(X,bipartite='Target')
g.add_nodes_from(fps,bipartite='Drugs')

# print(customers[0:10])
# g1.add_nodes_from(customers[:10],bipartite='customers')
# g1.add_nodes_from(products[:5],bipartite='products')

dataset=[]
# print(g.nodes(data=True))
#g.add_edges_from([(0, 'a')
# , (0, 'b'), (0, 'c'), (1, 'a'),(1, 'b'),(1,'c'),(2,'a')])
for i in range(len(X)):
    #print (customers[i],"--",res1[i])
    for j in fps:
        #print("asas",j.bits)
        g.add_edge(X[i],j)
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
nx.draw(g,pos ,with_labels=False)
plt.show()

dataset=[[x,y.bits] for x in X for y in fps]
#print(dataset) 
#print(len(dataset))


#--------------   + ve Samples   --------------------


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
    X1.append(str(result))
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

#------- Negative Samples------------
#--convert list to set



print(len(dataset))
print("--------------------------")
print(len(positive_dataset))
negative_dataset=[]
c=0
for i in dataset:
    c=0
    c=positive_dataset.count(i)
    #print(c)
    if(c==0):
        negative_dataset.append(i)
print(len(negative_dataset))


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
#print(len(X),len(Y))
#print(X.shape,Y.shape)

X_train,X_test,Y_train,Y_test=train_test_split(X,Y, test_size=0.2)
print(X_train.shape)
print(X_test.shape)
print(Y_train.shape)
print(Y_test.shape)
lasso=Lasso()
#lasso.fit(X_train,Y_train)
#train_score=lasso.score(X_train,Y_train)
print(X[0][0])
print(X[0][1])
for i in X[0][1]:
    X[0][0].append(i)
print(X[0][0])
