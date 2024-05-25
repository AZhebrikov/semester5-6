import numpy as np
from numpy.linalg import inv

import matplotlib.pyplot as plt




with open("ROVER2.GK","r") as f1:
    lines = f1.read()
lines = lines.split("\n")

X1_list=[]
Y1_list=[]
T1_list=[]

for line in lines:
    coordinates = line.split()
    if len(coordinates)==0:
        break
    X1_list.append(float(coordinates[1]))
    Y1_list.append(float(coordinates[2]))

L=len(X1_list)
i=0
for i in range(0,L):
    T1_list.append(int(i))     


     
     
fig1,ax1=plt.subplots(2)
 
ax1[0].scatter(T1_list,X1_list,color='gray',marker='.')
ax1[1].scatter(T1_list,Y1_list,color='blue',marker='.')

ax1[0].set_xlabel("Ось t")
ax1[0].set_ylabel("Ось x")

ax1[1].set_xlabel("Ось t")
ax1[1].set_ylabel("Ось y")

ax1[0].grid(True)
ax1[1].grid(True)

plt.show()


fig2,ax2=plt.subplots()
 
ax2.scatter(X1_list,Y1_list,color='gray',marker='.',label='Результаты')

ax2.set_xlabel("Ось x")
ax2.set_ylabel("Ось y")
plt.grid(True)
plt.legend()
plt.show()
