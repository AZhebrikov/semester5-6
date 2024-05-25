import numpy as np
import matplotlib.pyplot as plt
import math as mt

from numpy.linalg import inv




with open("rv1_2.txt","r") as f1:
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
    Y1_list.append(float(coordinates[0]))

L=len(X1_list)
i=0
for i in range(0,L):
    T1_list.append(int(i))     

#standart metod

X=np.zeros((2,2))
#MNK

x1_h=np.array(X1_list)
x2_h=np.array(Y1_list)

Dx1=np.array(x1_h)
Dx2=np.array(x2_h)

X[0,0]=x1_h.mean()
X[1,0]=x2_h.mean()

Dx1=Dx1-X[0,0]
Dx2=Dx2-X[1,0]

print(np.std(Dx1))
print(np.std(Dx2))

#differential metod

X1b=31.65630736842105 
X2b=42.06826210526315

x1_h=x1_h-X1b
x2_h=x2_h-X2b

Dxx1=np.array(x1_h)
Dxx2=np.array(x2_h)

X[0,1]=x1_h.mean()
X[1,1]=x2_h.mean()

Dxx1=Dxx1-X[0,1]
Dxx2=Dxx2-X[1,1]

print(np.std(Dxx1))
print(np.std(Dxx2))

X[0,1]=X1b+X[0,1]
X[1,1]=X2b+X[1,1]

print(X)
