import numpy as np
import matplotlib.pyplot as plt
import math as mt

from numpy.linalg import inv


X1b=31.65630736842105 
X2b=42.06826210526315

with open("rv1_1.txt","r") as f1:
    lines = f1.read()
lines = lines.split("\n")

X1_list=[]
Y1_list=[]
T1_list=[]

for line in lines:
    coordinates = line.split()
    if len(coordinates)==0:
        break
    X1_list.append(float(coordinates[1])-X1b)
    Y1_list.append(float(coordinates[0])-X2b)

L=len(X1_list)
i=0
for i in range(0,L):
    T1_list.append(int(i))     

#differential metod

A=np.zeros((3,2))
#MNK
e=np.ones(L)

x1_h=np.array(X1_list)
x1_hh=x1_h-x1_h.mean()*e

x2_h=np.array(Y1_list)
x2_hh=x2_h-x2_h.mean()*e

A[0,0]=(-1)*(x1_hh.dot(x2_hh.T))/((x1_hh**2).sum())
A[0,1]=-x1_h.mean()*A[0,0]-x2_h.mean()
#RNK
C=np.empty((L,2))
for i in range(0,L):
    C[i]=[x1_hh[i],x2_hh[i]]

U,S,V=np.linalg.svd(C)

A[1,0]=V[1,0]/V[1,1]
A[1,1]=-((x1_h*A[1,0]+x2_h).mean())

#MNK
maxim=1000000
for i in range(0,L):
    for j in range(0,L):
        if x1_h[i]>x1_h[j] or x1_h[i]<x1_h[j]:
            P=np.array([[x1_h[i],1],[x1_h[j],1]])
            G=np.array([-x2_h[i],-x2_h[j]])
            P = inv(P)
            N=P.dot(G)
            Delta=0
            for k in range(0,L):
                Delta=Delta+abs(N[0]*x1_h[k]+N[1]+x2_h[k])
            if Delta < maxim:
                A[2]=N
                maxim=Delta

x1_1 =np.array([x1_h[0],x1_h[L-1]])
x2_1 = -A[0,1]-A[0,0]*x1_1
x1_2 =np.array([x1_h[0],x1_h[L-1]])
x2_2 = -A[1,1]-A[1,0]*x1_1
x1_3 =np.array([x1_h[0],x1_h[L-1]])
x2_3 = -A[2,1]-A[2,0]*x1_1

str1='MНК'
str2='РНК'
str3='MНМ'

fig2,ax2=plt.subplots()
 
ax2.plot(X1_list,Y1_list,color='gray',marker='.',label='Результаты')

ax2.plot(x1_1,x2_1,color='green',marker='.',label=str1)
ax2.plot(x1_2,x2_2,color='red',marker='.',label=str2)
ax2.plot(x1_3,x2_3,color='blue',marker='.',label=str3)

ax2.set_xlabel("Ось x")
ax2.set_ylabel("Ось y")
plt.grid(True)
plt.legend()
plt.show()

print("Матрица коэффициентов прямой для соответствующих методов.") 
print(A)

print("Выражение наклона в радианах, градусах, минутах и секундах.") 
for s in range(0,3):
    print( mt.atan(A[s,0]) , mt.degrees(mt.atan(A[s,0])), 60*mt.degrees(mt.atan(A[s,0])) , 3600*mt.degrees(mt.atan(A[s,0])) )

Dx1=abs(x1_h*A[0,0]+x2_h+A[0,1])/mt.sqrt(A[0,0]**2+1)
Dx2=abs(x1_h*A[1,0]+x2_h+A[1,1])/mt.sqrt(A[1,0]**2+1)
Dx3=abs(x1_h*A[2,0]+x2_h+A[2,1])/mt.sqrt(A[2,0]**2+1)
print('Точность оценки')
print('МНК: ',np.std(Dx1))
print('РНК: ',np.std(Dx2))
print('МНМ: ',np.std(Dx3))
