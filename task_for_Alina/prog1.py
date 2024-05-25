import numpy as np
import matplotlib.pyplot as plt
import math as mt

from numpy.linalg import inv

with open("t13.txt") as f1:
    lines = f1.read()
lines = lines.split("\n")

x=[[],[],[],[],[],[],[],[],[],[],[],[]]
y=[[],[],[],[],[],[],[],[],[],[],[],[]]

i=0
j=0
for line in lines:
    if i==0:
        coordinates=line.split()
        for s in coordinates:
             x[j].append(float(s))
        i=i+1
    elif i==1:
        coordinates=line.split()
        for s in coordinates:
             y[j].append(float(s))
        i=0
        j=j+1
    elif len(coordinates)==0:
        break

z=[-1,-3./8,-3./8 + 0.0001,1]
Z=[0,0,1,1]

# Z=[]
# L=len(x[0])
# for i in range(L):
#     if x[0][i]<0.:
#         Z.append(float(y[0][i]))
#     elif x[0][i]>1.:
#         Z.append(float(y[0][i]-1.))
#     else: 
#         Z.append(float(y[0][i]-x[0][i]))
        

# fig,ax=plt.subplots(5,2)

fig,ax=plt.subplots()
 
# print(len(x[4]))
# print(len(y[4]))
# print(len(Z))


# ax.scatter(x[1],y[1],color='blue',marker='.')
# ax.scatter(x[2],y[2],color='red',marker='.')

# ax.plot(x[4],Z,color='gray',marker='.',label='Численное решение')

ax.plot(x[0],y[0],color='gray',marker='.',label='Численное решение')
ax.plot(z,Z,color='blue',marker='.',label='Истинное решение')

# ax[0,1].plot(x[2],y[2],color='gray',marker='.',label='Численное решение')
# ax[0,1].plot(z,Z,color='blue',marker='.',label='Истинное решение')


# ax[0,0].plot(x[1],y[1],color='gray',marker='.',label='0.100000')
# ax[0,0].plot(z,Z,color='blue',marker='.')
# ax[0,0].grid(True)
# ax[0,0].legend(loc='upper left')

# ax[0,1].plot(x[2],y[2],color='gray',marker='.',label='0.076923')
# ax[0,1].plot(z,Z,color='blue',marker='.')
# ax[0,1].grid(True)
# ax[0,1].legend(loc='upper left')

# ax[1,0].plot(x[3],y[3],color='gray',marker='.',label='0.059172')
# ax[1,0].plot(z,Z,color='blue',marker='.')
# ax[1,0].grid(True)
# ax[1,0].legend(loc='upper left')

# ax[1,1].plot(x[4],y[4],color='gray',marker='.',label='0.045517')
# ax[1,1].plot(z,Z,color='blue',marker='.')
# ax[1,1].grid(True)
# ax[1,1].legend(loc='upper left')

# ax[2,0].plot(x[5],y[5],color='gray',marker='.',label='0.035013')
# ax[2,0].plot(z,Z,color='blue',marker='.')
# ax[2,0].grid(True)
# ax[2,0].legend(loc='upper left')

# ax[2,1].plot(x[6],y[6],color='gray',marker='.',label='0.026933')
# ax[2,1].plot(z,Z,color='blue',marker='.')
# ax[2,1].grid(True)
# ax[2,1].legend(loc='upper left')

# ax[3,0].plot(x[7],y[7],color='gray',marker='.',label='0.015937')
# ax[3,0].plot(z,Z,color='blue',marker='.')
# ax[3,0].grid(True)
# ax[3,0].legend(loc='upper left')

# ax[3,1].plot(x[8],y[8],color='gray',marker='.',label='0.012259')
# ax[3,1].plot(z,Z,color='blue',marker='.')
# ax[3,1].grid(True)
# ax[3,1].legend(loc='upper left')

# ax[4,0].plot(x[9],y[9],color='gray',marker='.',label='0.009430')
# ax[4,0].plot(z,Z,color='blue',marker='.')
# ax[4,0].grid(True)
# ax[4,0].legend(loc='upper left')

# ax[4,1].plot(x[9],y[9],color='gray',marker='.',label='0.006223')
# ax[4,1].plot(z,Z,color='blue',marker='.')
# ax[4,1].grid(True)
# ax[4,1].legend(loc='upper left')


# ax.plot(x[0],Z,color='gray',marker='.')

# ax.set_xlabel("Ось X")
# ax.set_ylabel("Ось V-U")

# fig1,ax1=plt.subplots()
 
# ax1.scatter(x[3],y[3],color='gray',marker='.')
# ax1.scatter(x[4],y[4],color='blue',marker='.')
# ax1.scatter(x[5],y[5],color='red',marker='.')

# ax1.set_xlabel("Ось X1")
# ax1.set_ylabel("Ось X2")

# fig2,ax2=plt.subplots()
 
# ax2.scatter(x[6],y[6],color='gray',marker='.')
# ax2.scatter(x[7],y[7],color='blue',marker='.')
# ax2.scatter(x[s8],y[8],color='red',marker='.')

# ax2.set_xlabel("Ось X1")
# ax2.set_ylabel("Ось X2")
# plt.ylim(-0.5, 1.25)
# plt.xlim(-0.6,0.6)


plt.legend()
plt.grid(True)
plt.show()

