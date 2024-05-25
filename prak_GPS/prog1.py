import numpy as np
import matplotlib.pyplot as plt
import math as mt

from numpy.linalg import inv

# with open("base.txt") as f1:
#     lines = f1.read()
# lines = lines.split("\n")

# x=[]
# y=[]

# for line in lines:
#     coordinates = line.split()
#     if len(coordinates)==0:
#             break
#     x.append(float(coordinates[1]))
#     y.append(float(coordinates[0]))


# L=len(x)
# xb1_h=np.array(x)
# xb2_h=np.array(y)

# print(xb1_h.mean(),xb2_h.mean())
# print(np.median(xb1_h),np.median(xb2_h))

# xb1=np.median(xb1_h)
# xb2=np.median(xb2_h)

xb1=0.
xb2=0.

with open("ROVER2.GK") as f2:
    lines = f2.read()
lines = lines.split("\n")

X=[]
Y=[]

for line in lines:
    coordinates = line.split()
    if len(coordinates)==0:
            break
    X.append(float(coordinates[1]))
    Y.append(float(coordinates[2]))

print(X)
print(Y)

L=len(X)

xc1=np.array(X)
xc2=np.array(Y)
Xc1=np.median(X)
Xc2=np.median(Y)
print(Xc1,Xc2)

r1=[]
r2=[]
for j in range(0,L):
    r1.append(float(pow((X[j]-xb1)*(X[j]-xb1)+(Y[j]-xb2)*(Y[j]-xb2),1/2)))
    r2.append(float(pow((X[j]-Xc1)*(X[j]-Xc1)+(Y[j]-Xc2)*(Y[j]-Xc2),1/2)))
    
r1_h=np.array(r1)
r2_h=np.array(r2)

R_b1=r1_h.mean()
R_b2=np.median(r1_h)
R_c1=r2_h.mean()
R_c2=np.median(r2_h)

print(R_b1)
print(R_b2)
print(R_c1)
print(R_c2)

fig,ax=plt.subplots()
Drawing_uncolored_circle1 = plt.Circle((xb1,xb2), R_b1,color='green',label='МНК', fill = False)
Drawing_uncolored_circle2 = plt.Circle((xb1,xb2), R_b2,color='lightgreen',label='МНМ', fill = False)
Drawing_uncolored_circle3 = plt.Circle((Xc1,Xc2), R_c1,color='blue',label='МНК', fill = False)
Drawing_uncolored_circle4 = plt.Circle((Xc1,Xc2), R_c2,color='lightblue',label='МНМ', fill = False)
 
ax.add_artist(Drawing_uncolored_circle1)
ax.add_artist(Drawing_uncolored_circle2)
ax.add_artist(Drawing_uncolored_circle3)
ax.add_artist(Drawing_uncolored_circle4)

# ax.scatter(x,y,color='gray',marker='.',label='данные Базы')
# ax.scatter(np.median(xb1_h),np.median(xb2_h),color='darkgreen',marker='o',label='Коорд. Базы (МНК)')
# ax.scatter(np.median(xc1),np.median(xc2),color='darkblue',marker='o',label='Коорд. центра (МНМ)')

ax.scatter(X,Y,color='salmon',marker='.',label='данные ровера')

ax.set_xlabel("Ось X1")
ax.set_ylabel("Ось X2")

plt.legend()
plt.grid(True)
plt.show()
