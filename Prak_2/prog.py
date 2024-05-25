import matplotlib.pyplot as plt
import numpy as np

arrt=[]
for i in range(0,17):
    arrt.append([])

arrx=[]
for i in range(0,17):
    arrx.append([])

arry=[]
for i in range(0,17):
    arry.append([])

arrpx=[]
for i in range(0,17):
    arrpx.append([])

arrpy=[]
for i in range(0,17):
    arrpy.append([])

arru=[]
for i in range(0,17):
    arru.append([])

arr=[]
for i in range(0,17):
    arr.append([])


with open("t1.txt") as f:
    lines = f.read()
lines = lines.split("\n")

i=0
for line in lines:
    if len(line)== 0:
        break
    elif line=="END":
        i=i+1
    else:
        coordinates = line.split()
        arrt[i].append(float(coordinates[0]))
        arrx[i].append(float(coordinates[1]))
        arry[i].append(float(coordinates[2]))
        arrpx[i].append(float(coordinates[3]))
        arrpy[i].append(float(coordinates[4]))
        arru[i].append(float(coordinates[5]))
        arr[i].append(int(coordinates[6]))

fig,axes1=plt.subplots(2,1)
fig,axes2=plt.subplots(2,1)
fig,axes3=plt.subplots(1,1)



axes1[0].legend(loc='upper left')
axes1[0].set_xlabel("Ось t")
axes1[0].set_ylabel("Ось x")
axes1[1].legend(loc='upper left')
axes1[1].set_xlabel("Ось t")
axes1[1].set_ylabel("Ось y")

axes2[0].legend(loc='upper left')
axes2[0].set_xlabel("Ось t")
axes2[0].set_ylabel("Ось px")
axes2[1].legend(loc='upper left')
axes2[1].set_xlabel("Ось t")
axes2[1].set_ylabel("Ось py")

axes3.set_xlabel("Ось t")
axes3.set_ylabel("Ось u")

colors=['b','g','r','c','m','y','k','w']
i=14
axes1[0].plot(arrt[i],arrx[i],'grey',marker=',')
axes1[1].plot(arrt[i],arry[i],'grey',marker=',')
axes2[0].plot(arrt[i],arrpx[i],'grey',marker=',')
axes2[1].plot(arrt[i],arrpy[i],'grey',marker=',')
axes3.plot(arrt[i],arru[i],'grey',marker=',')

# lest=[2.893731,2.964275,9.957065]
# lesx=[-4.186838e+00,-4.391799e+00,-1.766136e-01]
# lesy=[-2.893731e+00,-2.893605e+00,4.099185e+00]
# lespx=[-2.559574e+01,-2.529316e+01,1.917492e+01]
# lespy=[-1.000000e+00,1.000000e+00,1.000000e+00]
# lesu=[-1.000000e+00,1.000000e+00,1.000000e+00]

lest=[5.853006,5.862724,19.993953]
lesx=[-1.712884e+01,-1.718574e+01,-5.007457e-02]
lesy=[-5.853006e+00,-5.853005e+00,8.278224e+00]
lespx=[-2.000278e+02,-1.998611e+02,1.570778e+02]
lespy=[-1.000000e+00,1.000000e+00,1.000000e+00]
lesu=[-1.000000e+00,1.000000e+00,1.000000e+00]


axes1[0].scatter(lest,lesx,color='green',marker='o')
axes1[1].scatter(lest,lesy,color='green',marker='o')
axes2[0].scatter(lest,lespx,color='green',marker='o')
axes2[1].scatter(lest,lespy,color='green',marker='o')
axes3.scatter(lest,lesu,color='green',marker='o')

plt.show()
