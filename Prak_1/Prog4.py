import matplotlib.pyplot as plt
import numpy as np

arrx=[]
for i in range(0,6):
    arrx.append([])

arry=[]
for i in range(0,6):
    arry.append([])

arr=[[],[]]

with open("DEF_4.txt") as f:
    lines = f.read()
lines = lines.split("\n")

i=-1
flag=0
for line in lines:
    if len(line)== 0:
        break
    elif line=="END":
        i=i+1
        flag=0
    elif line=="ENDEND":
        1
    elif line=="---":
        flag=1
    elif flag==0:
        coordinates = line.split()
        arrx[i].append(float(coordinates[0]))
        arry[i].append(float(coordinates[1]))
    else:
        coordinates = line.split()
        arr[0].append(float(coordinates[0]))
        arr[1].append(float(coordinates[1]))


fig,axes=plt.subplots(2,1)

axes[0].plot(arrx[0],arry[0],color='orange',marker='.',
label=('h= '+str(arr[0][0])+',y(10)-y*(10)= '+str(arr[1][0]))
)
axes[0].plot(arrx[1],arry[1],color='green',marker=',',
label=('h= '+str(arr[0][1])+',y(10)-y*(10)= '+str(arr[1][1]))
)
axes[0].plot(arrx[2],arry[2],color='blue',marker=',',
label=('h= '+str(arr[0][2])+',y(10)-y*(10)= '+str(arr[1][2]))
)
axes[0].legend(loc='upper left')
axes[0].set_xlabel("Ось t")
axes[0].set_ylabel("Ось F_1")

axes[1].plot(arrx[3],arry[3],color='orange',marker='.',
label=('h= '+str(arr[0][3])+',y(10)-y*(10)= '+str(arr[1][3]))
)
axes[1].plot(arrx[4],arry[4],color='green',marker=',',
label=('h= '+str(arr[0][4])+',y(10)-y*(10)= '+str(arr[1][4]))
)
axes[1].plot(arrx[5],arry[5],color='blue',marker=',',
label=('h= '+str(arr[0][5])+',y(10)-y*(10)= '+str(arr[1][5]))
)
axes[1].legend(loc='upper left')
axes[1].set_xlabel("Ось t")
axes[1].set_ylabel("Ось F_2")

plt.show()
