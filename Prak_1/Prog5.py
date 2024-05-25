import matplotlib.pyplot as plt
import numpy as np

arrt=[]
for i in range(0,4):
    arrt.append([])

arrx=[]
for i in range(0,4):
    arrx.append([])

arry=[]
for i in range(0,4):
    arry.append([])

arr=[]
for k in range(0,4):
    arr.append([[],[],[]])
    

with open("DEF_5.txt") as f:
    lines = f.read()
lines = lines.split("\n")


i=-1
flag=0
for line in lines:
    if len(line)==0:
        break
    elif line=="END":
        i=i+1
        flag=7
    elif line=="---":
        flag=0
    elif flag==7:
        coordinates = line.split()
        arrt[i].append(float(coordinates[0]))
        arrx[i].append(float(coordinates[1]))
        arry[i].append(float(coordinates[2]))
    else:
        coordinates = line.split()
        try:
            arr[i][0].append(float(coordinates[0]))
        except ValueError:
            arr[i][0].append(str(coordinates[0]))
        try:
            arr[i][1].append(float(coordinates[1]))
        except ValueError:
            arr[i][1].append(str(coordinates[1]))
        try:
            arr[i][2].append(float(coordinates[2]))
        except ValueError:
            arr[i][2].append(str(coordinates[2]))
        flag=flag+1


fig1,axes1=plt.subplots()
axes1.plot(arrx[0],arry[0],color='green',marker=',')

fig2,axes2=plt.subplots()
axes2.plot(arrx[1],arry[1],color='green',marker=',')

fig3,axes3=plt.subplots()
axes3.plot(arrx[2],arry[2],color='green',marker=',')

fig4,axes4=plt.subplots()
axes4.plot(arrx[3],arry[3],color='green',marker=',')

plt.show()
