import matplotlib.pyplot as plt
from numpy import append

with open("DEF_7_1.txt") as f:
    lines = f.read()
lines = lines.split("\n")

x=[]
y=[]

for line in lines:
    coordinates = line.split()
    if len(coordinates)==0:
        break
    x.append(float(coordinates[0]))
    y.append(float(coordinates[1]))


fig,ax=plt.subplots()

ax.plot(x,y,color='orange',marker=",")
ax.scatter(x,y,color='green',marker='.')
ax.scatter(x[0],y[0],color='blue',marker='o')

ax.set_xlabel("Ось t")
ax.set_ylabel("Ось x")
plt.show()