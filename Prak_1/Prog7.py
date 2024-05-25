import matplotlib.pyplot as plt
from numpy import append

with open("DEF_8.txt") as f:
    lines = f.read()
lines = lines.split("\n")

x=[]
y=[]

for line in lines:
    if line=="END":
        continue
    coordinates = line.split()
    if len(coordinates)==0:
        break
    x.append(float(coordinates[0]))
    y.append(float(coordinates[1]))


fig,ax=plt.subplots()

ax.scatter(x,y,color='blue',marker=",")

ax.scatter(x[0],y[0],color='blue',marker='o')

ax.set_xlabel("Ось x")
ax.set_ylabel("Ось y")
plt.show()