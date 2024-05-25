import matplotlib.pyplot as plt
from numpy import append

with open("DEF_3.txt") as f:
    lines = f.read()
lines = lines.split("\n")

x=list()
y=list()

for line in lines:
    coordinates = line.split()
    if len(coordinates)==0:
        break
    x.append(float(coordinates[0]))
    y.append(float(coordinates[1]))


fig,ax=plt.subplots(figsize=(5,1))

ax.plot(x,y,color='blue',marker='.')
ax.set_xlabel("Ось h")
ax.set_ylabel("Ось R")
plt.show()