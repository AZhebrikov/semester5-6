import matplotlib.pyplot as plt
from numpy import append

with open("DEF_8.txt") as f:
    lines = f.read()
lines = lines.split("\n")

x=[[],[]]
y=[[],[]]

i=0

for line in lines:
    if line=="END":
        i=i+1
        continue
    coordinates = line.split()
    if len(coordinates)==0:
        break
    x[i].append(float(coordinates[0]))
    y[i].append(float(coordinates[1]))


fig,ax=plt.subplots()

ax.scatter(x[0],y[0],color='blue',marker='.')
ax.scatter(x[1],y[1],color='orange',marker='.')

ax.scatter(x[0][0],y[0][0],color='blue',marker='o')
ax.scatter(x[1][0],y[1][0],color='orange',marker='o')

ax.set_xlabel("Ось x")
ax.set_ylabel("Ось y")
plt.show()