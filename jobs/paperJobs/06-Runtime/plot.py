#!/usr/bin/python
from matplotlib import pyplot as plt
import numpy as np

print "Loading data"
data = np.loadtxt('runtime.dat')
print "Data loaded"

nodes = np.int32(data[:,0])
parallelCPU = data[:,1]
serialCPU = np.mean(data[:,2]) * np.ones_like(parallelCPU)
parallelGPU = data[:,3]
serialGPU = np.mean(data[:,4]) * np.ones_like(parallelGPU)

fig = plt.figure()
ax = fig.add_axes([.1, .08, .85, .85])

plt.plot(nodes, parallelCPU, axes=ax, linestyle='--', linewidth=1, marker='^', color='b', label='Parareal CPU')
plt.plot(nodes, serialCPU  , axes=ax, linestyle='-' , linewidth=1,             color='b', label='Time-serial CPU')
plt.plot(nodes, parallelGPU, axes=ax, linestyle='--', linewidth=1, marker='^', color='r', label='Parareal GPU')
plt.plot(nodes, serialGPU  , axes=ax, linestyle='-' , linewidth=1,             color='r', label='Time-serial GPU')

xmin=nodes.min()*3./4.
xmax=nodes.max()*4./3.
ymin=data[:,1:].min()*3./4.
ymax=data[:,1:].max()*4./3.

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_xticks(nodes)
ax.set_xticklabels(nodes)
#ax.set_yticks((50, 100, 200, 500, 1000))
#ax.set_yticklabels((50, 100, 200, 500, 1000))
plt.xlabel('Number of nodes')
plt.ylabel('Runtime (seconds)')
plt.grid(True)
plt.title('Runtime of parareal and time-serial solver')
plt.legend()

fig.savefig('06-Runtime.eps')
fig.savefig('06-Runtime.png')

plt.show()

