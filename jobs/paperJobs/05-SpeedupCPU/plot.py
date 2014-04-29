#!/usr/bin/python
from matplotlib import pyplot as plt
import numpy as np

data = np.loadtxt('speedup.dat')

nodes = np.int32(data[:,0])

fig = plt.figure()
ax = fig.add_axes([.1, .08, .85, .85])

plt.plot(nodes, data[:,1], axes=ax, linestyle='-' , linewidth=1, marker='^', color='b', label='Speedup')
plt.plot(nodes, data[:,2], axes=ax, linestyle='--', linewidth=1,             color='b', label='Max speedup')

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
ax.set_yticks((1, 2, 4, 8, 16, 32))
ax.set_yticklabels((1, 2, 4, 8, 16, 32))
plt.xlabel('Number of nodes')
plt.ylabel('Speedup')
plt.grid(True)
plt.title('Comparison of parareal VS serial solver')
plt.legend()

fig.savefig('05-SpeedupCPU.eps')
fig.savefig('05-SpeedupCPU.png')

plt.show()

