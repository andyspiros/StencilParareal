#!/usr/bin/python
import numpy as np
from matplotlib import pyplot as plt

data = np.loadtxt('runtime.dat')

nodes = np.int32(data[:,0])
threads = np.int32(2.**np.r_[:data.shape[1]-1])
speedup = (np.repeat(data[:,1:2], 4, axis=1) / data[:,1:]).transpose()

fig = plt.figure()
ax = fig.add_axes([.1, .08, .85, .85])

lines = plt.plot(threads, speedup, marker='^')
for i,l in enumerate(lines):
    l.set_label("%d nodes" % nodes[i])

ymin=speedup.min()*3./4.
ymax=speedup.max()*4./3.

ax.set_xscale('log')
ax.set_yscale('log')
#ax.set_xlim(xmin, xmax)
ax.set_ylim(1., 6.)
ax.set_xticks(threads)
ax.set_xticklabels(threads)
ax.set_yticks((1, 2, 3, 4, 5, 6))
ax.set_yticklabels((1, 2, 3, 4, 5, 6))
plt.xlabel('Number of threads')
plt.ylabel('Speedup')
plt.grid(True)
plt.title('Effect of threads on parareal')
plt.legend(loc='lower right')

fig.savefig('02-SpeedupThreads.eps')
fig.savefig('02-SpeedupThreads.png')

plt.show()

