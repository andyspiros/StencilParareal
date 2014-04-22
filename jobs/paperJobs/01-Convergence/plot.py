#!/usr/bin/python

import numpy as np
from matplotlib import pyplot as plt

np.loadtxt("result.dat")
data = np.loadtxt("result.dat")
nodes = data[:,0]
nufreq = data[:,1]
data = data[:,2:]
k = np.r_[1:9]

f = plt.figure()
ax = f.add_axes([.1, .08, .85, .85])

plt.plot(k, data[0,:], axes=ax, linestyle='--', linewidth=2, marker='^', color='r', label='8 nodes, freq=0')
plt.plot(k, data[1,:], axes=ax, linestyle=':' , linewidth=2, marker='s', color='r', label='8 nodes, freq=100')
plt.plot(k, data[2,:], axes=ax, linestyle='--', linewidth=2, marker='^', color='g', label='32 nodes, freq=0')
plt.plot(k, data[3,:], axes=ax, linestyle=':' , linewidth=2, marker='s', color='g', label='32 nodes, freq=100')
plt.plot(k, data[4,:], axes=ax, linestyle='--', linewidth=2, marker='^', color='b', label='128 nodes, freq=0')
plt.plot(k, data[5,:], axes=ax, linestyle=':' , linewidth=2, marker='s', color='b', label='128 nodes, freq=100')

ax.set_yscale('log')
plt.xlabel('Number of iterations k')
plt.ylabel('Relative error')
plt.title('Commparison of parareal VS serial solver')
plt.legend(True)

f.savefig('01-Comparison.eps')
f.savefig('01-Comparison.png')
