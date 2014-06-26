#!/usr/bin/python
from matplotlib import pyplot as plt
import numpy as np

data = np.loadtxt('energy.dat')

nodes = np.int32(data[:,0])
dataCPU = data[:,1:4]
dataGPU = data[:,4:8]
nnodes = data.shape[0]

# First plot: CPU
figCPU = plt.figure()
axCPU = figCPU.add_axes([.1, .08, .85, .85])
width = 0.7

left = np.r_[:nnodes]+1-(width/2)

height = dataCPU[:,0]
bottom = np.zeros_like(height)
plt.bar(left=left, height=height, bottom=bottom, width=width, axes=axCPU, color='r')

bottom += height
height = dataCPU[:,1]
plt.bar(left=left, height=height, bottom=bottom, width=width, axes=axCPU, color='g')

bottom += height
height = dataCPU[:,2]
plt.bar(left=left, height=height, bottom=bottom, width=width, axes=axCPU, color='b')

axCPU.set_xlim((0.5, nnodes+.5))
axCPU.set_xticks(np.r_[:nnodes]+1)
axCPU.set_xticklabels(nodes)

plt.xlabel('Number of nodes')
plt.ylabel('Energy (Joules)')
plt.grid(axis='y')
plt.title('Energy consumption of CPU parareal')
plt.legend({'Node', 'Network', 'Blowers'}, loc='upper left')

# First plot: GPU
figGPU = plt.figure()
axGPU = figGPU.add_axes([.1, .08, .85, .85])
width = 0.7

left = np.r_[:nnodes]+1-(width/2)

height = dataGPU[:,0]
bottom = np.zeros_like(height)
plt.bar(left=left, height=height, bottom=bottom, width=width, axes=axGPU, color='r')

bottom += height
height = dataGPU[:,1]
plt.bar(left=left, height=height, bottom=bottom, width=width, axes=axGPU, color='g')

bottom += height
height = dataGPU[:,2]
plt.bar(left=left, height=height, bottom=bottom, width=width, axes=axGPU, color='b')

axGPU.set_xlim((0.5, nnodes+.5))
axGPU.set_xticks(np.r_[:nnodes]+1)
axGPU.set_xticklabels(nodes)

plt.xlabel('Number of nodes')
plt.ylabel('Energy (Joules)')
plt.grid(axis='y')
plt.title('Energy consumption of GPU parareal')
plt.legend({'Node', 'Network', 'Blowers'}, loc='upper left')

# Set ylim
ymax = max(axCPU.get_ylim()[1], axGPU.get_ylim()[1])
axCPU.set_ylim((0, ymax))
axGPU.set_ylim((0, ymax))

# Save figures
figCPU.savefig('08-EnergyCPU.eps')
figCPU.savefig('08-EnergyCPU.png')

figGPU.savefig('08-EnergyGPU.eps')
figGPU.savefig('08-EnergyGPU.png')


plt.show()
