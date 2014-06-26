#!/usr/bin/python
from matplotlib import pyplot as plt
import numpy as np

def fillFigure(fig, data, target):
    ax = fig.add_axes([.1, .08, .85, .85])
    width = 0.3
    data3 = data[:, :3]
    data6 = data[:, 3:]

    left = np.r_[:nnodes]+1-width

    height = data3[:,0]
    bottom = np.zeros_like(height)
    plt.bar(left=left, height=height, bottom=bottom, width=width, axes=ax, color='r')

    bottom += height
    height = data3[:,1]
    plt.bar(left=left, height=height, bottom=bottom, width=width, axes=ax, color='g')

    bottom += height
    height = data3[:,2]
    plt.bar(left=left, height=height, bottom=bottom, width=width, axes=ax, color='b')

    ax.set_xlim((0.5, nnodes+.5))
    ax.set_xticks(np.r_[:nnodes]+1)
    ax.set_xticklabels(nodes)

    left = np.r_[:nnodes]+1

    height = data6[:,0]
    bottom = np.zeros_like(height)
    plt.bar(left=left, height=height, bottom=bottom, width=width, axes=ax, color='r')

    bottom += height
    height = data6[:,1]
    plt.bar(left=left, height=height, bottom=bottom, width=width, axes=ax, color='g')

    bottom += height
    height = data6[:,2]
    plt.bar(left=left, height=height, bottom=bottom, width=width, axes=ax, color='b')

    plt.xlabel('Number of nodes')
    plt.ylabel('Energy (Joules)')
    plt.grid(axis='y')
    plt.title('Energy consumption of %s parareal' % target)
    plt.legend({'Node', 'Network', 'Blowers'}, loc='upper left')

# CPU
dataCPU = np.loadtxt('energyCPU.dat')

nodes = np.int32(dataCPU[:,0])
data3 = dataCPU[:,1:4]
data6 = dataCPU[:,4:8]
nnodes = dataCPU.shape[0]

figCPU = plt.figure()
fillFigure(figCPU, dataCPU[:, 1:], 'CPU')


# GPU
dataGPU = np.loadtxt('energyGPU.dat')

nodes = np.int32(dataGPU[:,0])
data3 = dataGPU[:,1:4]
data6 = dataGPU[:,4:8]
nnodes = dataGPU.shape[0]

figGPU = plt.figure()
fillFigure(figGPU, dataGPU[:, 1:], 'GPU')

## Save figures
figCPU.savefig('08-EnergyCPU.eps')
figCPU.savefig('08-EnergyCPU.png')

figGPU.savefig('08-EnergyGPU.eps')
figGPU.savefig('08-EnergyGPU.png')


plt.show()
