import sys, os
from os.path import join as pjoin

import matplotlib
matplotlib.use('Agg')

import numpy as np
from matplotlib import pyplot as plt
plt.rc('text', usetex=True)


with open('result.dat') as fs:
    lines = fs.readlines()[2:]

executables = set([l.split()[0] for l in lines])

for exe in executables:
    data = np.array([l.split()[1:] for l in lines if l.split()[0] == exe],
                    dtype=np.double)

    # First plot: Runtime of parareal
    f1 = plt.figure()
    plt.hold(True)
    ks = np.unique(data[:,1])
    for k in ks:
        idx = data[:,1] == k
        nodes = data[idx, 0]
        runtime = data[idx, 6]
        label = 'k = {0:.0f}'.format(k)
        plt.plot(nodes, runtime, label=label, marker='^')
    serial = np.mean(data[:,5]).repeat(len(nodes))
    plt.plot(nodes, serial, label="Serial", color='k', linewidth=2)
    a = f1.axes[0]
    a.set_title("Runtime of {0}".format(exe))
    a.set_xscale('log')
    a.set_yscale('log')
    a.set_xlim((nodes.min()*3./4., nodes.max()*4./3.))
    a.set_xticks(nodes)
    a.set_xticklabels(['{0:.0f}'.format(t) for t in nodes])
    a.set_xlabel("Nodes")
    a.set_ylabel("Time (msec)")
    a.grid(True, which='major')
    a.legend(loc="lower left")
    f1.savefig("runtime_{0}.pdf".format(exe[-3:]), format="pdf", bbox_inches="tight")

    # Second plot: Speedup of parareal
    f2 = plt.figure()
    plt.hold(True)
    ks = np.unique(data[:,1])
    for k in ks:
        idx = data[:,1] == k
        nodes = data[idx, 0]
        speedup = data[idx, 2]
        maxspeedup = data[idx, 3]
        label = 'k = {0:.0f}'.format(k)
        line = plt.plot(nodes, speedup, label=label, marker='^')[0]
        plt.plot(nodes, maxspeedup, linestyle='--', color=line.get_color())
    a = f2.axes[0]
    a.set_title("Speedup of {0}".format(exe))
    a.set_xscale('log')
    a.set_yscale('log')
    a.set_xlim((nodes.min()*3./4., nodes.max()*4./3.))
    a.set_xticks(nodes)
    a.set_xticklabels(['{0:.0f}'.format(t) for t in nodes])
    a.set_xlabel("Nodes")
    a.set_ylabel("Speedup")
    a.grid(True, which='major')
    a.legend(loc="lower right")
    f2.savefig("speedup_{0}.pdf".format(exe[-3:]), format="pdf", bbox_inches="tight")

    # Third plot: error
    f3 = plt.figure()
    plt.hold(True)
    ks = np.unique(data[:,1])
    for k in ks:
        idx = data[:,1] == k
        nodes = data[idx, 0]
        error = data[idx, 4]
        label = 'k = {0:.0f}'.format(k)
        plt.plot(nodes, error, label=label, marker='^')[0]
    a = f3.axes[0]
    a.set_title("Relative error of {0}".format(exe))
    a.set_xscale('log')
    a.set_yscale('log')
    a.set_xlim((nodes.min()*3./4., nodes.max()*4./3.))
    a.set_xticks(nodes)
    a.set_xticklabels(['{0:.0f}'.format(t) for t in nodes])
    a.set_xlabel("Nodes")
    a.set_ylabel("Error")
    a.grid(True, which='major')
    a.legend(loc="lower right")
    f3.savefig("error_{0}.pdf".format(exe[-3:]), format="pdf", bbox_inches="tight")

    # Fourth plot: energy increase
    f4 = plt.figure()
    plt.hold(True)
    ks = np.unique(data[:,1])
    for k in ks:
        idx = data[:,1] == k
        nodes = data[idx, 0]
        energyInc = data[idx, 8] / data[idx, 7]
        label = 'k = {0:.0f}'.format(k)
        plt.plot(nodes, energyInc, label=label, marker='^')[0]
    a = f4.axes[0]
    a.set_title("Increase in energy consumption of {0}".format(exe))
    a.set_xscale('log')
    a.set_yscale('log')
    a.set_xlim((nodes.min()*3./4., nodes.max()*4./3.))
    a.set_xticks(nodes)
    a.set_xticklabels(['{0:.0f}'.format(t) for t in nodes])
    a.set_xlabel("Nodes")
    a.set_ylabel(r'$\frac{E(\mathrm{parareal})}{\mathrm{E(serial)}}$', fontsize=18)
    a.grid(True, which='major')
    a.legend(loc="lower right")
    f4.savefig("energyNodes_{0}.pdf".format(exe[-3:]), format="pdf", bbox_inches="tight")

    # Fifth plot: energy increase VS speedup
    f5 = plt.figure()
    plt.hold(True)
    ks = np.unique(data[:,1])
    for k in ks:
        idx = data[:,1] == k
        nodes = data[idx, 0]
        energyInc = data[idx, 8] / data[idx, 7]
        speedup = data[idx, 2]
        label = 'k = {0:.0f}'.format(k)
        plt.plot(speedup, energyInc, label=label, marker='^')[0]
    a = f5.axes[0]
    a.set_title("Increase in energy consumption of {0}".format(exe))
    a.set_yscale('log')
    a.set_xlabel("Speedup")
    a.set_ylabel(r'$\frac{E(\mathrm{parareal})}{\mathrm{E(serial)}}$', fontsize=18)
    a.grid(True, which='major')
    a.legend(loc="lower right")
    f5.savefig("energySpeedup_{0}.pdf".format(exe[-3:]), format="pdf", bbox_inches="tight")


plt.show()

