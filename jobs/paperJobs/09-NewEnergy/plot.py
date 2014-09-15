#!/usr/bin/python
import sys, os, re, glob
import numpy as np
import matplotlib as mpl

mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.size'] = 8.0

from matplotlib import pyplot as plt

keys = ['Node', 'Network', 'Blower', 'Device']
colors = ['r', 'b', 'm', 'g']

def getFiles(directory):
    if not directory.endswith('/'):
        directory += '/'

    files = [(1, directory + 'serial_0.log')]

    # Find parallel files
    pattern = re.compile('.*/parallel_([0-9]+)_0.log')
    for f in glob.glob(directory + 'parallel_*_0.log'):
        match = pattern.match(f)
        if match is not None:
            n = int(match.groups()[0])
            files.append((n, f))

    files.sort(key=lambda f: f[0])
    return files

def parseFileEnergy(fname):
    pattern = re.compile('^([a-zA-Z]+) energy *: ([0-9\.]+) J *\(([0-9e\+\.]+) W/node\)$')

    # Match lines
    result = {}
    with open(fname) as fs:
        for line in fs:
            match = pattern.match(line)
            if match is not None:
                [name, energy, power] = match.groups()
                result[name] = (float(energy)/1000., float(power))

    # Check that all data is parsed
    if not 'Node' in result \
    or not 'Device' in result \
    or not 'Network' in result \
    or not 'Blower' in result:
        print("Error in file", fname)
        sys.exit(1)

    return result

def plotEnergy(directory):
    files = getFiles(directory)
    nfiles = len(files)

    if 'CPU' in directory:
        CG = 'C'
    elif 'GPU' in directory:
        CG = 'G'
    else:
        print('Error: don\'t know whether directory is GPU or CPU')

    # Organize data
    data = {}
    if CG == 'C':
        mykeys = [k for k in keys if k != 'Device']
    else:
        mykeys = keys

    for k in mykeys:
        data[k] = []

    for f in files:
        d = parseFileEnergy(f[1])
        for k in data.keys():
            data[k].append((f[0], d[k]))

    # Create bars
    fig = plt.figure(figsize=(8./2.54, 8./2.54), dpi=300)
    ax = fig.add_axes((0.16, 0.11, 0.82, 0.82))

    width=0.7
    wh = width/2.
    left = np.r_[:nfiles]-wh
    bottom = np.zeros(nfiles)
    for k,c in zip(mykeys, colors):
        height = np.array([e[1][0] for e in data[k]])
        plt.bar(axes=ax, left=left, width=width, bottom=bottom, height=height, color=c, label=k)
        bottom = bottom + height

    # XTicks and grid
    labels = ['Serial'] + [str(f[0]) for f in files[1:]]
    ax.set_xticks(np.r_[:nfiles])
    ax.set_xticklabels(labels)

    for tic in ax.xaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False

    ax.yaxis.grid(True)

    # Axis labels and title
    ax.set_xlabel('Nodes')
    ax.set_ylabel('Energy (kJ)')

    # Legend
    ax.legend(loc='upper left')

    fig.savefig('Energy' + CG + 'PU.pdf', dpi=300)

if __name__ == '__main__':
    plotEnergy('CPU')
    plotEnergy('GPU')
