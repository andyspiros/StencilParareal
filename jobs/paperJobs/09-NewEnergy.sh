#!/bin/bash -l

usage()
{
    echo "Usage: 09-NewEnergy executableCPU executableGPU"
}

# Executables
exeCPU=$(readlink -e $1)
if [ ! -x $exeCPU ]; then
    echo "Error: executable not found. Aborting"
    exit 2
fi

exeGPU=$(readlink -e $2)
if [ ! -x $exeGPU ]; then
    echo "Error: executable not found. Aborting"
    exit 2
fi

# Environment
pushd `dirname $0` > /dev/null
scriptpath=`pwd`
popd > /dev/null
jobscriptCPU="$scriptpath/09-NewEnergyCPU.job"
jobscriptGPU="$scriptpath/09-NewEnergyGPU.job"

echo "Using executables:"
echo " - $exeCPU"
echo " - $exeGPU"
echo
echo "Running job scripts: $jobscriptCPU and $jobscriptGPU"
echo " - $jobscriptCPU"
echo " - $jobscriptGPU"
echo

# Start jobs CPU
mkdir -p CPU
pushd CPU
for nodes in 1 4 8 16 32 64 128; do
    sbatch --nodes $nodes -J "CPU-$nodes" --time 2:00:00 $jobscriptCPU $exeCPU 8 3
done
popd

# Start jobs GPU
mkdir -p GPU
pushd GPU
for nodes in 1 4 8 16 32 64 128; do
    sbatch --nodes $nodes -J "GPU-$nodes" --time 1:00:00 $jobscriptGPU $exeGPU 3
done
popd

