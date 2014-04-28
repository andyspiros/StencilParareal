#!/bin/bash -l

usage()
{
    echo "Usage: 04-speedupGPU executable target [kmax]"
    echo "  target is either CPU or GPU"
    echo "  kmax is by default 6"
}

# Executable
exe=$(readlink -e $1)
if [ ! -x $exe ]; then
    echo "Error: executable not found. Aborting"
    exit 2
fi

# Get target
if [ $2 == 'GPU' ]; then
    time="30:00"
elif [ $2 == 'CPU' ]; then
    time="6:00:00"
else
    echo "Target $2 not recognized"
    exit 3
fi

# Environment
pushd `dirname $0` > /dev/null
scriptpath=`pwd`
popd > /dev/null
jobscript="$scriptpath/04-Speedup.job"

echo "Using executable $exe"
echo "Running job script $jobscript"

# Start jobs
for nodes in 4 8 16 32 64 128; do
    sbatch --nodes $nodes --time $time $jobscript $exe $2 $3
done

