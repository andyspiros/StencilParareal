#!/bin/bash -l

usage()
{
    echo "Usage: 01-Convergence executable target"
    echo "  target is either 'CPU' or 'GPU'"
}

# Parse parameters
if [ $# -ne 2 ]; then
    usage
    exit 1
fi

exe=$(readlink -e $1)
if [ ! -x $exe ]; then
    echo "Error: executable not found. Aborting"
    exit 2
fi

target="GPU"
if [ $2 == 'GPU' ]; then
    target="GPU"
elif [ $2 == 'CPU' ]; then
    target='CPU'
else
    echo "Target $2 not recognized"
    exit 3
fi

# Environment
pushd `dirname $0` > /dev/null
scriptpath=`pwd`
popd > /dev/null
jobscript="$scriptpath/01-Convergence.job"

echo "Using executable $exe with target $target"
echo "Running job script $jobscript"

# Check that directories do not exist
for nodes in 8 32 128; do
    rundir="run_${nodes}"
    if [ -d $rundir ]; then
        echo "Error: the directory ${rundir} exists alredy"
        echo "Removing"
        rm -rf $rundir
    fi
    mkdir $rundir
done

# Submit jobs
for nodes in 8 32 128; do
    rundir="run_${nodes}"
    pushd $rundir > /dev/null
    sbatch --time=2:30:00 --nodes=$nodes $jobscript $exe $target
    popd > /dev/null
done

