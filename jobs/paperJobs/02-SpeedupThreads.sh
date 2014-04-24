#!/bin/bash -l

usage()
{
    echo "Usage: 02-SpeedupThreads executable [k]"
    echo "  k defaults to 6"
}

# Retrieve executable
if [ $# -lt 1 ] || [ $# -gt 2 ]; then
    usage
    exit 1
fi

exe=$(readlink -e $1)
if [ ! -x $exe ]; then
    echo "Error: executable not found. Aborting"
    exit 2
fi

# Get kmax
if [ $# -lt 2 ]; then
    kmax=6
else
    kmax=$2
fi

# Environment
pushd `dirname $0` > /dev/null
scriptpath=`pwd`
popd > /dev/null
jobscript="$scriptpath/02-SpeedupThreads.job"

echo "Using executable $exe"
echo "Running job script $jobscript"

# Check that directories do not exist
for nodes in 8 32 128; do
    rundir="run_${nodes}"
    if [ -d $rundir ]; then
        echo "Warning: the directory ${rundir} exists alredy"
        echo "Removing"
        rm -rf $rundir
    fi
    mkdir $rundir
done

# Submit jobs
for nodes in 8 32 128; do
    rundir="run_${nodes}"
    pushd $rundir > /dev/null
    sbatch --time=30:00 --nodes=$nodes $jobscript $exe $kmax
    popd > /dev/null
done

