#!/bin/bash -l

usage()
{
    echo "Usage: 04-SpeedupGPU executable [kmax]"
    echo "  kmax is by default 6"
}

# Executable
exe=$(readlink -e $1)
if [ ! -x $exe ]; then
    echo "Error: executable not found. Aborting"
    exit 2
fi

# Environment
pushd `dirname $0` > /dev/null
scriptpath=`pwd`
popd > /dev/null
jobscript="$scriptpath/04-SpeedupGPU.job"

echo "Using executable $exe"
echo "Running job script $jobscript"

# Start jobs
for nodes in 4 8 16 32 64 128; do
    echo sbatch --nodes $nodes --time 30:00 $jobscript $exe $2
    sbatch --nodes $nodes --time 30:00 $jobscript $exe $2
done

