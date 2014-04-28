#!/bin/bash -l

usage()
{
    echo "Usage: 05-SpeedupCPU executable threads [kmax]"
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
jobscript="$scriptpath/05-SpeedupCPU.job"

echo "Using executable $exe"
echo "Running job script $jobscript"

# Start jobs
for nodes in 4 8 16 32 64 128; do
    sbatch --nodes $nodes --time 6:00:00 $jobscript $exe $2 $3
done

