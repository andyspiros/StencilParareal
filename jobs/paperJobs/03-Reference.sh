#!/bin/bash -l

usage()
{
    echo "Usage: 03-Reference executable target"
    echo "  target is either 'CPU' or 'GPU'"
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
jobscript="$scriptpath/03-Reference.job"

echo "Using executable $exe"
echo "Running job script $jobscript"

# Start jobs
for nodes in 4 8 16 32 64 128; do
    sbatch --nodes $nodes --time 4:00:00 $jobscript $exe $2
done

