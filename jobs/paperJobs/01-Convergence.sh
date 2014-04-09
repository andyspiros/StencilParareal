#! /bin/bash -l
#SBATCH --account=h03
#SBATCH --ntasks-per-node=1

# Compare convergence of Parareal against serial fine solution for
# omega = 0 and omega = 100 to verify minor impact of time-dependent
# coefficient on Parareal convergence.

# x: number of iterations k = 1, ..., 8
# y: relative difference Parareal to fine serial solution in iteration k
# lines  for nodes = 8, 32, 128

usage()
{
    echo "Usage: 01-Convergence executable [target]"
    echo "  target is either 'CPU' or 'GPU'"
}

# Retrieve executable
if [ $# -lt 1 ] || [ $# -gt 2 ]; then
    usage
    exit 1
fi

if [ -x $1 ]; then
    exe=$1
else
    echo "Error: executable not found. Aborting"
    exit 2
fi

target="GPU"
if [ $# -gt 1 ]; then
    if [ $2 == 'GPU' ]; then
        target="GPU"
    elif [ $2 == 'CPU' ]; then
        target='CPU'
    else
        echo "Target $2 not recognized"
        exit 3
    fi
fi

if [ $target == 'GPU' ]; then
    export MPICH_RDMA_ENABLED_CUDA=1
    async=--noasync
else
    async=
fi

nodes="${SLURM_JOB_NUM_NODES}"

gridSize=128
endTime=0.1
timeStepsCoarse=2048
timeStepsFine=32768
nu0=0.1
nufreq=100

args="--gridSize $gridSize --endTime $endTime --timeStepsCoarse $timeStepsCoarse --timeStepsFine $timeStepsFine --nu0 $nu0 --nufreq $nufreq $async"

resultfile="result_$nodes.dat"
truncate -s 0 $resultfile

# Perform execution
for k in $(seq 1 8); do
    outfile="run_$nodes_$k.log"
    aprun -n $nodes -N 1 $exe $args --kmax $k > $outfile

    error=$(sed -rn 's/^Error at end ([0-9.]+)$/\1/p' < $outfile)
    printf '%.7e  ' $error >> $resultfile
done

