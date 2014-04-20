#! /bin/bash -l
#SBATCH --account=c01
#SBATCH --ntasks-per-node=1

# Compare convergence of Parareal against serial fine solution for
# omega = 0 and omega = 100 to verify minor impact of time-dependent
# coefficient on Parareal convergence.

# x: number of iterations k = 1, ..., 8
# y: relative difference Parareal to fine serial solution in iteration k
# lines  for nodes = 8, 32, 128

usage()
{
    echo "Usage: 01-Convergence executable target"
    echo "  target is either 'CPU' or 'GPU'"
}

# Retrieve executable
if [ $# -ne 2 ]; then
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
if [ $2 == 'GPU' ]; then
    target="GPU"
elif [ $2 == 'CPU' ]; then
    target='CPU'
else
    echo "Target $2 not recognized"
    exit 3
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

args="--gridSize $gridSize --endTime $endTime --timeStepsCoarse $timeStepsCoarse --timeStepsFine $timeStepsFine --nu0 $nu0 $async"

resultfile_0="result_${nodes}_0.dat"
resultfile_100="result_${nodes}_100.dat"
truncate -s 0 $resultfile_100 $resultfile_0

# Perform execution
for k in $(seq 1 8); do
    outfile0="run_${nodes}_${k}_0.log"
    aprun -n $nodes -N 1 $exe $args --nufreq 0 --kmax $k > $outfile0

    outfile100="run_${nodes}_${k}_100.log"
    aprun -n $nodes -N 1 $exe $args --nufreq 100 --kmax $k > $outfile100

    error=$(sed -rn 's/^Error at end: ([0-9\.]+)$/\1/p' < $outfile0)
    printf '%.7e  ' $error >> $resultfile0

    error=$(sed -rn 's/^Error at end: ([0-9\.]+)$/\1/p' < $outfile100)
    printf '%.7e  ' $error >> $resultfile100
done

