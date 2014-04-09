#! /bin/bash -l

basedir=$(pwd)

# Check directory of executables
if [ -z "$1" ]
then
    exedir=$basedir
else
    exedir=$(readlink -en $1)
fi
exe=$exedir/pararealCPU

if [ -z "$2" ]
then
    gs=32
else
    gs=$2
fi

# Just for debugging
if [ -z $SLURM_JOB_NUM_NODES ]
then
    SLURM_JOB_NUM_NODES=32
fi

procs=4
while [ $procs -le $SLURM_JOB_NUM_NODES ]
do
    coarseTS=$((gs*4))
    coarseTS=$(($coarseTS < $procs ? $procs : $coarseTS))
    fineTS=$(($gs*128))
    printf "%4d   " $procs

    for async in "" "--noasync"
    do
        for k in 2 3 4
        do
            aprun -N 1 -n $procs "${exe}" --gridSize=$gs --timeStepsFine=$fineTS --timeStepsCoarse=$coarseTS --nu=0.1 --endTime=0.10 --mode=parallel "--kmax=${k}" $async > outputParallel.log
            parallel=$(sed -rn 's/Parallel run time: ([0-9]*)/\1/p' < outputParallel.log)

            printf "%.7e " $parallel
        done
        printf "   "
    done
    printf "\n"

    procs=$(($procs*2))
done
