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

truncate result.dat

procs=4
while [ $procs -le $SLURM_JOB_NUM_NODES ]
do
    coarseTS=$((gs*4))
    coarseTS=$(($coarseTS < $procs ? $procs : $coarseTS))
    fineTS=$(($gs*128))
    printf "%4d   " $procs

    for threadnum in 1 2 4 8
    do
        for k in 2 4
        do

            export OMP_NUM_THREADS=$threadnum
            aprun -N 1 -n $procs -d 8 "${exe}" --gridSize=$gs --timeStepsFine=$fineTS --timeStepsCoarse=$coarseTS --nu=0.1 --endTime=0.10 --mode=parallel "--kmax=${k}" > outputParallel.log
            parallel=$(sed -rn 's/Parallel run time: ([0-9]*)/\1/p' < outputParallel.log)

            echo "Run aprun -N 1 -n ${procs} -d 8 ${exe} --gridSize=${gs} --timeStepsFine=${fineTS} --timeStepsCoarse=${coarseTS} --nu=0.1 --endTime=0.10 --mode=parallel --kmax=${k}"

            printf "%.7e " $parallel >> result.dat
        done
        printf "   "  >> result.dat
    done
    printf "\n" >> result.dat

    procs=$(($procs*2))
done
