#!/bin/bash -l
#SBATCH --nodes=32
#SBATCH --time=20:00

# Get execution directory
basedir=$(pwd)

# Get script directory
pushd $(dirname $0) > /dev/null
scriptdir=$(pwd)
popd > /dev/null

export OMP_NUM_THREADS=8

# Check directory of executables
if [ -z "$1" ]
then
    exedir=$basedir
else
    exedir=$(readlink -en $1)
fi


# Create result file
truncate -s 0 result.dat

echo " Executable       procs     kmax       Speedup         Max. speedup      Relative error     Serial runtime     Parallel runtime   Energy serial   Energey parallel" > result.dat
echo >> result.dat

for executable in pararealCPU pararealGPU
do
    if [ $executable == "pararealGPU" ]
    then
        async=--noasync
        export MPICH_RDMA_ENABLED_CUDA=1
    else
        async=
        export MPICH_RDMA_ENABLED_CUDA=0
        if [[ -z ${OMP_NUM_THREADS} ]]
        then
            export OMP_NUM_THREADS=8
        fi
    fi

    exe="${exedir}/${executable}"

    procs=4
    while [ $procs -le $SLURM_JOB_NUM_NODES ]
    do
        for k in 2 3 4
        do
            resdir="${basedir}/output_${executable}_${procs}_${k}"
            mkdir -p $resdir
            pushd $resdir > /dev/null

            echo "Running ${exe} with ${procs} nodes, kmax = ${k}"

            # Run application parallel
            aprun -N 1 -n $procs -d 8 "${exe}" --gridSize=32 --timeStepsFine=4096 --timeStepsCoarse=128 --nu=0.1 --endTime=0.10 --mode=parallel "--kmax=${k}" $async > outputParallel.log
            apid=$(sed -rn 's/^Application ([0-9]*).*/\1/p' < outputParallel.log)
            parallel=$(sed -rn 's/Parallel run time: ([0-9]*)/\1/p' < outputParallel.log)
            sleep 4
            energyParallel=$(sed -rn "s/.*apid: ${apid}.*\['energy_used', ([0-9]*)\]/\1/p" < "/scratch/daint/RUR/rur-$(date +%Y%m%d)")
            energyParallel=$(bc -l <<< "$energyParallel / .95 + 39.21*$procs*$parallel")

            # Run application serial
            aprun -N 1 -n 1 -d 8 "${exe}" --gridSize=32 --timeStepsFine=4096 --nu=0.1 --endTime=0.10 --mode=serial "--kmax=${k}" $async > outputSerial.log
            apid=$(sed -rn 's/^Application ([0-9]*).*/\1/p' < outputSerial.log)
            serial=$(sed -rn 's/Serial run time: ([0-9]*)/\1/p' < outputSerial.log)
            sleep 4
            energySerial=$(sed -rn "s/.*apid: ${apid}.*\['energy_used', ([0-9]*)\]/\1/p" < "/scratch/daint/RUR/rur-$(date +%Y%m%d)")
            energySerial=$(bc -l <<< "$energySerial / 0.95 + 39.21*$serial")

            # Run application compare
            aprun -N 1 -n $procs -d 8 "${exe}" --gridSize=32 --timeStepsFine=4096 --timeStepsCoarse=128 --nu=0.1 --endTime=0.10 --mode=compare "--kmax=${k}" $async > outputCompare.log
            speedup=$(sed -rn 's/Speedup: ([0-9]*)/\1/p' < outputCompare.log)
            maxspeedup=$(sed -rn 's/Maximal speedup: ([0-9]*)/\1/p' < outputCompare.log)
            error=$(sed -rn 's/Error at end: ([0-9]*)/\1/p' < outputCompare.log)

            # Print data
            printf "%12s      %4d      %2d      %.7e      %.7e      %.7e      %.7e      %.7e      %.7e      %.7e\n" $executable $procs $k $speedup $maxspeedup $error $serial $parallel $energySerial $energyParallel >> result.dat

            popd > /dev/null
        done
        procs=$(($procs*2))
    done
done

# Success.  Copy plotting script
plotscript="${scriptdir}/daintMultiplot.py"
cp $plotscript plot.py

