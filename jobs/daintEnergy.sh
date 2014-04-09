#!/bin/bash -l

basedir=$(pwd)

# Check directory of executables
if [ -z "$1" ]
then
    exedir=$basedir
else
    exedir=$(readlink -en $1)
fi

# Open file with data
truncate -s 0 energy.dat

echo " Executable       procs     kmax     Parallel comp     Parallel comm     Parallel blower    Serial comp       Serial comm       Serial blower"
echo

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

            echo "Execute $executable on $procs nodes with kmax=$k"

            # Run application parallel
            aprun -N 1 -n $procs "${exe}" --gridSize=32 --timeStepsFine=4096 --timeStepsCoarse=128 --nu=0.1 --endTime=0.10 --mode=parallel "--kmax=${k}" $async > outputParallel.log
            apid=$(sed -rn 's/^Application ([0-9]*).*/\1/p' < outputParallel.log)
            parallelRT=$(sed -rn 's/Parallel run time: ([0-9]*)/\1/p' < outputParallel.log)
            sleep 4

            echo " - apid of parallel: $apid"
            echo " - runtime of parallel: $parallelRT"

            energyParallel=$(sed -rn "s/.*apid: ${apid}.*\['energy_used', ([0-9]*)\]/\1/p" < "/scratch/daint/RUR/rur-$(date +%Y%m%d)")
            eParComp=$(bc -l <<< "$energyParallel / 0.95")
            eParComm=$(bc -l <<< "$procs * 25 * $parallelRT")
            eParBlow=$(bc -l <<< "$procs * 14.21 * $parallelRT")

            echo " - Energy of parallel: $energyParallel"

            # Run application serial
            aprun -N 1 -n 1 "${exe}" --gridSize=32 --timeStepsFine=4096 --nu=0.1 --endTime=0.10 --mode=serial > outputSerial.log
            apid=$(sed -rn 's/^Application ([0-9]*).*/\1/p' < outputSerial.log)
            serialRT=$(sed -rn 's/Serial run time: ([0-9]*)/\1/p' < outputSerial.log)
            sleep 4

            echo " - apid of serial: $apid"
            echo " - runtime of serial: $serialRT"

            energySerial=$(sed -rn "s/.*apid: ${apid}.*\['energy_used', ([0-9]*)\]/\1/p" < "/scratch/daint/RUR/rur-$(date +%Y%m%d)")
            eSerComp=$(bc -l <<< "$energySerial / 0.95")
            eSerComm=$(bc -l <<< "25 * $serialRT")
            eSerBlow=$(bc -l <<< "14.21 * $serialRT")

            echo " - Energy of serial: $energySerial"

            popd > /dev/null

            # Print data
            printf "%12s      %4d      %2d      %.7e      %.7e      %.7e      %.7e      %.7e      %.7e\n" $executable $procs $k $eParComp $eParComm $eParBlow $eSerComp $eSerComm $eSerBlow >> energy.dat

        done
        procs=$(($procs*2))
    done
done


