#!/bin/bash -l
exe=$1

nu0=0.1
gridSize=128
endTime=0.1

runidx=0

echo "Order freq tSteps    AbsErr          RelErr"

#for order in 1 4
for order in 4
do
    for nufreq in 100 10 1
    do
        timeSteps=1024

        while [ "${timeSteps}" -le 8192 ]
        do
            
            runidx=$(($runidx+1))
            outfile="run_${runidx}.log"
            
            # Run
            aprun $exe --nu0 $nu0 --nufreq $nufreq --gridSize $gridSize --endTime $endTime --timeSteps $timeSteps --order $order > $outfile

            # Parse
            errorabs=$(sed -rn 's/^ - Process .* error=(.*)$/\1/p' < $outfile)
            errorrel=$(sed -rn 's/^ Total relative error: (.*)$/\1/p' < $outfile)

            # Output
            printf "%5d %4d %6d %.7e %.7e\n" $order $nufreq $timeSteps $errorabs $errorrel

            # Next timeSteps
            timeSteps=$((timeSteps*2))
        done
    done
done

