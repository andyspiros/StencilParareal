#!/bin/bash

for i in 4 8 16 32 64 128; do
    s=$(sed -rn 's/^Speedup: (.*)$/\1/p' < "run_${i}_0.log");
    m=$(sed -rn 's/^Maximal speedup: (.*)$/\1/p' < "run_${i}_0.log");
    printf "%3d  %e  %e\n" $i $s $m;
done > speedup.dat

for i in 4 8 16 32 64 128; do
    s=$(sed -rn 's/^Serial run time: (.*)$/\1/p' < "run_${i}_0.log");
    p=$(sed -rn 's/^Parallel run time: (.*)$/\1/p' < "run_${i}_0.log");
    printf "%3d  %e  %e\n" $i $s $p;
done > runtime.dat
