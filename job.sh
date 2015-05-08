#!/bin/bash

begin=-1
end=1
step=1120
P=7
for I in {559..560}; do
    #dc=$(echo "scale=4;$begin*(-1)+($end-$begin*(-1))/$step*$I" | bc | awk '{printf "%f", $0}')
    #dc=$(echo "scale=10;$begin+($end- $begin)/$step*$I" | bc | awk '{printf "%f", $0}')
    #echo $dc
    #sed '7s/=0.0/=1' parameter0.txt > parameter.txt
    #awk '{if (NR==7) print "dc='$dc'"; else print $0}' parameter0.txt > parameter.txt
echo "#PBS -l nodes=1:ppn=1 
#PBS -j oe
#PBS -o log
#PBS -l walltime=240:00:00
cd /home/glwang/STT/

./a.out $I $P $begin $end $step
txt2mat nx$I.txt nx$I nx$I.mat
txt2mat ny$I.txt ny$I ny$I.mat
txt2mat nz$I.txt nz$I nz$I.mat
rm nx$I.txt ny$I.txt nz$I.txt
" > single.sh
done
