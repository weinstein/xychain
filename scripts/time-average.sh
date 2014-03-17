#!/bin/bash


E_SUM=0
N_MEAS=0
curtemp=$(head -n 1 ${1} | cut -d" " -f 3)
prevtemp=$curtemp
while read -r line
do
	curtemp=$(echo $line | cut -d" " -f 3)
	if [ $curtemp != $prevtemp ]
	then
		energy_avg=$(echo "$E_SUM/$N_MEAS" | bc -l)
		echo "$prevtemp $energy_avg"
		E_SUM=0
		N_MEAS=0
		prevtemp=$curtemp
   	fi
	energy=$(echo $line | cut -d" " -f 4)
	E_SUM=$(echo "$E_SUM + $energy" | bc -l)
	N_MEAS=$(($N_MEAS + 1))
done < "${1}"
