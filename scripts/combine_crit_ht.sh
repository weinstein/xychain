#!/bin/bash

for temp in {1..9}
do
   outfile=${2}/T_0.${temp}.txt
   echo > $outfile
   for fieldA in {0..2..1}
   do
      for fieldB in {0..9..1}
      do
         for n in {0..10..1}
         do
            infile=${1}/h_0.${fieldA}${fieldB}_T_0.${temp}_${n}.txt
            echo T=0.${temp} h=0.${fieldA}${fieldB} n=${n}
            echo $infile into $outfile
            echo 0.${fieldA}${fieldB} $(tail -n 1 ${infile}) >> $outfile
         done
      done
   done
done
