#!/bin/bash

for temp in {1..9}
do
   for fieldA in {0..2..1}
   do
      for fieldB in {0..9..1}
      do
         for n in {0..10..1}
         do
            echo T=0.${temp} h=0.${fieldA}${fieldB} n=${n}
            bin/xychain annealing --flagfile=./scripts/2dxy_flags.txt --initial_config=twisted --island_field_config=uniform --overrelax_freq=1 --annealing_steps=1 --mc_sweeps=2001 --novisuals --print_freq=2000 --temp_range=0.${temp},0.${temp} --island_local_field=0.${fieldA}${fieldB} > ${1}/h_0.${fieldA}${fieldB}_T_0.${temp}_${n}.txt
         done
      done
   done
done
