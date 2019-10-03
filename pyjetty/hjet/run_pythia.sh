#!/bin/bash

for trig in "--py-minbias" "--t-min 6 --t-max 7 --py-pthatmin 6" "--t-min 12 --t-max 22 --py-pthatmin 12" "--t-min 20 --t-max 30 --py-pthatmin 20"
do
	for us in " " "--py-noMPI --py-noISR"
	do
		./hjet_simple.py ${trig} --charged --nev 10000 ${us} $@
	done
done
