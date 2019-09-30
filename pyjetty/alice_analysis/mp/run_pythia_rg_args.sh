#!/bin/bash

for jetR in 0.4 0.6
do
	for pthatmin in 20 40 60 80
	do
		for mpilevel in "       " "--noMPI"
			do
				for isrlevel in "       " "--noISR"
				do
					# echo "./pythia_rg.py --nev 10000 --jetR ${jetR} --pthatmin ${pthatmin} --ecm 5020 ${mpilevel} ${isrlevel}"
					./pythia_rg.py --nev 10000 --jetR ${jetR} --pthatmin ${pthatmin} --ecm 5020 ${mpilevel} ${isrlevel}
				done
			done
	done
done