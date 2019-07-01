#!/bin/bash

sfiles=$(find ./epps16-pp-5000-bias -name "*.dat")
for fn in ${sfiles}
do
	echo "[i] processing ${fn}"
	./reco_jets_epps16.py --read ${fn}
done
