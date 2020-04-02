#!/bin/bash

# --py-eic-lowQ2 
# ./pythia_dc.py --py-eic --py-eic-cgamma --nev 1000 $@ --pythiaopts 421:mayDecay=false
# ./pythia_dc.py --py-eic --py-eic-cgamma --nev 1000 $@ 

./pythia_dc.py --py-eic --py-eic-cgamma --nev 1000 $@ --pythiaopts 421:mayDecay=false,Photon:Q2max=1.,Photon:Wmin=10.,PhaseSpace:pTHatMin=2.,Photon:ProcessType=0 --output dc_dnodecay.root
./pythia_dc.py --py-eic --py-eic-lowQ2  --nev 1000 $@ --output dc_all.root

