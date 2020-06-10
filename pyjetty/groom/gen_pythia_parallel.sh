#!/bin/bash

seeds="1111 2222 3333 4444 5555 6666 7777 8888 9999 121212"

pthatmin=50
parallel --bar --joblog pythia_write_groomed_parallel.log --keep-order --tag ./pythia_write_groomed.py --fj-R 0.4 --py-pthatmin ${pthatmin} --nev 5000 --py-seed {}  --output pythia_out_pthatmin${pthatmin}_{}.root ::: ${seeds}

#./pythia_write_groomed.py --py-pthatmin ${pthatmin} --nev 1000 --py-seed -1 --output pythia_out_pthatmin${pthatmin}_-1.root

parallel --bar --joblog pythia_write_groomed_parallel.log --keep-order --tag ./pythia_write_groomed.py --fj-R 0.4 --npart-min 300 --py-PbPb --nev 1000 --py-seed {} --output pythia_out_PbPb_{}.root ::: ${seeds}
