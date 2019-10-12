#!/bin/bash

# parallel --keep-order --jobs 8 --tag ./cstoy.py --nev 10000 --noue --pthatmin 100 --alpha {1} --dRmax {2} $@ ::: 0 0.5 1 2 4 ::: 0.2 0.4 0.8 2.0

# parallel --keep-order --jobs 8 --tag ./cstoy_tw.py --embed /Users/ploskon/data/alice/LHC18qr/local_PbPb_file_list.txt --overwrite --nev 50000 --dRmax {1} --zcut {2} :: 0.1 0.2 0.3 0.4 0.8 ::: 0.1 0.15 0.2 0.25 0.35

# parallel --keep-order --jobs 8 --tag ./cstoy_tw.py --embed /Users/ploskon/data/alice/LHC18qr/local_PbPb_file_list.txt --nev 50000 --dRmax {} --zcut 0.2 ::: 0.1 0.2 0.25 0.3 0.4 0.8
parallel --keep-order --jobs 8 --tag ./cstoy_tw.py --embed /Users/ploskon/data/alice/LHC18qr/local_PbPb_file_list.txt --nev 50000 --dRmax 0.25 --zcut {} $@ ::: 0.1 0.15 0.2 0.25 0.35
