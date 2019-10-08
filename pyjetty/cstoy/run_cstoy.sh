#!/bin/bash

parallel --keep-order --tag ./cstoy_tw.py --overwrite --py-pthatmin 80 --nev 5000 --dRmax 0.25 --py-seed {} ::: 22 33 44 55 66 77
hadd -f output_alpha_0_dRmax_0.25_SDzcut_0.1_seed_merged.root output_alpha_0_dRmax_0.25_SDzcut_0.1_seed_??.root
