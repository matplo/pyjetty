#!/bin/bash

parallel --keep-order --jobs 8 --tag ./cstoy.py --nev 10000 --noue --pthatmin 100 --alpha {1} --dRmax {2} $@ ::: 0 0.5 1 2 4 ::: 0.2 0.4 0.8 2.0
