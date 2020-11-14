#!/bin/bash

BASE_DIR=/home/james/pyjetty/pyjetty/alice_analysis/generation/herwig

for BIN in $(seq 1 20);
do
    echo "Generating bin: $BIN"
    cd $BASE_DIR/run/$BIN
    Herwig read $BASE_DIR/config/$BIN/LHC_5020_MPI.in
done
