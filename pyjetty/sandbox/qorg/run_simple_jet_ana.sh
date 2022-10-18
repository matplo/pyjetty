#!/bin/bash

find ${PWD}/v2_z_jets -name "*.root" | tee z_jet_files.txt
./simple_jet_ana.py -i ./z_jet_files.txt &
./simple_jet_ana.py -i ./z_jet_files.txt --ptree tree_Particle_gen_mDT0.04 --output z_jet_files_jets_mDT0.04.root &

find ${PWD}/v2_q_jets -name "*.root" | tee quark_jet_files.txt
./simple_jet_ana.py -i ./quark_jet_files.txt &

find ${PWD}/v2_g_jets -name "*.root" | tee glue_jet_files.txt
./simple_jet_ana.py -i ./glue_jet_files.txt &


