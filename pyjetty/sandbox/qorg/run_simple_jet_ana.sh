#!/bin/bash

./simple_jet_ana.py -i ./z_jet_files.txt
./simple_jet_ana.py -i ./quark_jet_files.txt
./simple_jet_ana.py -i ./glue_jet_files.txt
./simple_jet_ana.py -i ./z_jet_files.txt --ptree tree_Particle_gen_mDT0.04 --output z_jet_files_jets_mDT0.04.root