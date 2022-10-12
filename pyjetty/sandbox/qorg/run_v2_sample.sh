#!/bin/bash

nev=10000

this_dir=${PWD}

mkdir v2_z_jets
cd v2_z_jets
${this_dir}/pythia_gen_jets_ml_v2.py --py-noue --py-seed 123456 --Zjet --ZjetR 0.8 --mDT --mDTzcut 0.04 --py-cmnd ${this_dir}/pythia_gen_qorg_Zjet_master.cmnd --jet-ptmin 500 --jet-ptmax 550 --nev ${nev}
cd ${this_dir}

mkdir v2_g_jets
cd v2_g_jets
${this_dir}/pythia_gen_jets_ml_v2.py --py-noue --py-seed 123456 --py-hardQCDgluons --py-pthatmin 500 --jet-ptmin 500 --jet-ptmax 550 --nev ${nev}
cd ${this_dir}

mkdir v2_q_jets
cd v2_q_jets
${this_dir}/pythia_gen_jets_ml_v2.py --py-noue --py-seed 123456 --py-hardQCDquarks --py-pthatmin 500 --jet-ptmin 500 --jet-ptmax 550 --nev ${nev}
cd ${this_dir}
