#!/bin/bash

$HEPPY_DIR/heppy/examples/pythia_gen_write_hepmc.py --nev 1000 
$PYJETTY_DIR/pyjetty/sandbox/hepmc2antuple.py -i pythia_gen_test_hepmc2.dat  -o test_hepmc_convert.root --nev 1000 --as-data
echo "./test_hepmc_convert.root" > test_hepmc_convert_list.txt
$PYJETTY_DIR/pyjetty/cstoy/csdata.py test_hepmc_convert_list.txt --overwrite --output test_analyze_hepmc2root.root

$PYJETTY_DIR/pyjetty/sandbox/hepmc2antuple_tn.py -i pythia_gen_test_hepmc2.dat  -o test_hepmc_convert_tn.root --nev 1000 --as-data
echo "./test_hepmc_convert_tn.root" > test_hepmc_convert_list.txt
$PYJETTY_DIR/pyjetty/cstoy/csdata.py test_hepmc_convert_list.txt --overwrite --output test_analyze_hepmc2root_tn.root
