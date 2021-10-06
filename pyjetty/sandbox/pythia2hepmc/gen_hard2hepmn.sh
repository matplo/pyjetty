SPATH=$PYJETTY_DIR/pyjetty/sandbox/pythia2hepmc

$SPATH/pythia_gen_write_hepmc.py --py-cmnd ./common_pythia.cmnd ./pthat_bias.cmnd  ./hardQCD.cmnd --py-cmnd-out hard2hepmc.cmnd $@
