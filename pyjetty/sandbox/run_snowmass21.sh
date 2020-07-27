#!/bin/bash

cmnd="$PYJETTY_DIR/pyjetty/sandbox/snomass21_sim.py --nev 30000 --py-biasref 10 --py-biaspow 4 "
$cmnd --output sm21_inc.root

$cmnd --output sm21_glue.root --py-hardQCDgluons

$cmnd --output sm21_quark.root --py-hardQCDquarks

$cmnd --output sm21_charm.root --py-hardQCDcharm

$cmnd --output sm21_beauty.root --py-hardQCDbeauty
