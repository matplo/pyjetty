#!/bin/bash

files=$(find "$HOME/data/aleph/LEP1/1995" -name "*.aleph")
time parallel --bar ./invariant_mass.py -i {} -n -1 --oxo ::: ${files}
