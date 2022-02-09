#!/bin/bash

if [ -z "${PYJETTY_DIR}" ]; then
    echo "PYJETTY_DIR is not set - stop here."
else
    $PYJETTY_DIR/pyjetty/alephana/lep_gen_alice_data_struct.py --nev 100000 --py-pthatmin 5 --output lep_pthatmin5_nev100k.root
fi

