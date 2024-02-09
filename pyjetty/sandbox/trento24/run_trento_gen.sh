#!/bin/bash

nev=20000

./eec_lund.py --py-pthatmin 100 --jet-ptmin 100 --jet-ptmax 120 --nev 10000 --output eec_lund_inc_100GeV_${nev}.parquet

./eec_lund.py --py-pthatmin 100 --jet-ptmin 100 --jet-ptmax 120 --nev 10000 --py-hardQCDgluons --output eec_lund_gluons_100GeV_${nev}.parquet
./eec_lund.py --py-pthatmin 100 --jet-ptmin 100 --jet-ptmax 120 --nev 10000 --py-hardQCDquarks --output eec_lund_quarks_100GeV_${nev}.parquet
./eec_lund.py --py-pthatmin 100 --jet-ptmin 100 --jet-ptmax 120 --nev 10000 --py-hardQCDusd    --output eec_lund_uds_100GeV_${nev}.parquet

./eec_lund.py --py-pthatmin 100 --jet-ptmin 100 --jet-ptmax 120 --nev 10000 --py-hardQCDcharm  --output eec_lund_charm_100GeV_${nev}.parquet
./eec_lund.py --py-pthatmin 100 --jet-ptmin 100 --jet-ptmax 120 --nev 10000 --py-hardQCDbeauty --output eec_lund_beauty_100GeV_${nev}.parquet
# ./eec_lund.py --py-pthatmin 60 --jet-ptmin 60 --jet-ptmax 80 --nev 1000
