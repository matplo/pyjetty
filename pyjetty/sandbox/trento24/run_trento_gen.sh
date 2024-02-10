#!/bin/bash

nev=50000

./eec_lund.py --py-pthatmin 100 --jet-ptmin 100 --jet-ptmax 120 --nev ${nev} --output eec_lund_inc_100GeV_${nev}.parquet

./eec_lund.py --py-pthatmin 100 --jet-ptmin 100 --jet-ptmax 120 --nev ${nev} --py-hardQCDgluons --output eec_lund_gluons_100GeV_${nev}.parquet
./eec_lund.py --py-pthatmin 100 --jet-ptmin 100 --jet-ptmax 120 --nev ${nev} --py-hardQCDquarks --output eec_lund_quarks_100GeV_${nev}.parquet
./eec_lund.py --py-pthatmin 100 --jet-ptmin 100 --jet-ptmax 120 --nev ${nev} --py-hardQCDuds    --output eec_lund_uds_100GeV_${nev}.parquet

./eec_lund.py --py-pthatmin 100 --jet-ptmin 100 --jet-ptmax 120 --nev ${nev} --py-hardQCDcharm  --output eec_lund_charm_100GeV_${nev}.parquet
./eec_lund.py --py-pthatmin 100 --jet-ptmin 100 --jet-ptmax 120 --nev ${nev} --py-hardQCDbeauty --output eec_lund_beauty_100GeV_${nev}.parquet
# ./eec_lund.py --py-pthatmin 60 --jet-ptmin 60 --jet-ptmax 80 --nev 1000


./eec_lund.py --py-pthatmin 100 --jet-ptmin 100 --jet-ptmax 120 --nev ${nev} --py-hardQCDcharm  --stable-charm --output eec_lund_stable_charm_100GeV_${nev}.parquet
./eec_lund.py --py-pthatmin 100 --jet-ptmin 100 --jet-ptmax 120 --nev ${nev} --py-hardQCDbeauty --stable-beauty --output eec_lund_stable_beauty_100GeV_${nev}.parquet


