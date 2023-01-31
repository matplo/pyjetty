#!/bin/bash

R0=0.2
nev=10000
./test_jaxis_ue.py /data/rstorage/alice_data_qr.txt --nev ${nev} \
	--jet-pt-min 20 --jet-pt-max 40 \
	--jet-R0 ${R0} \
	--output jaxis_test_20_40_${R0}.root &
./test_jaxis_ue.py /data/rstorage/alice_data_qr.txt --nev ${nev} \
	--jet-pt-min 40 --jet-pt-max 60 \
	--jet-R0 ${R0} \
	--output jaxis_test_40_60_${R0}.root &
./test_jaxis_ue.py /data/rstorage/alice_data_qr.txt --nev ${nev} \
	--jet-pt-min 60 --jet-pt-max 80 \
	--jet-R0 ${R0} \
	--output jaxis_test_60_80_${R0}.root &

R0=0.4
nev=10000
./test_jaxis_ue.py /data/rstorage/alice_data_qr.txt --nev ${nev} \
	--jet-pt-min 20 --jet-pt-max 40 \
	--jet-R0 ${R0} \
	--output jaxis_test_20_40_${R0}.root &
./test_jaxis_ue.py /data/rstorage/alice_data_qr.txt --nev ${nev} \
	--jet-pt-min 40 --jet-pt-max 60 \
	--jet-R0 ${R0} \
	--output jaxis_test_40_60_${R0}.root &
./test_jaxis_ue.py /data/rstorage/alice_data_qr.txt --nev ${nev} \
	--jet-pt-min 60 --jet-pt-max 80 \
	--jet-R0 ${R0} \
	--output jaxis_test_60_80_${R0}.root &
