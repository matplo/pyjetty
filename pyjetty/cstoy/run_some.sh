#!/bin/bash

./cstoy.py --nev 10000 --noue --pthatmin 100 --alpha 0.5 --dRmax 0.25
./cstoy.py --nev 10000 --noue --pthatmin 100 --alpha 0.0 --dRmax 0.25
./cstoy.py --nev 10000 --noue --pthatmin 100 --alpha 1.0 --dRmax 0.25

./cstoy.py --nev 10000 --noue --pthatmin 100 --alpha 0.5 --dRmax 0.4
./cstoy.py --nev 10000 --noue --pthatmin 100 --alpha 0.0 --dRmax 0.4
./cstoy.py --nev 10000 --noue --pthatmin 100 --alpha 1.0 --dRmax 0.4

./cstoy.py --nev 10000 --noue --pthatmin 100 --alpha 0.5 --dRmax 1.0
./cstoy.py --nev 10000 --noue --pthatmin 100 --alpha 0.0 --dRmax 1.0
./cstoy.py --nev 10000 --noue --pthatmin 100 --alpha 1.0 --dRmax 1.0

