#!/bin/bash

# for flavor in gluons uds charm beauty
# do
# 	./eec_flavor_dep_example.py --nev 1000 --py-pthatmin 20 --jet-ptmin 20 --jet-ptmax 40 --py-hardQCD${flavor} --output eec_${flavor}.root &
# done
# 
parallel ./eec_flavor_dep_example.py --nev 10000 --py-pthatmin 20 --jet-ptmin 20 --jet-ptmax 40 --py-hardQCD{1} --output eec_{1}.root ::: gluons uds charm beauty


