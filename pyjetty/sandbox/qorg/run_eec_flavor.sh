#!/bin/bash

# for flavor in gluons uds charm beauty
# do
# 	./eec_flavor_dep_example.py --nev 1000 --py-pthatmin 20 --jet-ptmin 20 --jet-ptmax 40 --py-hardQCD${flavor} --output eec_${flavor}.root &
# done
# 
# parallel ./eec_flavor_dep_example.py --nev 10000 --py-pthatmin 20 --jet-ptmin 20 --jet-ptmax 40 --py-hardQCD{1} --output eec_{1}.root ::: gluons uds charm beauty

parallel --progress --jobs 8 ./eec_flavor_dep_example.py --nev 100 --py-pthatmin {2} --jet-ptmin {2} --jet-ptmax '$(expr {2} + 50)' --py-hardQCD{1} --output eec_{1}_{2}.root "2>&1 >" eec_{1}_{2}.log ::: gluons uds charm beauty ::: 20 100 500

