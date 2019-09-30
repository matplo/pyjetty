#!/bin/bash

cdir=${PWD}
for dd in pthatmin20 pthatmin40 pthatmin60 pthatmin80
do
	cd ${dd}
	tdraw_cfg.py --clean ../tdraw_thetag.cfg 
	cd ${cdir}
done
cd ${cdir}
