#!/bin/bash

if [ "$1" != "" ]; then
  TASK_ID=$1
  echo "Task ID: $TASK_ID"
else
  echo "Wrong command line arguments"
fi

cd
echo $PWD
source .bashrc

# Enter aliroot environment (assumes I already have a valid token) and call download script
alidock exec alienv setenv AliRoot/latest -c /mnt/pyjetty/pyjetty/alice_analysis/slurm/utils/james/downloadData.sh $TASK_ID
