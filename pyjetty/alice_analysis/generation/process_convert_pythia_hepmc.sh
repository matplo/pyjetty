#! /bin/bash

if [ "$1" != "" ]; then
  FILEPATH=$1
  #echo "Input file: $INFILE_NAME"
else
  echo "Wrong command line arguments"
fi

# Set output location based on input filename
INFILE_NAME="${FILEPATH##*/}"
echo "Input file: $INFILE_NAME"

OUTFILE_NAME="${INFILE_NAME%%.*}.root"
OUTPUT_DIR=/rstorage/james/groom/$OUTFILE_NAME
echo "Output file: $OUTPUT_DIR"

python hepmc2antuple_tn.py -i $FILEPATH -o $OUTPUT_DIR -g pythia --no-progress-bar
