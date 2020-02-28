#! /bin/bash
#
# Script to merge output ROOT files from all pt-hat bins together, in stages
JOB_ID=35198
OUTPUT_DIR="/rstorage/u/alice/AnalysisResults/ang/$JOB_ID"

# Merge all output files from each pt-hat bin
# Stage 2
N=10
mkdir -p $OUTPUT_DIR/Stage2
for K in $(seq 1 $N);
do
  INDEX1=$((2*K-1))
  INDEX2=$((2*K))
  hadd -f $OUTPUT_DIR/Stage2/AnalysisResults$K.root $OUTPUT_DIR/Stage1/$INDEX1/*.root $OUTPUT_DIR/Stage1/$INDEX2/*.root
done

# Stage 3
N=5
mkdir -p $OUTPUT_DIR/Stage3
for K in $(seq 1 $N);
do
  INDEX1=$((2*K-1))
  INDEX2=$((2*K))
  hadd -f $OUTPUT_DIR/Stage3/AnalysisResults$K.root $OUTPUT_DIR/Stage2/AnalysisResults$INDEX1.root $OUTPUT_DIR/Stage2/AnalysisResults$INDEX2.root
done

# Stage 4
mkdir -p $OUTPUT_DIR/Stage4
hadd -f $OUTPUT_DIR/Stage4/AnalysisResults1.root $OUTPUT_DIR/Stage3/AnalysisResults1.root $OUTPUT_DIR/Stage3/AnalysisResults2.root
hadd -f $OUTPUT_DIR/Stage4/AnalysisResults2.root $OUTPUT_DIR/Stage3/AnalysisResults3.root $OUTPUT_DIR/Stage3/AnalysisResults4.root $OUTPUT_DIR/Stage3/AnalysisResults5.root

# Final merge
hadd -f $OUTPUT_DIR/AnalysisResultsFinal.root $OUTPUT_DIR/Stage4/AnalysisResults1.root $OUTPUT_DIR/Stage4/AnalysisResults2.root
