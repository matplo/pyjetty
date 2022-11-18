OUTPUT_DIR='/Users/jamesmulligan/Analysis_Angularity/plot-angularity/hepdata'

DIRS="AngR02_ptbin1 AngR02_ptbin2 AngR02_ptbin3 AngR02_ptbin4 AngR04_ptbin1 AngR04_ptbin2 AngR04_ptbin3 AngR04_ptbin4"

for DIR in $DIRS ;
do
  NEW_DIR=${OUTPUT_DIR}/${DIR}
  mkdir -p $NEW_DIR
  mkdir -p $NEW_DIR/final_results
  mkdir -p $NEW_DIR/systematics
  scp -r james@hic.lbl.gov:/rstorage/alice/AnalysisResults/ang/$DIR/ang/final_results/*.root $NEW_DIR/final_results
  scp -r james@hic.lbl.gov:/rstorage/alice/AnalysisResults/ang/$DIR/ang/systematics/*.root $NEW_DIR/systematics
done

#scp -r james@hic.lbl.gov:/home/james/plot-angularity /Users/jamesmulligan/Analysis_Angularity/plot-angularity/hepdata
