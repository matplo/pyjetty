### Analysis framework

This directory contains a simple analysis framework to analyze ROOT
trees using the python package heppy.

The framework has two stages:
1. Processing: Loop over the data, performing all desired event-by-event
      calculations, and save output histograms.
2. Post-processing: Perform analysis on the ensemble distributions,
      and plot final histograms.

The processing machinery is contained in the `process`:
  - `process/base`: Contains several shared common classes
      - the root base class: `base.py`
      - the processing base class: `process_base.py`
      - an IO class: `process_io.py`
      - a process utility class: `process_utils.py`
  - `process/user`: Contains user processing classes, which should derive
    from `process_base.py`. Each user should implement
    their own analysis here, leveraging the functionalities
    in `process/base`.

The directory `slurm` contains scripts to run the processing stage
with the batch system (`slurm/sbatch`), and some utility scripts
that are useful for merging/scaling data after processing (`slurm/utils`).

The directory `config` contains yaml configs that are used for both the
processing and post-processing codes.

The directory `analysis` contains the post-processing machinery:
  - `analysis/base`: Contains several shared common classes
  - the root base class: `base.py`
  - the analysis base class: `analysis_base.py`
  - an analysis utility class: `analysis_utils.py`
  - `analysis/user`: Contains user analysis classes, which should derive
    from `analysis/base.py`. Each user should implement
    their own analysis here, leveraging the functionalities in `analysis/base`.

### Example: 

To process a single data file:

```
python process/user/rg_pp/analysis_rg_data.py -f /path/to/AnalysisResults.root -c /path/to/my/config.yaml -o /path/to/my/outputdir
```


--------------------------------------------------------------------
Contact: James Mulligan (james.mulligan@berkeley.edu)
