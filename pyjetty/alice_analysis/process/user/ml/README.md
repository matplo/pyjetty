# Running the code

## Setup on hiccup

Please note that in order to run the code on the hiccup cluster, all but the lightest code tests should use the slurm batch system. This means that every code you run should either:
1. Be called through `sbatch`, which will launch your jobs across the 8 hiccup nodes (see steps 0 and 1 below for examples)
2. Be run on an interactive node requested via slurm (this should always be done for steps 2-3 below):
   ```
   srun -N 1 -n 20 -t 24:00:00 -p std --pty bash
   ``` 
   which requests 1 full node (20 cores) for 24 hours in the std partition. You can choose the time, queue, and number of cores to suite your needs. When you’re   done with your session, just type `exit`.

## pp-AA (full events) from James

The events are located on hiccup at:
- JEWEL (vacuum + medium) with thermal background, full events: `/rstorage/ml/egml/data/files.txt`

The analysis workflow is as follows:

0. Generate JEWEL events: 
   ```
   cd ~/pyjetty/pyjetty/alice_analysis/generation/jewel/vac
   sbatch slurm_jewel_pp.sh
   ```
   and 
   ```
   cd ~/pyjetty/pyjetty/alice_analysis/generation/jewel/medium
   sbatch slurm_jewel_PbPb.sh
   ```
   Then create a filelist manually in `/rstorage/ml/egml/data`, to be used in the next step:
   ```
   find "$PWD" -name "*.root" >files.txt
   ...
   cat <pp>/files.txt >> <AA>/files.txt
   ```
   Note that we do not need to shuffle anything yet.
   
1. Compute Nsubjettiness arrays from input events, and write them to file, along with labels and four-vectors: 
   ```
   cd ~/pyjetty/pyjetty/alice_analysis/process/user/ml/slurm
   sbatch slurm_compute_nsubjettiness.sh
   ```
   You should update the skimmed filelist in `slurm_compute_nsubjettiness.sh` if necessary.
   
   This step use a common config file at `alice_analysis/config/ml/ppAA.yaml`.
   
   This writes output to `/rstorage/ml/egml/nsubjettiness/<job_id>/files.txt`
   
   Shuffle the file list: `shuf --output=files_randomized.txt files.txt`

2. Aggregate the results from each file's processing output
   ```
   cd ~/pyjetty/pyjetty/alice_analysis/process/user/ml
   python aggregate_nsubjettiness.py -o /rstorage/ml/egml/nsubjettiness/<process_job_id>
   ```
   The `-o` path should point to the directory containing `files.txt` from Step 2. This is the location that the output file, `nsubjettiness.h5`, will be written to. 
   
3. Fit model and make plots:
   ```
   cd ~/pyjetty/pyjetty/alice_analysis/analysis/user/ml
   python analyze_ppAA.py -c <config> -o <output_dir>
   ```
   The `-o` path should point to the directory containing `nsubjettiness.h5` from Step 3. This is the location that the output plots will be written to. 
   
   This step use a common config file at `alice_analysis/config/ml/ppAA.yaml`.

## Quark-gluon

The analysis consists of two steps:
1. Compute Nsubjettiness arrays from input events, and write them to file: 
   ```
   cd alice_analysis/process/user/ml
   python process_qg.py -c <config> -o <output_dir>
   ```
2. Read in file, and fit model:
   ```
   cd alice_analysis/analysis/user/ml
   python analyze_qg.py -c <config> -o <output_dir>
   ```

Both steps use a common config file at `alice_analysis/config/ml/qg.yaml`.

# Pull requests

To contribute code to pyjetty, you should make a pull request from your fork of the repo.

For example, working with a single branch:

1. First, create your fork of `matplo/pyjetty` by clicking the "Fork" tab in the upper right of the `matplo/pyjetty` webpage.

2. Next, go to your local pyjetty repository (doing `git clone` from the original `matplo/pyjetty`). We want to add a "remote" to your fork. To see your current remotes, you can type `git remote -v`. 

- To add a remote to your fork, type `git remote add myFork git@github.com:<github-username>/pyjetty.git` (or using https, if you use that)

3. Suppose you now have a commit on your local repository that you want to propose to add to the main `pyjetty` repository.

- First, we should update our local code to the latest changes in `matplo/pyjetty` and put our commit on top of that: `git pull --rebase`
- Next, commit to your fork: `git push myFork master`
- Then, go to your fork and create a pull request.
