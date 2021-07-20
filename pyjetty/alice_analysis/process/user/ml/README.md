# Running the code

## pp-AA (full events) from James

The events are located on hiccup at:
- JEWEL (vacuum + medium) with thermal background, full events: `/rstorage/ml/egml/data/files.txt`

The analysis workflow is as follows:

0. Generate JEWEL events: 
   ```
   cd pyjetty/alice_analysis/generation/jewel/vac
   sbatch slurm_jewel_pp.sh
   ```
   and 
   ```
   cd pyjetty/alice_analysis/generation/jewel/medium
   sbatch slurm_jewel_PbPb.sh
   ```
   Then create a filelist manually in `/rstorage/ml/egml/data`, to be used in the next step.
   
1. Compute Nsubjettiness arrays from input events, and write them to file, along with labels and four-vectors: 
   ```
   cd pyjetty/alice_analysis/process/user/ml/slurm
   sbatch slurm_compute_nsubjettiness.sh
   ```
   You should update the skimmed filelist in `slurm_compute_nsubjettiness.sh` if necessary.
   
   This step use a common config file at `alice_analysis/config/ml/ppAA.yaml`.
   
   This writes output to `/rstorage/ml/egml/nsubjettiness/<job_id>/files.txt`

2. Aggregate the results from each file's processing output
   ```
   cd pyjetty/alice_analysis/process/user/ml
   python aggregate_nsubjettiness.py -o /rstorage/ml/egml/nsubjettiness/<process_job_id>
   ```
   The `-o` path should point to the directory containing `files.txt` from Step 2. This is the location that the output file, `nsubjettiness.h5`, will be written to. 
   
3. Fit model and make plots:
   ```
   cd alice_analysis/analysis/user/ml
   python analyze_ppAA.py -c <config> -o <output_dir>
   ```
   The `-o` path should point to the directory containing `nsubjettiness.h5` from Step 3. This is the location that the output plots will be written to. 
   
   This step use a common config file at `alice_analysis/config/ml/ppAA.yaml`.

## pp-AA (full events) from Yue Shi

The events are located on hiccup at:
- PYTHIA/JEWEL + HYDJET, full events: `/rstorage/ml/egml/data/v3/files_hydjet.txt`

The analysis workflow is as follows:

1. Skim the events into numpy arrays of four-vectors and labels
   ```
   cd pyjetty/alice_analysis/process/user/ml/slurm
   sbatch slurm_skim.sh
   ```
   This writes the output to `/rstorage/ml/egml/skim/<job_id>/files.txt` (for both PYTHIA and JEWEL, and both hard event + background).
   
   This step likely does not have to be repeated, unless we get new event samples.

2. Compute Nsubjettiness arrays from input events, and write them to file, along with labels and four-vectors: 
   ```
   cd pyjetty/alice_analysis/process/user/ml/slurm
   sbatch slurm_compute_nsubjettiness.sh
   ```
   You should update the skimmed filelist in `slurm_compute_nsubjettiness.sh` if necessary (if step 1 is repeated).
   
   This step use a common config file at `alice_analysis/config/ml/ppAA.yaml`.
   
   This writes output to `/rstorage/ml/egml/nsubjettiness/<job_id>/files.txt`

3. Aggregate the results from each file's processing output
   ```
   cd pyjetty/alice_analysis/process/user/ml
   python aggregate_nsubjettiness.py -o /rstorage/ml/egml/nsubjettiness/<process_job_id>
   ```
   The `-o` path should point to the directory containing `files.txt` from Step 2. This is the location that the output file, `nsubjettiness.h5`, will be written to. 
   
4. Fit model and make plots:
   ```
   cd alice_analysis/analysis/user/ml
   python analyze_ppAA.py -c <config> -o <output_dir>
   ```
   The `-o` path should point to the directory containing `nsubjettiness.h5` from Step 3. This is the location that the output plots will be written to. 
   
   This step use a common config file at `alice_analysis/config/ml/ppAA.yaml`.
   
## pp-AA (jets only)

There are also earlier versions of events where only the jets are included:
- PYTHIA/JEWEL + HYDJET, jets: `/rstorage/ml/egml/data/v2/with_hydjet/files.txt`
- PYTHIA/JEWEL (no UE), jets: `/rstorage/ml/egml/data/v2/no_ue/files.txt`
- (The initial events, now superceded, are located at `/rstorage/ml/egml/data/v1/files.txt`)

You would need to check out an earlier version of the code (c. April 2021) to use this. I dont' foresee this will be used, but for completeness I keep the documentation here:

1. Skim the events into numpy arrays of four-vectors and labels
   ```
   cd pyjetty/alice_analysis/process/user/ml/slurm
   sbatch slurm_skim.sh                                # (for PYTHIA/JEWEL + HYDJET)
   sbatch slurm_skim_no_ue.sh                          # (for no UE)
   ```
   This writes the output to `/rstorage/ml/egml/skim/<job_id>/files.txt` (for both PYTHIA and JEWEL).

2. Compute Nsubjettiness arrays from input events, and write them to file, along with labels and four-vectors: 
   ```
   cd pyjetty/alice_analysis/process/user/ml/slurm
   sbatch slurm_compute_nsubjettiness.sh               # (for PYTHIA/JEWEL + HYDJET)
   sbatch slurm_compute_nsubjettiness_no_ue.sh         # (for no UE)
   ```
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
