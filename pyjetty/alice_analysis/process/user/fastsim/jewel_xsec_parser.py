''' jewel_xsec_parser.py
Recursively iterate through jewel.log files, find cross sections,
sum event numbers, and create scaleFactors.yaml file to use for
pT-hat bin scaling.

Author: Ezra Lesser (elesser@berkeley.edu), Winter 2022
'''

from __future__ import division, print_function
from os import path
import numpy as np

######################################################################
### USER-DEFINED CONSTANTS HERE:

# Input directory where the pT-hat bins are stored
IN_DIR = "/rstorage/generators/jewel_alice/tree_gen/851894/"

# Output directory where the scaleFactors.yaml files should be created
OUT_DIR = "./"

# Number of pT-hat bins
N_PTHAT = 20

# Number of statistics bins
N_STATS = 400

# Multiply all scale factors by a constant so that first bin is 10
# (Useful to prevent very small scale factors)
scale_to_first_bin = True


######################################################################
######################################################################

print("Calculating the scale factors per each pT-hat bin...")

scale_factors = np.zeros(N_PTHAT)

# Loop on pT-hat bins (1-N_PTHAT)
for i in range(1, N_PTHAT+1):
  print("pT-hat bin %i..." % i, end='\r')
  N_counts = 0
  sigma = 0
  # Loop over the different statistics for each bin
  for j in range(1, N_STATS+1):
    filename = path.join(IN_DIR, str(i), str(j), "jewel.log")
    with open(filename, 'r') as f:
      for line in f:
        if line[1:14] == "cross section":
          words = [word for word in line.strip().split(' ') if len(word)]
          sigma += float(words[-2])
        elif line[1:21] == "sum of event weights":
          words = [word for word in line.strip().split(' ') if len(word)]
          N_counts += float(words[-1])
  # Average sigma per file
  sigma /= N_STATS
  # Save sigma / N as the scale factor for this pT-hat bin
  scale_factors[i-1] = sigma / N_counts

if scale_to_first_bin:
  val = scale_factors[0] / 10.
  scale_factors = [sf/val for sf in scale_factors]

print("Creating scaleFactors.yaml file in %s" % OUT_DIR)

with open(path.join(OUT_DIR, "scaleFactors.yaml"), "w") as f:
  for i in range(1, N_PTHAT+1):
    f.write("%i: %.15e\n" % (i, scale_factors[i-1]))

print("Done!")
