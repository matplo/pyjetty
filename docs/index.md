# Synopsis

[pyjetty](https://github.com/matplo/pyjetty) uses functionality implemented in [heppy](https://github.com/matplo/heppy) (public)

# A short story - to compile "all" what's needed

```
this_directory=$PWD
git clone git@github.com:matplo/heppy.git
cd ./heppy

# we will build in virtualenv
./scripts/pipenv_heppy.sh shell
pipenv install numpy

# now build commands here... (root takes sometime)
./external/fastjet/build.sh
./external/lhapdf6/build.sh
./external/hepmc/build.sh
./external/hepmc3/build.sh
./external/root/build.sh
./external/pythia8/build.sh
./cpptools/build.sh

# run a test (but before install deps)
pipenv install tqdm argparse 
./heppy/examples/pythia_gen_fastjet_lund_test.py

# now install pyjetty
cd ${this_directory}

git clone git@github.com:matplo/pyjetty.git
export HEPPY_DIR=${this_directory}/heppy
./cpptools/build.sh
```

**warning:** this was tested on Mojave (brew python3) and CENTOS-7 (with devtools-7); we will make it running for Catalina; we do not expect issues for other Linux-like OS'

## alternative with pre-installed packages

The installation on heppy will look for packages in the following directories defined by environment variables:
```
FASTJET_DIR
LHAPDF_DIR xor LHAPDF6_DIR
HEPMC_DIR xor HEPMC2_DIR
HEPMC3_DIR xor HEPMC2_DIR
ROOTSYS
PYTHIA8_DIR
```

**important note:** all those should use the same python3 (we do not support python 2 anymore)... 

## some other notable dependencies - python packages only

Most of those can be installed by simply using `pip install`
- `numpy`
- `yaml`
- `argparse` - argument parsing
- `tqdm` (Mateusz like's it - some others may be annoyed by it) - a progress bar
- `uproot` - ROOT I/O in pure ROOT (in principle some of the deps on the compiled ROOT could be mitigated with this package alone)
- `pyhepmc_ng` support for pythonic HEPMC

Less frequently used (but still)
- `math`
- `ctypes`
- `shutil`
- `array`
- `joblib` - a pickle I/O

Some third party - can ignore for the moment (MP uses it a lot):
- `dlist` from https://github.com/matplo/rootutils (some more deps thererin)

### install all one line

... except the third party rootutils

```
pip install numpy yaml argparse tqdm uproot pyhepmc_ng math ctypes shutil array joblib
```

# Contributing

Please fork and make a pull request.
Please let us know if you are using this code - we are very much open for collaboration - email us at heppy@lbl.gov
