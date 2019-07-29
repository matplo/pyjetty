# pyjetty

- some work on HEP things in python (with quite a bit of help on c++ side)
- this uses CMake, SWIG
- Python3 required
- recommendation: work within pipenv (or virtual env)

# example python script

 - `cpptools/tests/pythia_gen_fj_lund_test.py`

# recommended build/setup

```
pipenv shell
source setup.sh --install
```

- to 'redo the install' of the externals (all) or individually
```
source ./external/setup_env.sh --install
```
- useful debuging option `--configure-only`

- to rebuild cpptools
```
source setup.sh --re
```
- useful debuging option `--configure-only`

OR

```
./cpptools/scripts/build_cpptools.sh
```
-useful options `--clean` (only build) or '--cleanall' (build and install)

## Notes: 
- this will download and install PYTHIA, HepMC2, HepMC3, LHAPDF6, FASTJET into the `external` subdirectory. This behavior can be controlled by `.pyjetty_config_external` file (sourced as a shell script) - you can control what version packages to use by building those libs yourself... (no or empty `.pyjetty_config_external` is fine)
- the `.pyjetty_config_external` in a local directory takes precedence (default is to take one from the downloaded/git directory)
- for some options `./scripts/build_cpptools.sh --help`

# running

- example 'workflow' (note no `--install`)

```
pipenv shell
source setup.sh
./cpptools/tests/pythia_gen_fj_lund_test.py
...
```

## alternative

- build PYTHIA [HepMC2, HepMC3, LHAPDF6], FASTJET and set appropriate environment variables (LHAPDF6 is completely optional; HEPMC2 also but pythiaext module will not be built);
- for PYTHIA and/or FASTJET `pythia8-config` and/or `fastjet-config` are expected to be accessible (in PATH)
- export $PYJETTY_SETUP_EXTERNAL to point to a shell script setting up the environment or set it to a value - for example: `export ${PYJETTY_SETUP_EXTERNAL}=mysetup.sh`

# add/extend c++ (swig) to python

- in the `cpptools/src` directory create your code/directory
- edit `cpptools/src/pyjetty.i`, `cpptools/src/CMakeLists.txt` as needed
- run `./cpptools/scripts/build_cpptools.sh`
