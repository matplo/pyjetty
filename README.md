# pyjetty

- meant as an extension of [https://github.com/matplo/heppy](https://github.com/matplo/heppy)

# quick start with heppy in a docker

- an example (you probably want to start in a different directory)...
```
$ cd /tmp
$ git clone https://github.com/matplo/heppy.git
$ ./heppy/docker/heppy/run.sh
$ cd /host/tmp
$ git clone https://github.com/matplo/pyjetty.git
$ heppy_load
$ ./pyjetty/cpptools/build.sh --tenngen --tglaubermc
```

- run an example within the docker `./pyjetty/pyjetty/examples/pythia_gen_fastjet_lund_test.py`

# recommended build/setup from source

 - install heppy and `export HEPPY_DIR=<path to heppy>`
 - execute `$HEPPY_DIR/scripts/pipenv_heppy.sh run <where pyjetty>/cpptools/build.sh` - that's it...

 - alternatively load the heppy module (from in heppy/modules)

# modules

- modules/pyjetty contains an env module

# example python script

```
module use <where pyjetty>/modules
module load pyjetty/1.0
pyjettypython $PYJETTY_DIR/cpptools/tests/pythia_gen_fj_lund_test.py
```
notes/tips: use a new shell after building; no need to load or set heppy the pyjetty module should take care of things

# add/extend c++ (swig) to python

- in the `cpptools/src` directory create your code/directory
- edit `cpptools/src/pyjetty.i`, `cpptools/src/CMakeLists.txt` as needed
- run `./cpptools/scripts/build_cpptools.sh`

# more functionality?

- TennGen background generator (via https://github.com/matplo/TennGen)
-- install with `$HEPPY_DIR/scripts/pipenv_heppy.sh run <where pyjetty>/cpptools/build.sh --tenngen`
- TGlauberMC - a Glauber MC implementation (via https://github.com/matplo/TGlauberMC)
-- install with `$HEPPY_DIR/scripts/pipenv_heppy.sh run <where pyjetty>/cpptools/build.sh --tglaubermc`

Note: to install both use both flags: `--tenngen --tglaubermc`

# contributing

Please fork and make a pull request.
Please let us know if you are using this code - we are very much open for collaboration - email us at heppy@lbl.gov
