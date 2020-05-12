# Synopsis

[pyjetty](https://github.com/matplo/pyjetty) uses functionality implemented in [heppy](https://github.com/matplo/heppy) (public)

# A short story - to compile "all" what's needed

```
this_directory=$PWD
git clone git@github.com:matplo/heppy.git
cd ./heppy

# we will build in virtualenv
./scripts/pipenv_heppy.sh install numpy
./scripts/pipenv_heppy.sh --three shell
# now build commands here... (root takes sometime)
./external/fastjet/build.sh
./external/lhapdf6/build.sh
./external/hepmc2/build.sh
./external/hepmc3/build.sh
./external/root/build.sh
./external/pythia8/build.sh
./cpptools/build.sh
```

```
module use ${this_directory}/heppy/modules
# run a test (but before install deps)
pipenv install tqdm argparse 
${HEPPY_DIR}/heppy/examples/pythia_gen_fastjet_lund_test.py
```

```
# now install pyjetty
cd ${this_directory}

git clone git@github.com:matplo/pyjetty.git
export HEPPY_DIR=${this_directory}/heppy
cd pyjetty
./cpptools/build.sh
```

**warning:** this was tested on Mojave (brew python3) and CENTOS-7 (with devtools-7); we will make it running for Catalina; we do not expect issues for other Linux-like OS'

**note:** fastjet/build.sh now takes an optional --cgal flag that enables an attempt to compile with cgal (we've encountered trouble using CGAL header only installation with cmake on MacOSX - cgal is disabled by default - try the flag and see whether the configure step is successful...)

# using pyjetty (not a compile mode)

The build generated an envoronment module (http://modules.sourceforge.net) setting up things for normal (python coding / running) use so you can:

```
module use ${this_directory}/modules
module load pyjetty/1.0
```

and you can try a test code

```
$PYJETTY_DIR/pyjetty/examples/pythia_gen_fastjet_lund_test.py
```

## I like to have a few bash functions in my login scripts (.profile or .bashrc for example)

```
function pyjetty_load()
{
    export PS1="(pyj)$PS1"
    if [ -d $HOME/devel/pyjetty/modules ]; then
        module use $HOME/devel/pyjetty/modules
        module load pyjetty/1.0
    fi

}
export -f pyjetty_load
```

such than from the command line I can type

```
pyjetty_load
```

and the module is loaded.

Note, on CENTOS-7 we want devtools - so for example (insert before `module load`):

```
source scl_source enable devtoolset-7
```

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
- `pyyaml`
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
pip install numpy pyyaml argparse tqdm uproot pyhepmc_ng math ctypes shutil array joblib
```

# Contributing

Please fork and make a pull request.
Please let us know if you are using this code - we are very much open for collaboration - email us at heppy@lbl.gov
