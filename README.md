# pyjetty

- meant as an extension of [https://github.com/matplo/heppy](https://github.com/matplo/heppy)

# example python script

 - `cpptools/tests/pythia_gen_fj_lund_test.py`

# recommended build/setup

```
pipenv shell
./scripts/setup.sh --buildext
module use <somedir>/heppy/modules
module load heppy/main_python
module use $PWD/modules
module load pyjetty/main_python
./cpptools/tests/pythia_gen_fj_lund_test.py
```

- to wipe out and rebuild everything
```
./scripts/cleanup.sh
./scripts/setup.sh --rebuild
```

- useful debuging option `--configure-only`

- to rebuild cpptools only
```
./scripts/setup.sh --rebuild
```

- one can also use ./cpptools/scripts/build_cpptools.sh directly with [--rebuild] [--install] [--clean] [--cleanall]
- for some options `./cpptools/scripts/build_cpptools.sh --help`

# modules

- ./modules/pyjetty contains modules - use the one with 'main' to load everything


# add/extend c++ (swig) to python

- in the `cpptools/src` directory create your code/directory
- edit `cpptools/src/pyjetty.i`, `cpptools/src/CMakeLists.txt` as needed
- run `./cpptools/scripts/build_cpptools.sh`
