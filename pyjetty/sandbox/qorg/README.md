# generate PYTHIA8 hardQCD events 

- just run ./gen_ml_sample.sh --nev=<number of events>

## script runs

````

  # quarks and gluons
  ${THISD}/pythia_gen.py                     --py-noue --py-pthatmin 100. --py-seed ${seed} --ml --nev ${nev}

  # only ff->gg
	${THISD}/pythia_gen.py --py-hardQCDgluons --py-noue --py-pthatmin 100. --py-seed ${seed} --ml --nev ${nev}

  # only ff->qq or ff->qqbar
	${THISD}/pythia_gen.py --py-hardQCDquarks --py-noue --py-pthatmin 100. --py-seed ${seed} --ml --nev ${nev}

```

## notes

- by default the split is into 1k events - change with `--nperfile <number>` option to the python script
- MPI's and ISR are OFF
- in general as usual use `pythia_gen.py -h` for options...
- after the run the pythia settings are saved to a `.cmnd` file - see examples below

### pythia settings

- any q or glue
```
! List of all modified PYTHIA 8.244 settings.
HardQCD:all = on
Next:numberCount = 0
Next:numberShowEvent = 0
Next:numberShowInfo = 0
Next:numberShowProcess = 0
PartonLevel:ISR = off
PartonLevel:MPI = off
PhaseSpace:pTHatMin = 100.00000
Random:seed = 123456
Random:setSeed = on
```

- only gluons 
```
! List of all modified PYTHIA 8.244 settings.
HardQCD:gg2gg = on
HardQCD:qqbar2gg = on
Next:numberCount = 0
Next:numberShowEvent = 0
Next:numberShowInfo = 0
Next:numberShowProcess = 0
PartonLevel:ISR = off
PartonLevel:MPI = off
PhaseSpace:pTHatMin = 100.00000
Random:seed = 123456
Random:setSeed = on
```

- only quarks 
```
! List of all modified PYTHIA 8.244 settings.
HardQCD:gg2qqbar = on
HardQCD:hardbbbar = on
HardQCD:hardccbar = on
HardQCD:qq2qq = on
HardQCD:qqbar2qqbarNew = on
Next:numberCount = 0
Next:numberShowEvent = 0
Next:numberShowInfo = 0
Next:numberShowProcess = 0
PartonLevel:ISR = off
PartonLevel:MPI = off
PhaseSpace:pTHatMin = 100.00000
Random:seed = 123456
Random:setSeed = on
```