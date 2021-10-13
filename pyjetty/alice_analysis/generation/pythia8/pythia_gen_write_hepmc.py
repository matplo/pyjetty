#!/usr/bin/env python

from __future__ import print_function

import fastjet
import fjcontrib
import fjext

import tqdm
import argparse
import os

import pythia8
import pythiaext

# experimental hepmc3
# import pythiahepmc3

from heppy.pythiautils import configuration as pyconf


def main():
    parser = argparse.ArgumentParser(description='pythia8 in python', prog=os.path.basename(__file__))
    parser.add_argument('-o', '--output_dir', help='output file location', default='.', type=str, required=True)
    pyconf.add_standard_pythia_args(parser)
    args = parser.parse_args()	

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    output_file = os.path.join(args.output_dir, 'pythia8.hepmc')

    pyhepmc2writer = pythiaext.Pythia8HepMC2Wrapper(output_file)
    # pyhepmc3writer = pythiahepmc3.Pythia8HepMC3Wrapper(output_file)

    mycfg = []
    pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
    if args.nev < 10:
        args.nev = 10
    for i in range(args.nev):
        if not pythia.next():
            continue
        if i%100 == 0:
            print(f'event {i}/{args.nev}')

        pyhepmc2writer.fillEvent(pythia)
        # pyhepmc3writer.fillEvent(pythia)

    pythia.stat()
    pythia.settings.writeFile(os.path.join(args.output_dir, 'modified_settings.hepmc'))
    print(f"[i] file written: {output_file}")

if __name__ == '__main__':
	main()