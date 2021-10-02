import argparse
import sys

def get_args_from_settings(ssettings):
    sys.argv=[' '] + ssettings.split()
    parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly')
    pyconf.add_standard_pythia_args(parser)
    parser.add_argument('--output', default="test_ang_ue.root", type=str)
    parser.add_argument('--user-seed', help='pythia seed', default=1111, type=int)
    args = parser.parse_args()
    return args

from heppy.pythiautils import configuration as pyconf

def pythia_init_from_string(ssettings, mycfg = []):
    # example: ssettings = "--py-ecm 5000 --user-seed=100000 --nev 1000"
    mycfg = []
    args = get_args_from_settings(ssettings)
    pythia8 = pyconf.create_and_init_pythia_from_args(args, mycfg)
    return pythia8

