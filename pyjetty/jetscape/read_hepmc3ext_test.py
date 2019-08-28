import sys
import os
import hepmc3ext

def main():
	input_file="$HOME/data/jetscape/test_out.hepmc"
	if len(sys.argv) > 1:
	    input_file = sys.argv[1]

	input_file = os.path.expandvars(input_file)

	print('[i] reading from:', input_file)
	hepmc3ext.test_loop(input_file)

if __name__ == '__main__':
	main()
