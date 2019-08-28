import sys
import os
import pyhepmc_ng
import tqdm

def main():
	input_file="$HOME/data/jetscape/test_out.hepmc"
	if len(sys.argv) > 1:
	    input_file = sys.argv[1]

	input_file = os.path.expandvars(input_file)

	print('[i] reading from:', input_file)
	input = pyhepmc_ng.ReaderAscii(input_file)
	if input.failed():
		print ("[error] unable to read from {}".format(input_file))
		return

	event = pyhepmc_ng.GenEvent()
	while not input.failed():
		e = input.read_event(event)
		if input.failed():
			break

if __name__ == '__main__':
	main()
