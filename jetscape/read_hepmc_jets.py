import sys
import os
import fastjet as fj
import pyhepmc_ng
import tqdm

def main():
	input_file="$HOME/data/jetscape/test_out.hepmc"
	if len(sys.argv) > 1:
	    input_file = sys.argv[1]

	input_file = os.path.expandvars(input_file)

	print('[i] reading from:', input_file)
	# input = pyhepmc_ng.ReaderAsciiHepMC2(input_file)
	input = pyhepmc_ng.ReaderAscii(input_file)
	if input.failed():
		print ("[error] unable to read from {}".format(input_file))
		return
	nevents = 1000

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	jet_R0 = 0.4
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	jet_selector = fj.SelectorPtMin(0.0) & fj.SelectorPtMax(200.0) & fj.SelectorAbsEtaMax(3)

	all_jets = []
	event = pyhepmc_ng.GenEvent()
	pbar = tqdm.tqdm(range(nevents))
	while not input.failed():
		e = input.read_event(event)
		if input.failed():
			break
		fjparts = []
		for i,p in enumerate(event.particles):
			if p.status == 1 and not p.end_vertex:
				psj = fj.PseudoJet(p.momentum.px, p.momentum.py, p.momentum.pz, p.momentum.e)
				psj.set_user_index(i)
				fjparts.append(psj)
		jets = jet_selector(jet_def(fjparts))
		all_jets.append([ [j.pt(), j.eta()] for j in jets])
		pbar.update()
		for j in jets:
			hjetpt.Fill(j.perp())
		if pbar.n >= nevents:
			break
	pbar.close()
if __name__ == '__main__':
	main()