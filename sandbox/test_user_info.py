import fastjet as fj
import pythia8
from recursivetools import pyrecursivetools as rt
from pythiafjtools import pypythiafjtools as pyfj
#from tqdm import tnrange, tqdm_notebook
from tqdm import tqdm


def create_and_init_pythia(config_strings=[]):
	pythia = pythia8.Pythia()
	for s in config_strings:
		pythia.readString(s)
	for extra_s in ["Next:numberShowEvent = 0", "Next:numberShowInfo = 0", "Next:numberShowProcess = 0"]:
		pythia.readString(extra_s)
	if pythia.init():
		return pythia
	return None


def print_jets(jets):
	for ijet in range(len(jets)):
		jet = jets[ijet]
		constituents = jet.constituents()
		print("{0:>5s} {1:>10s} {2:>10s} {3:>10s} {4:>12s}".format(
			"jet #", "pt", "rap", "phi", "N particles"))
		print("{0:5d} {1:10.3f} {2:10.4f} {3:10.4f} {4:5d}".format(
			ijet, jet.pt(), jet.rap(), jet.phi(), len(constituents)))
		for c in fj.sorted_by_pt(constituents):
			_p = pyfj.getPythia8Particle(c)
			if _p:
				print(" - {0:>10s} id={1:5d} status={2:5d} pT={3:10.3f}".format(_p.name(), _p.id(), _p.status(), _p.pT()))


def main():
	fj.ClusterSequence.print_banner()
	print()

	sconfig_pythia = ["Beams:eCM = 8000.", "HardQCD:all = on", "PhaseSpace:pTHatMin = 20."]
	pythia = create_and_init_pythia(sconfig_pythia)
	jet_R0 = 0.4
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	jet_selector = fj.SelectorPtMin(20.0) & fj.SelectorPtMax(40.0) & fj.SelectorAbsEtaMax(1)
	# sd = rt.SoftDrop(0, 0.1, 1.0)

	all_jets = []
	for iEvent in tqdm(range(1000), 'event'):
		if not pythia.next():
			continue
		parts = pyfj.vectorize(pythia, True, -1, 1, True)
		jets = jet_selector(jet_def(parts))
		all_jets.extend(jets)

	print_jets(all_jets)

if __name__ == '__main__':
	main()
