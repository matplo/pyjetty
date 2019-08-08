import pythia8


def create_and_init_pythia(config_strings=[]):
	pythia = pythia8.Pythia()
	extra_s = ["Next:numberShowEvent = 0", "Next:numberShowInfo = 0", "Next:numberShowProcess = 0"]
	config_strings.extend(extra_s)
	for s in config_strings:
		pythia.readString(s)
	if pythia.init():
		print ('[i] pythia initialized with', config_strings)
		return pythia
	return None


def add_standard_pythia_args(parser):
	parser.add_argument('--ecms', help='low or high sqrt(s) GeV', default='low', type=str)
	parser.add_argument('--ecm', help='sqrt(s) GeV', default=13000, type=float)
	parser.add_argument('--pthatmin', help='minimum hat{pT}', default=-1, type=float)
	parser.add_argument('--biaspow', help='power of the bias (hard)', default=4, type=float)
	parser.add_argument('--biasref', help='reference pT for the bias', default='50', type=float)
	parser.add_argument('--noue', help="no underlying event - equivalend to no ISR and MPIs set to off", default=False, action='store_true')
	parser.add_argument('--noISR', help="ISR set to off", default=False, action='store_true')
	parser.add_argument('--noMPI', help="MPIs set to off", default=False, action='store_true')
	parser.add_argument('--hard', help="enable hardQCD (ON if no other process selected)", default=False, action='store_true')
	parser.add_argument('--charm', help="enable hardccbar", default=False, action='store_true')
	parser.add_argument('--beauty', help="enable hardbbbar", default=False, action='store_true')
	parser.add_argument('--photon', help="enable prompt photon production",  default=False, action='store_true')
	parser.add_argument('--nev', help='number of events', default=1, type=int)
	parser.add_argument('-n', '--nevents', help='number of events', default=1000, type=int)


def pythia_config_from_args(args):
	sconfig_pythia = []
	procsel = 0
	if args.photon:
		sconfig_pythia.append("PromptPhoton:all = on")
		procsel += 1
	if args.charm:
		sconfig_pythia.append("HardQCD:hardccbar = on")
		procsel += 1
	if args.beauty:
		sconfig_pythia.append("HardQCD:hardbbbar = on")
		procsel += 1
	if args.hard:
		sconfig_pythia.append("HardQCD:all = on")
		procsel += 1
	if procsel == 0:
		sconfig_pythia.append("HardQCD:all = on")
	if args.pthatmin < 0:
		_extra = [	"PhaseSpace:bias2Selection=on",
					"PhaseSpace:bias2SelectionPow={}".format(args.biaspow),
					"PhaseSpace:bias2SelectionRef={}".format(args.biasref)]
		sconfig_pythia.extend(_extra)
	else:
		sconfig_pythia.append("PhaseSpace:pTHatMin = {}".format(args.pthatmin))
	if args.noue:
		sconfig_pythia.append("PartonLevel:ISR = off")
		sconfig_pythia.append("PartonLevel:MPI = off")
	if args.noISR:
		sconfig_pythia.append("PartonLevel:ISR = off")
	if args.noMPI:
		sconfig_pythia.append("PartonLevel:MPI = off")
	if args.ecm:
		sconfig_pythia.append("Beams:eCM = {}".format(args.ecm))
	return sconfig_pythia

def create_and_init_pythia_from_args(args, user_args=[]):
	_args = []
	from_args = pythia_config_from_args(args)
	_args.extend(from_args)
	_args.extend(user_args)
	pythia = create_and_init_pythia(_args)
	return pythia
