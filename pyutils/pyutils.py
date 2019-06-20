import pythia8


def create_and_init_pythia(config_strings=[]):
	pythia = pythia8.Pythia()
	for s in config_strings:
		pythia.readString(s)
	for extra_s in ["Next:numberShowEvent = 0", "Next:numberShowInfo = 0", "Next:numberShowProcess = 0"]:
		pythia.readString(extra_s)
	if pythia.init():
		return pythia
	return None
