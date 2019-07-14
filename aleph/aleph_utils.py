import pythia8
from mptools import pymptools as mpt
import os
from configobj import ConfigObj

def get_n_events(fname):
	ic = 0
	info_file = fname+'.info'
	cobj = ConfigObj(info_file)
	try:
		ic = int(cobj['nevents'])
	except KeyError:
		print('[i] counting events in', fname)
		reader = mpt.Reader(fname)
		while reader.read_next_event():
			ic += 1
		cobj['nevents']=ic
		cobj.write()
	print ('[i] info file', info_file, ic)
	return ic

#enum SIMPLEPWFLAG {ALEPH_CHARGED_TRACK, ALEPH_CHARGED_LEPTONS1, ALEPH_CHARGED_LEPTONS2, ALEPH_V0, ALEPH_PHOTON, ALEPH_NEUTRAL_HADRON}

class AlephWFLAG(object):
	ALEPH_CHARGED_TRACK		= 0
	ALEPH_CHARGED_LEPTONS1 	= 1
	ALEPH_CHARGED_LEPTONS2	= 2
	ALEPH_V0				= 3
	ALEPH_PHOTON			= 4
	ALEPH_NEUTRAL_HADRON	= 5
	def __init__(self, iflag):
		self.flag = iflag
