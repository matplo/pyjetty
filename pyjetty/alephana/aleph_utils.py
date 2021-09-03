import os
from configobj import ConfigObj
import aleph

def get_n_events(fname):
	ic = 0
	info_file = fname+'.info'
	cobj = ConfigObj(info_file)
	try:
		ic = int(cobj['nevents'])
	except KeyError:
		print('[i] counting events in', fname)
		reader = aleph.Reader(fname)
		while reader.read_next_event():
			ic += 1
		cobj['nevents']=ic
		cobj.write()
	print ('[i] info file', info_file, ic)
	return ic

def get_n_events_gzip(fname):
	import gzip
	ic = 0
	info_file = fname+'.info'
	cobj = ConfigObj(info_file)
	try:
		ic = int(cobj['nevents'])
	except KeyError:
		print('[i] counting events in', fname)
		with gzip.open(fname) as f:
			data = f.readlines()
			for l in data:
				lstr = l.decode("utf-8")
				if 'END_EVENT' in lstr:
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
