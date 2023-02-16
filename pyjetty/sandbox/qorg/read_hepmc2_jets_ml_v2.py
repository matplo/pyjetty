from __future__ import print_function
from ast import arg, parse
from flavor_tagger import FlavorTaggerUtil, FT_print_psj, FT_psj_to_str

import fastjet as fj
import fjcontrib
import fjext

import tqdm
import argparse
import os
import sys

import pythia8
import pythiaext
import pythiafjext
# pythia8+hepmc3 writing
# import pythiahepmc3

from heppy.pythiautils import configuration as pyconf

import importlib
root_avail = importlib.util.find_spec("ROOT")
if root_avail is not None:
	import ROOT

import pyhepmc_ng
import threading
import multiprocessing

def tau21_R(j, R, algo=fj.antikt_algorithm):
	fj.ClusterSequence.print_banner()
	jet_def = fj.JetDefinition ( algo, R )

	_j = fj.sorted_by_pt(jet_def(j.constituents()))[0]

	beta = 1.0
	nSub1_beta1 = fjcontrib.Nsubjettiness(1,   fjcontrib.OnePass_WTA_KT_Axes(), fjcontrib.UnnormalizedMeasure(beta))
	nSub2_beta1 = fjcontrib.Nsubjettiness(2,   fjcontrib.OnePass_WTA_KT_Axes(), fjcontrib.UnnormalizedMeasure(beta))

	tau1 = nSub1_beta1.result(_j)
	tau2 = nSub2_beta1.result(_j)
	tau_ratio = -1
	if tau1 != 0:
		tau_ratio = tau2/tau1
	return tau_ratio

def tau21(j):
	# fj.ClusterSequence.print_banner()
	# jet_def = fj.JetDefinition ( algo, R )
	beta = 1.0
	nSub1_beta1 = fjcontrib.Nsubjettiness(1,   fjcontrib.OnePass_WTA_KT_Axes(), fjcontrib.UnnormalizedMeasure(beta))
	nSub2_beta1 = fjcontrib.Nsubjettiness(2,   fjcontrib.OnePass_WTA_KT_Axes(), fjcontrib.UnnormalizedMeasure(beta))

	tau1 = nSub1_beta1.result(j)
	tau2 = nSub2_beta1.result(j)
	tau_ratio = -1
	if tau1 != 0:
		tau_ratio = tau2/tau1
	return tau_ratio


from pyjetty.mputils.generic_object import GenericObject

class UniqueRunNumber(GenericObject):
	def __init__(self, **kwargs):
		super(UniqueRunNumber, self).__init__(**kwargs)
		try:
			self.current_thread_id = threading.currentThread().native_id
		except:
			self.current_thread_id = threading.currentThread().ident
		self.fname = os.path.join(os.getcwd(), '.unique_run_numbers')

	def unique(self):
		run_number = None
		not_unique = self.load_run_numbers()
		if self.slurm:
			run_number = os.getenv('SLURM_JOB_ID')
		if run_number is None:
			run_number = self.current_thread_id			
		while run_number in not_unique:
			run_number += 1
		self.write_run_number(run_number)
		return int(run_number)

	def load_run_numbers(self):
		# load numbers from a file - one file per pid
		ids = []
		try:
			with open(self.fname, 'r') as f:
				ids = [int(l.strip()) for l in f.readlines()]
		except:
			if not os.path.exists(self.fname):
				with open(self.fname, 'w') as f:
					pass
		return ids

	def write_run_number(self, run_number):
		with open(self.fname, 'a') as f:
			f.writelines(['{}\n'.format(run_number)])

class JetFileOutput(GenericObject):
	def __init__(self, **kwargs):
		super(JetFileOutput, self).__init__(**kwargs)
		if self.args:
			if type(self.args) == dict:
				self.configure_from_dict(self.args, ignore_none=True)
			else:
				self.configure_from_dict(self.args.__dict__)
		if self.debug:
			print(self.args)
		self.sname = []
		self.args_d = self.args.__dict__
		if self.args_d['output']:
			self.sname = [self.args_d['output']]
		else:
			for opt in self.args_d:
				if self.args_d[opt]:
					self.sname.append(opt)
					if type(self.args_d[opt]) in [float, int]:
						self.sname.append(str(self.args_d[opt]))
		self.sname.append('.cmnd')
		self.cmnd_file_name = '_'.join(self.sname)
		self.nfile = 0

		self.jet_def_ca = fj.JetDefinition(fj.cambridge_algorithm, fj.JetDefinition.max_allowable_R)
		self.mDT = None
		try:
			if self.args_d['mDT']:
				self.mDT = fjcontrib.ModifiedMassDropTagger(self.args_d['mDTzcut'])
		except:
			pass
		self.urn = UniqueRunNumber(args=self.args)


	def write_file_partons(self, jets, pythia):
		root_file_name = '_'.join(self.sname).replace('.cmnd', '_{}.root'.format(self.nfile))
		root_file = ROOT.TFile(root_file_name, 'recreate')
		root_ntuple_parts = ROOT.TNtuple('tree_Particle_gen', 'particles from PYTHIA8 - jets', 'run_number:ev_id:ParticlePt:ParticleEta:ParticlePhi:ParticlePID:ParticleMass')
		root_ntuple_ev = ROOT.TNtuple('tree_Event_gen', 'event info from PYTHIA8 - jets', 'run_number:ev_id:xsec:code:partonID')
		root_ntuple_jet = ROOT.TNtuple('tree_jet_gen', 'jet kinematics','run_number:ev_id:jetPt:jetEta:jetPhi:jetMass:ParticlePt:ParticleEta:ParticlePhi:ParticlePID:ParticleMass:dR')
		if pythia:
			run_number = self.args_d['py_seed'] + self.nfile
		else:
			run_number = self.nfile
		if pythia:
			if self.args_d['py_seed'] < 0:
				run_number = self.nfile
		run_number = self.urn.unique()
		for i, _j in enumerate(jets):
			ev_number_stream = i
			j = _j[0]
			j_flavor = _j[1]
			j_flavor_psj = _j[2]
			jetDR = j_flavor_psj.delta_R(j)
			if pythia:
				root_ntuple_ev.Fill(run_number, ev_number_stream, pythia.info.sigmaGen(), pythia.info.code(), j_flavor)
				root_ntuple_jet.Fill(run_number, ev_number_stream, j.pt(), j.eta(), j.phi(), j.m(), j_flavor_psj.pt(), j_flavor_psj.eta(), j_flavor_psj.phi(), pythiafjext.getPythia8Particle(j_flavor_psj).id(), j_flavor_psj.m(), jetDR)
			else:
				root_ntuple_ev.Fill(run_number, ev_number_stream, 0, 0, 0)
				root_ntuple_jet.Fill(run_number, ev_number_stream, j.pt(), j.eta(), j.phi(), j.m(), j_flavor_psj.pt(), j_flavor_psj.eta(), j_flavor_psj.phi(), 0, j_flavor_psj.m(), jetDR)
			for p in j.constituents():
				if pythia:
					if pythiafjext.getPythia8Particle(p).isParton() is False:
						root_ntuple_parts.Fill(run_number, ev_number_stream, p.pt(), p.eta(), p.phi(), pythiafjext.getPythia8Particle(p).id(), p.m())
				else:
					# note: for herwig the user index stores the PID
					root_ntuple_parts.Fill(run_number, ev_number_stream, p.pt(), p.eta(), p.phi(), p.user_index(), p.m())
		self.nfile += 1
		root_file.Write()
		root_file.Close()
		# print('[i] writing', root_file.GetName(), 'jets:', len(jets))

	def write_file_Zjets(self, jets, pythia):
		root_file_name = '_'.join(self.sname).replace('.cmnd', '_{}.root'.format(self.nfile))
		root_file = ROOT.TFile(root_file_name, 'recreate')
		root_ntuple_parts = ROOT.TNtuple('tree_Particle_gen', 'particles from PYTHIA8 - jets', 'run_number:ev_id:ParticlePt:ParticleEta:ParticlePhi:ParticlePID:ParticleMass')
		root_ntuple_parts_mDT = ROOT.TNtuple('tree_Particle_gen_mDT{}'.format(self.args_d['mDTzcut']), 'particles from PYTHIA8 - jets', 'run_number:ev_id:ParticlePt:ParticleEta:ParticlePhi:ParticlePID:ParticleMass')
		root_ntuple_ev = ROOT.TNtuple('tree_Event_gen', 'event info from PYTHIA8 - jets', 'run_number:ev_id:xsec:code')
		root_ntuple_zjet = ROOT.TNtuple('tree_Zjet_gen', 'Zjet kinematics','run_number:ev_id:ZjetPt:ZjetEta:ZjetPhi:ZjetMass:ParticlePt:ParticleEta:ParticlePhi:ParticlePID:ParticleMass:dR:tau21:tau21R')
		if pythia:
			run_number = self.args_d['py_seed'] + self.nfile
		else:
			run_number = self.nfile
		if pythia:
			if self.args_d['py_seed'] < 0:
				run_number = self.nfile
		run_number = self.urn.unique()
		for i, _j in enumerate(jets):
			ev_number_stream = i
			Zjet = _j[0]
			psjZ = _j[1]
			ZjetDR = psjZ.delta_R(Zjet)
			if pythia:
				root_ntuple_ev.Fill(run_number, ev_number_stream, pythia.info.sigmaGen(), pythia.info.code())
				root_ntuple_zjet.Fill(run_number, ev_number_stream, Zjet.pt(), Zjet.eta(), Zjet.phi(), Zjet.m(), psjZ.pt(), psjZ.eta(), psjZ.phi(), pythiafjext.getPythia8Particle(psjZ).id(), psjZ.m(), ZjetDR, tau21(Zjet), tau21_R(Zjet, self.args.jet_R))
			else:
				root_ntuple_ev.Fill(run_number, ev_number_stream, 0, 0)
				root_ntuple_zjet.Fill(run_number, ev_number_stream, Zjet.pt(), Zjet.eta(), Zjet.phi(), Zjet.m(), psjZ.pt(), psjZ.eta(), psjZ.phi(), 23, psjZ.m(), ZjetDR, tau21(Zjet), tau21_R(Zjet, self.args.jet_R))

			for p in Zjet.constituents():
				if pythia:
					root_ntuple_parts.Fill(run_number, ev_number_stream, p.pt(), p.eta(), p.phi(), pythiafjext.getPythia8Particle(p).id(), p.m())
				else:
					# note: for herwig the user index stores the PID
					root_ntuple_parts.Fill(run_number, ev_number_stream, p.pt(), p.eta(), p.phi(), p.user_index(), p.m())

			# now mDT jets
			if self.mDT:
				ca_jets = fj.sorted_by_pt(self.jet_def_ca(Zjet.constituents()))
				if (len(ca_jets) > 1):
					print('[w] CA reclustering returns more than 1 subjet (?)')
				else:
					# Use modified mass-drop tagger to clean up jet.
					taggedJet = self.mDT.result(ca_jets[0])
					if taggedJet.has_constituents():
						for p in taggedJet.constituents():
							if pythia:
								root_ntuple_parts_mDT.Fill(run_number, ev_number_stream, p.pt(), p.eta(), p.phi(), pythiafjext.getPythia8Particle(p).id(), p.m())
							else:
								# note: for herwig the user index stores the PID
								root_ntuple_parts_mDT.Fill(run_number, ev_number_stream, p.pt(), p.eta(), p.phi(), p.user_index(), p.m())
		self.nfile += 1
		root_file.Write()
		root_file.Close()

def print_jet(j):
	print(f'[i] jet pt={j.perp()} m={j.m()} nc={len(j.constituents())}')
	for c in j.constituents():
		print(f'   pt={c.perp()} phi={c.phi()} m={c.m()} {c.user_index()}')

def count_events(fname):
	count = None
	try:
		with open(fname+'.count', 'r') as f:
			_s = f.readlines()[-1]
		count = int(_s)
	except:
		pass
	if count:
		return count
	input_hepmc = pyhepmc_ng.ReaderAsciiHepMC2(fname)
	if input_hepmc.failed():
		print ("[error] unable to read from {}".format(fname))
		sys.exit(1)
	event_hepmc = pyhepmc_ng.GenEvent()
	pbar = tqdm.tqdm(desc='counting events... {}'.format(fname), leave=False)
	while not input_hepmc.failed():
		_ = input_hepmc.read_event(event_hepmc)
		if input_hepmc.failed():
			pbar.update(0)
			break
		pbar.update()
	n = pbar.n
	pbar.close()
	# print('\n[i] number of events', n, 'in', fname)
	try:
		with open(fname+'.count', 'w') as f:
			f.writelines(['{}'.format(n)])
	except:
		pass
	return n


def analyze_file(fname, args):
	input_hepmc = pyhepmc_ng.ReaderAsciiHepMC2(fname)
	if input_hepmc.failed():
		print ("[error] unable to read from {}".format(fname))
		sys.exit(1)

	nev_file = count_events(fname)
	nev_nperfile = args.nperfile
	if args.nperfile >= nev_file:
		nev_nperfile = nev_file

	part_selector = fj.SelectorAbsRapMax(4.)
	part_selector_p = fj.SelectorAbsRapMax(4.)
	jet_selector = fj.SelectorPtMin(args.jet_ptmin) * fj.SelectorPtMax(args.jet_ptmax) * fj.SelectorAbsRapMax(2.5)

	z_selector = fj.SelectorPtMin(args.Z_ptmin) * fj.SelectorPtMax(args.Z_ptmax) * fj.SelectorAbsRapMax(2.5)

	# jet_def_akt = fj.JetDefinition ( fj.antikt_algorithm, args.ZjetR)
	jet_def_akt = fj.JetDefinition ( fj.antikt_algorithm, args.jet_R)

	fout = JetFileOutput(args=args)

	run_number = 1 # this number should change
	urn = UniqueRunNumber()
	run_number = urn.unique()

	event_number = 0
	n_total_jets_accepted = 0
	jets_accepted = []

	pbar_ev = tqdm.tqdm(range(nev_file), desc='events {}'.format(fname))
	# pbar = tqdm.tqdm(range(args.nev), desc='accepted jets')
	n_jets_to_accept = args.nev
	if n_jets_to_accept >= nev_file:
		n_jets_to_accept = nev_file
	pbar = tqdm.tqdm(range(n_jets_to_accept), desc='accepted jets {}'.format(fname))

	event_hepmc = pyhepmc_ng.GenEvent()
	count_Zboson = 0
	while not input_hepmc.failed():
		ev = input_hepmc.read_event(event_hepmc)
		if input_hepmc.failed() or n_total_jets_accepted > args.nev:
			break

		pbar_ev.update()

		event_number = event_number + 1

		zs = []
		parts_pythia_h = fj.vectorPJ()
		for i, p in enumerate(event_hepmc.particles):
			# if p.id == 23: # BUG!!!
			if p.pid == 23:
				# print(p)
				# print(p.momentum.perp())
				count_Zboson += 1
				psj = fj.PseudoJet(p.momentum.px, p.momentum.py, p.momentum.pz, p.momentum.e)
				if len(z_selector([psj])) > 0:
					psj.set_user_index(i)
					zs.append(psj)
			if p.status == 1 and not p.end_vertex:
				psj = fj.PseudoJet(p.momentum.px, p.momentum.py, p.momentum.pz, p.momentum.e)
				# psj.set_user_index(abs(p.pid))
				psj.set_user_index(p.pid)
				parts_pythia_h.append(psj)

		parts_pythia_h = part_selector(parts_pythia_h)
		# print('[i] number of particles selected', len(parts_pythia_h))
		if len(parts_pythia_h) < 1 :
				continue

		jets = fj.sorted_by_pt(jet_def_akt(parts_pythia_h))
		jets_selected = jet_selector(jets)

		if len(jets_selected) < 1:
			continue

		# needs a fix from here on
		n_total_jets_accepted += 1

		# quark or gluon tagging
		if args.Zjet is False:
			if len(jets_selected) < 1:
				continue
			n_jets_accepted = 0
			for j in jets_selected:
				jets_accepted.append([j, 0, j])
				n_jets_accepted += 1

			if n_jets_accepted < 1:
				continue

			if len(jets_accepted) >= nev_nperfile:
				fout.write_file_partons(jets_accepted, pythia=None)
				jets_accepted.clear()

			n_total_jets_accepted += n_jets_accepted
			pbar.update(n_jets_accepted)
			continue

		# Zjet generation
		if args.Zjet is True:
			n_jets_accepted = 0
			if len(zs) < 1:
				# print('[w] no Z in the event?')
				continue
			# now we got the Z
			psjZ = fj.sorted_by_pt(zs)[0]
			phepmc = event_hepmc.particles[zs[0].user_index()]
			if args.debug:
				print('[i] Zs')
				for _z in fj.sorted_by_pt(zs):
					_phepmc = event_hepmc.particles[_z.user_index()]
					print(' - {} Z(s) found {} pt={}, m={} uidx={} hepmcV={} hepmcPID={}'.format(event_number, len(zs), _z.perp(), _z.m(), _z.user_index(), _phepmc, _phepmc.pid))
			# Z should be close to either of two hardest jets (in R space).
			ZjetDR_min = 1e6
			Zjet = None
			for _j in jets_selected:
				ZjetDR = psjZ.delta_R(_j)
				#print_jet(_j)
				#print(' - deltaR Z-jet=', ZjetDR)
				if ZjetDR < ZjetDR_min:
					ZjetDR_min = ZjetDR
					Zjet = _j
			if Zjet is None:
				print('[w] unable to match Z and and any jet within R < 1.', ZjetDR)
				continue
			if ZjetDR_min > float(args.jet_R) / 2.:
				continue
			if args.debug:
				print('[selected]')
				print_jet(Zjet)
				print(' - deltaR Z-jet=', ZjetDR_min)
			jets_accepted.append([Zjet, psjZ])
			n_jets_accepted = 1

			if n_jets_accepted < 1:
				continue

			if len(jets_accepted) >= nev_nperfile:
				fout.write_file_Zjets(jets_accepted, pythia=None)
				jets_accepted.clear()

			n_total_jets_accepted += n_jets_accepted
			pbar.update(n_jets_accepted)
			continue

		break
 
	if len(jets_accepted) > 0:
		if args.Zjet is False:
			fout.write_file_partons(jets_accepted, pythia=None)
			jets_accepted.clear()
		else:
			fout.write_file_Zjets(jets_accepted, pythia=None)
			jets_accepted.clear()

	pbar.close()
	pbar_ev.close()


def count_threads_alive(threads):
	_count = len([thr for thr in threads if thr.is_alive()])
	return _count


def launch_with_threads(args, nthreads=multiprocessing.cpu_count() * 2):
	print('[i] opening', args.input, 'as list of files...')
	with open(args.input, 'r') as f:
		file_list = [fn.strip('\n') for fn in f.readlines()]
	base_output = args.output
	threads = list()
	for i, fn in enumerate(tqdm.tqdm(file_list, desc='files')):
		if base_output:
			args.output = base_output + '_{}'.format(i)
		_t = threading.Thread(target=analyze_file, args=(fn, args,))
		threads.append(_t)
		_t.start()
		while count_threads_alive(threads) >= nthreads:
			_ = [thr.join(0.1) for thr in threads if thr.is_alive()]


def main():
	parser = argparse.ArgumentParser(description='pythia8 in python', prog=os.path.basename(__file__))
	parser.add_argument('--jet-R', help='specify jet R for flavor tagging', default=0.8, type=float)
	parser.add_argument('-c', '--nperfile', help='nevents per output file', default=1000, type=int)
	parser.add_argument('--jet-ptmin', help='minimum pT cut on jets', default=500., type=float)
	parser.add_argument('--jet-ptmax', help='maximum pT cut on jets', default=550., type=float)
	parser.add_argument('--Z-ptmin', help='minimum pT cut on Z', default=50., type=float)
	parser.add_argument('--Z-ptmax', help='maximum pT cut on Z', default=550., type=float)
	parser.add_argument('--Zjet', help='force Zjet generation - will read pythia_gen_qorg_Zjet_master.cmnd from current directory', default=False, action='store_true')
	parser.add_argument('--mDT', help='mDT clean up with z_cut=0.04 (default) on the Zjets', default=False, action='store_true')
	parser.add_argument('--mDTzcut', help='mDT clean up with z_cut=0.04 (default) on the Zjets', default=0.04, type=float)
	# parser.add_argument('--ZjetR', help='specify the radius of the anti-kT jet for Z-match - default is R=1.0', default=1.0, type=float)
	parser.add_argument('-i', '--input', help='hepmc file input to analyze', default=None, type=str)
	parser.add_argument('--nev', help='number of events to process', default=1000, type=int)
	parser.add_argument('-l', '--list', help='treat input file as a file containing list of files to process', action='store_true', default=False)
	parser.add_argument('-g', '--debug', help='print some extras', default=False, action='store_true')
	parser.add_argument('-o', '--output', help='overwrite default naming of the output', default='', type=str)
	parser.add_argument('-t', '--threads', help='mutlithread', default=1, type=int)
	parser.add_argument('--slurm', help='get run number from slurm batch job id', action='store_true', default=False)

	args = parser.parse_args()

	# print(args)

	args.mDT = args.Zjet
	if args.Zjet is False:
		del args.mDT
		del args.mDTzcut

	fj.ClusterSequence.print_banner()

	if args.list:
		launch_with_threads(args, args.threads)
	else:
		fname = args.input
		analyze_file(fname, args)


if __name__ == "__main__":
	main()
