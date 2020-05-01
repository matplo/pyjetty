from pyjetty.mputils import MPBase, pinfo, perror, pwarning
import os
import psutil
import sys
import array
import time

import ROOT
ROOT.gROOT.SetBatch(True)

class MemTrace(MPBase):
	_instance = None

	def __init__(self):
		raise RuntimeError('Call instance() instead')

	@classmethod
	def instance(cls, **kwargs):
		if cls._instance is None:
			pinfo('Creating new MemTrace instance')
			cls._instance = cls.__new__(cls)
			super(MemTrace, cls._instance).__init__(**kwargs)
			cls._instance._process = psutil.Process(os.getpid())
			cls._instance.event_number = 0
			cls._instance.event_tree_name = 'mt'
			cls._instance.process = psutil.Process(os.getpid())
			cls._instance.trees = {}
			cls._instance.fout = None
			cls._instance.output_name='memtrace.root'
			cls._instance.toffset = time.time()
			cls._partial_write = False
		cls._instance.configure_from_args(**kwargs)
		return cls._instance

	def reset(self, **kwargs):
		self.configure_from_args(**kwargs)
		self.reset_output()

	def reset_output(self):
		cwd = ROOT.gDirectory.CurrentDirectory()
		ROOT.gDirectory.cd('/')
		if self.fout is None:
			self.fout = ROOT.TFile(self.output_name, 'recreate')
		else:
			self.fout.Close()
			self.fout = ROOT.TFile(self.output_name, 'recreate')
		pinfo('MemTrace.reset_output', self.output_name, 'path:', ROOT.gDirectory.GetPath())
		cwd.cd()

	def snapshot(self, label='mem'):
		cwd = ROOT.gDirectory.CurrentDirectory().GetPath()
		new_key = False
		try:
			self.fout.cd()
		except:
			self.reset_output()
		try: 
			self.trees[label]
		except KeyError:
			new_key = True
		if new_key:
			self.fout.cd()
			self.trees[label] = ROOT.TNtuple(label, label, 'rss:vms:t')
			pinfo('MemTrace - new tuple {}'.format(label), file=sys.stderr)
		rss = self.process.memory_info().rss
		vms = self.process.memory_info().vms
		# n = self.trees[label].GetEntries()
		ts = time.time() - self.toffset
		self.trees[label].Fill(rss, vms, ts)
		if self._partial_write:
			self._write_()
		ROOT.gDirectory.cd(cwd)
		# cwd.cd()
		
	def write_as_graph(self, tree):
		data_x = []
		data_yrss = []
		data_yvms = []
		to_Gb = 1. / (1024.*1024*1024.)
		last_x   = -1
		last_rss = 0
		last_vms = 0
		for i,e in enumerate(tree):
			# if i == 0:
			# 	data_x.append(e.t - 0.1)
			# 	data_yrss.append(0)
			# 	data_yvms.append(0)				
			data_x.append(e.t)
			data_yrss.append(e.rss * to_Gb)
			data_yvms.append(e.vms * to_Gb)
			last_rss = e.rss
			last_vms = e.vms
			last_x = e.t
		# data_x.append(e.t + 0.1)
		# data_yrss.append(0)
		# data_yvms.append(0)				

		f_x = array.array('f', data_x)
		f_rss = array.array('f', data_yrss)
		f_vms = array.array('f', data_yvms)

		gr_rss = ROOT.TGraph(len(data_x), f_x, f_rss)
		gr_rss.SetName('gr_rss_{}'.format(tree.GetName()))
		gr_rss.SetTitle('rss {}'.format(tree.GetName()))
		gr_rss.Write()

		gr_vms = ROOT.TGraph(len(data_x), f_x, f_vms)
		gr_vms.SetName('gr_vms_{}'.format(tree.GetName()))
		gr_vms.SetTitle('vms {}'.format(tree.GetName()))
		gr_vms.Write()

	def _write_(self):
		self.fout.cd()
		for t in self.trees:
			self.trees[t].Write()
		self.fout.Write()
		self.fout.Purge()

	def write(self, write_graphs = True):
		self.fout.cd()
		for t in self.trees:
			self.trees[t].Write()
			if write_graphs:
				self.write_as_graph(self.trees[t])
		self.fout.Write()			
		self.fout.Close()
		pinfo('MemTrace.write - write & close {}'.format(self.output_name))
		# except:
		# 	perror('MemTrace.write - unable to write & close {}'.format(self.output_name))

	def __del__(self):
		try:
			self.fout.Write()			
			self.fout.Close()
		except:
			None

		self._instance = None
